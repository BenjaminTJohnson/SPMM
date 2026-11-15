#!/usr/bin/env python3
"""Run spmm in a run folder and pack its frames for the melting-steps viewer.

usage: run_packed.py [--live] [--dir RUN_DIR] [--exe PATH] [--out NAME.json] [--clean]

The run folder must hold shape.dat (and, optionally, spmm.nml).  The program's output goes
to run.log there and is echoed to the terminal.  When the run ends, all frames are packed
into one JSON file (default packed.json in the run folder) with pack_melt_frames.py.

  --live    also repack each time a new frame has been written, so the JSON can be opened
            (and reopened) in the viewer while the run is going.  Each repack rebuilds the
            whole file, which takes up to a minute for a large particle.
  --exe     program to run; default is build/src/spmm next to this repository's scripts.
  --clean   delete frame files already in the run folder first.  Without it the script
            stops if any are there, so an old run is never mixed into the new one.

Needs numpy and scipy (for the packer).
"""
import argparse
import glob
import os
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PACKER = os.path.join(HERE, "pack_melt_frames.py")
FRAMES = "*_f[0-9][0-9][0-9][0-9][0-9][0-9].out"


def frames(run_dir):
    return sorted(glob.glob(os.path.join(run_dir, FRAMES)))


def pack(files, out):
    """Pack to a temporary file, then rename, so the viewer never reads a half-written JSON."""
    if not files:
        return False
    tmp = out + ".tmp"
    r = subprocess.run([sys.executable, PACKER, tmp] + files, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write("packing failed:\n" + r.stdout + r.stderr)
        if os.path.exists(tmp):
            os.remove(tmp)
        return False
    os.replace(tmp, out)
    print("[run_packed] packed %d frames into %s" % (len(files), out), flush=True)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--live", action="store_true", help="repack whenever a new frame has been written")
    ap.add_argument("--dir", default=".", help="run folder (default: current folder)")
    ap.add_argument("--exe", default=os.path.join(HERE, "..", "build", "src", "spmm"))
    ap.add_argument("--out", default="packed.json", help="packed file name, inside the run folder unless a path is given")
    ap.add_argument("--clean", action="store_true", help="delete frame files already in the run folder first")
    a = ap.parse_args()

    run_dir = os.path.abspath(a.dir)
    exe = os.path.abspath(a.exe)
    out = a.out if os.path.isabs(a.out) else os.path.join(run_dir, a.out)
    if not os.path.isfile(os.path.join(run_dir, "shape.dat")):
        sys.exit("no shape.dat in " + run_dir)
    if not (os.path.isfile(exe) and os.access(exe, os.X_OK)):
        sys.exit("program not found or not executable: %s (build it, or pass --exe)" % exe)
    old = frames(run_dir)
    if old:
        if not a.clean:
            sys.exit("%d frame files are already in %s; move them away or pass --clean to delete them" % (len(old), run_dir))
        for f in old:
            os.remove(f)

    log = open(os.path.join(run_dir, "run.log"), "w")
    t0 = time.time()
    proc = subprocess.Popen([exe], cwd=run_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    def echo():
        for line in proc.stdout:
            log.write(line)
            log.flush()
            sys.stdout.write(line)
            sys.stdout.flush()

    th = threading.Thread(target=echo, daemon=True)
    th.start()

    packed = 0
    while proc.poll() is None:
        time.sleep(2.0)
        if not a.live:
            continue
        # every frame except the newest is complete; the newest is too once it has stopped growing
        fs = frames(run_dir)
        if fs and time.time() - os.path.getmtime(fs[-1]) < 3.0:
            fs = fs[:-1]
        if len(fs) > packed and pack(fs, out):
            packed = len(fs)
    th.join()
    log.close()
    print("[run_packed] %s finished with exit code %d after %.0f s" % (os.path.basename(exe), proc.returncode, time.time() - t0), flush=True)
    ok = pack(frames(run_dir), out)
    sys.exit(proc.returncode if proc.returncode else (0 if ok else 1))


if __name__ == "__main__":
    main()
