#!/usr/bin/env python3
"""Write DDSCAT target files (the FROM_FILE shape.dat format) from SPMM frames.

usage: frames_to_ddscat.py [--out DIR] [--spacing MICRONS] [--main-body] FRAME [FRAME ...]
   e.g. frames_to_ddscat.py --spacing 20 run/agg_f*.out

For each frame this writes DIR/<frame name>/shape.dat (default DIR is "ddscat" beside the
first frame), so each frame has its own DDSCAT run folder, and one table, DIR/frames.txt,
listing for every frame the melted fraction, the dipole counts and the effective radius
to put in ddscat.par.

Composition 1 is ice and composition 2 is liquid water in every file, so ddscat.par needs
NCOMP = 2 with the ice refractive-index file first and the water file second, including
for frames that are all ice or all water.

The target axes are a1 = x and a2 = y of the SPMM lattice.  The target origin is put at the
centre of the occupied dipoles (line 6 of the file), which does not affect cross sections.

DDSCAT takes the dipole spacing d from the effective radius, aeff = (3 N / 4 pi)^(1/3) d.
The table gives aeff/d; with --spacing it also gives aeff in micrometres.

  --main-body   keep only the largest connected body (needs the body column, column 6)

Frame lines are "index x y z type [body]"; type 1 is ice, -1 or 2 is liquid.  Needs numpy.
"""
import argparse
import math
import os
import re
import sys

import numpy as np


def read_frame(path):
    rows = []
    with open(path, "rb") as fh:
        data = fh.read()
    lines = data.split(b"\n")
    if not data.endswith(b"\n"):
        lines = lines[:-1]            # an unfinished last line in a frame that is still being written
    for line in lines:
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            vals = [int(p) for p in parts[:6]]
        except ValueError:
            continue                  # header or other text
        rows.append(vals + [1] * (6 - len(vals)))
    if not rows:
        raise ValueError("no dipole lines")
    return np.array(rows, dtype=np.int64)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("frames", nargs="+")
    ap.add_argument("--out", help='output folder (default: "ddscat" beside the first frame)')
    ap.add_argument("--spacing", type=float, help="dipole spacing in micrometres, to report aeff in micrometres")
    ap.add_argument("--main-body", action="store_true", help="keep only the largest connected body")
    a = ap.parse_args()

    out = a.out or os.path.join(os.path.dirname(os.path.abspath(a.frames[0])), "ddscat")
    os.makedirs(out, exist_ok=True)
    table = ["# SPMM frames as DDSCAT targets.  Composition 1 = ice, 2 = liquid water (NCOMP = 2).",
             "# aeff = aeff/d times the dipole spacing d" +
             (";  spacing d = %g micrometres" % a.spacing if a.spacing else ""),
             "# %-28s %9s %9s %9s %9s %10s%s" % ("folder", "melted_%", "NAT", "ice", "liquid", "aeff/d",
                                                 "  aeff_um" if a.spacing else "")]
    for path in a.frames:
        try:
            t = read_frame(path)
        except (OSError, ValueError) as e:
            sys.exit("%s: %s" % (path, e))
        if a.main_body:
            t = t[t[:, 5] == 1]
        liquid = (t[:, 4] == -1) | (t[:, 4] == 2)
        comp = np.where(liquid, 2, 1)
        xyz = t[:, 1:4]
        n = len(t)
        if n == 0:
            sys.exit("%s: no dipoles left" % path)
        x0 = -xyz.mean(axis=0)
        stem = os.path.splitext(os.path.basename(path))[0]
        m = re.search(r"f(\d{6})$", stem)
        melted = int(m.group(1)) / 1000.0 if m else float("nan")
        aeff_d = (3.0 * n / (4.0 * math.pi)) ** (1.0 / 3.0)
        d = os.path.join(out, stem)
        os.makedirs(d, exist_ok=True)
        desc = "SPMM %s: %.3f%% melted, %d ice, %d water" % (stem, melted, int((~liquid).sum()), int(liquid.sum()))
        with open(os.path.join(d, "shape.dat"), "w") as fh:
            fh.write(" >%s\n" % desc[:66])
            fh.write("%10d = NAT\n" % n)
            fh.write("  1.000000  0.000000  0.000000 = A_1 vector\n")
            fh.write("  0.000000  1.000000  0.000000 = A_2 vector\n")
            fh.write("  1.000000  1.000000  1.000000 = lattice spacings (d_x,d_y,d_z)/d\n")
            fh.write("%12.5f %12.5f %12.5f = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d for dipole 0 0 0\n" % tuple(x0))
            fh.write("     JA  IX  IY  IZ ICOMP(x,y,z)\n")
            body = np.column_stack([np.arange(1, n + 1), xyz, comp, comp, comp])
            np.savetxt(fh, body, fmt="%7d %5d %5d %5d %d %d %d")
        table.append("%-30s %9.3f %9d %9d %9d %10.4f%s" % (stem, melted, n, int((~liquid).sum()), int(liquid.sum()), aeff_d,
                                                          "  %7.2f" % (aeff_d * a.spacing) if a.spacing else ""))
    with open(os.path.join(out, "frames.txt"), "w") as fh:
        fh.write("\n".join(table) + "\n")
    print("%d shape files written under %s; see frames.txt there" % (len(a.frames), out))


if __name__ == "__main__":
    main()
