#!/usr/bin/env python3
"""Pack SPMM melting frames into one small JSON file for the melting-steps viewer.

usage: pack_melt_frames.py OUT.json FRAME [FRAME ...]
   e.g. pack_melt_frames.py melt.json run/domain5_f*.out

Works on a run that is still going: pack whatever frames exist so far, and an
unfinished last line in a frame that is being written is ignored.

Frame lines are "index x y z type [body]" or "x y z type"; type 1 is ice, -1 or
2 is liquid.  Header lines are skipped.  Connected bodies (26 neighbours, ice
and liquid together) are recomputed here, so the body column is not required.
The melted fraction comes from six digits in a name like domain5_f049477.out;
for other names it is the share of the largest frame that is no longer ice.
"""
import base64
import json
import os
import re
import sys
import zlib

import numpy as np
from scipy import ndimage


def read_frame(path):
    rows, width = [], None
    with open(path, "rb") as fh:
        data = fh.read()
    lines = data.split(b"\n")
    complete = data.endswith(b"\n")
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            vals = [int(p) for p in parts[:6]]
        except ValueError:
            continue
        last = i == len(lines) - 1
        if last and not complete and width is not None and len(parts) != width:
            continue                      # line still being written
        width = len(parts)
        rows.append(vals[:4] if len(parts) == 4 else vals[1:5])
    return np.array(rows, dtype=np.int64).reshape(-1, 4)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out, paths = sys.argv[1], sys.argv[2:]
    frames = []
    for path in sorted(paths):
        a = read_frame(path)
        if len(a) == 0:
            print("skipped (no dipole lines): %s" % path)
            continue
        frames.append((os.path.basename(path), a))
    if not frames:
        sys.exit("no frames with dipole lines")
    lo = np.min([a[:, :3].min(0) for _, a in frames], 0)
    hi = np.max([a[:, :3].max(0) for _, a in frames], 0)
    dims = (hi - lo + 1).astype(int)
    n0 = max(len(a) for _, a in frames)
    structure = np.ones((3, 3, 3), dtype=bool)
    packed = []
    for name, a in frames:
        p = a[:, :3] - lo
        liquid = (a[:, 3] == -1) | (a[:, 3] == 2)
        g = np.zeros(dims, np.uint8)
        g[p[:, 0], p[:, 1], p[:, 2]] = np.where(liquid, 2, 1)
        lab, nbody = ndimage.label(g > 0, structure=structure)
        sizes = np.bincount(lab.ravel())[1:]
        detached = (lab > 0) & (lab != sizes.argmax() + 1)
        g[detached] += 2                  # 3 detached ice, 4 detached liquid
        nice, nliq = int((~liquid).sum()), int(liquid.sum())
        ev = np.sqrt(np.maximum(np.sort(np.linalg.eigvalsh(np.cov((p - p.mean(0)).T))), 0))
        m = re.search(r"f(\d{6})\.out$", name)
        packed.append(dict(
            name=name,
            label=int(m.group(1)) / 1000.0 if m else 100.0 * (n0 - nice) / n0,
            ice=nice, liq=nliq, removed=n0 - nice - nliq,
            bodies=int(nbody), detached=int(detached.sum()),
            ax=[round(float(ev[0] / ev[2]), 3), round(float(ev[1] / ev[2]), 3)],
            ext=(p.max(0) - p.min(0) + 1).tolist(),
            z=base64.b64encode(zlib.compress(g.tobytes(), 9)).decode()))
    packed.sort(key=lambda f: f["label"])
    with open(out, "w") as fh:
        json.dump(dict(dims=dims.tolist(), lo=lo.tolist(), n0=int(n0), frames=packed), fh, separators=(",", ":"))
    print("%d frames, grid %s, %d dipoles in the largest frame -> %s (%.0f kB)" % (
        len(packed), "x".join(map(str, dims)), n0, out, os.path.getsize(out) / 1e3))


if __name__ == "__main__":
    main()
