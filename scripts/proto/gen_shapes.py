"""Generate test shapes: rows 'i j k m' (m=2 water, m=0 ice), one header line.
usage: gen_shapes.py blob|slab|film OUT [args]"""
import sys, numpy as np
kind, out = sys.argv[1], sys.argv[2]; a = [float(x) for x in sys.argv[3:]]
rows = []
if kind == "blob":            # water box nx ny nz
    nx, ny, nz = map(int, a[:3])
    rows = [(i, j, k, 2) for i in range(1, nx+1) for j in range(1, ny+1) for k in range(1, nz+1)]
elif kind == "slab":          # ice slab L x L x T with a water cube of edge c on top
    L, T, c = map(int, a[:3]); o = (L - c)//2
    rows = [(i, j, k, 0) for i in range(1, L+1) for j in range(1, L+1) for k in range(1, T+1)]
    rows += [(o+i, o+j, T+k, 2) for i in range(1, c+1) for j in range(1, c+1) for k in range(1, c+1)]
elif kind == "film":          # ice sphere radius r, water shell thickness t, optional water blob radius b on +x
    r, t, b = a[0], a[1], (a[2] if len(a) > 2 else 0.0)
    M = int(r + t + 2*b + 3); cb = r + t + 0.6*b
    for i in range(-M, M+1):
        for j in range(-M, M+1):
            for k in range(-M, M+1):
                d = (i*i + j*j + k*k) ** 0.5
                if d <= r: rows.append((i, j, k, 0))
                elif d <= r + t or (b > 0 and ((i-cb)**2 + j*j + k*k) ** 0.5 <= b): rows.append((i, j, k, 2))
with open(out, "w") as f:
    f.write("x y z t\n")
    for r_ in rows: f.write("%d %d %d %d\n" % r_)
print(kind, "N=%d ice=%d water=%d" % (len(rows), sum(1 for r_ in rows if r_[3] == 0), sum(1 for r_ in rows if r_[3] == 2)))
