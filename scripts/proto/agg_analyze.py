"""usage: agg_analyze.py LABEL GLOB [targets...]  -- per-frame water distribution metrics"""
import sys, glob, numpy as np
from scipy import ndimage
label, pat = sys.argv[1], sys.argv[2]; want = [float(x) for x in sys.argv[3:]] or [0.05, 0.10, 0.20, 0.50]
K = np.ones((3, 3, 3)); fs = sorted(glob.glob(pat))
def grids(f):
    t = np.loadtxt(f, dtype=int); p = t[:, 1:4]; ty = t[:, 4]
    q = p - p.min(0) + 2; sh = q.max(0) + 3
    ice = np.zeros(sh, bool); wat = np.zeros(sh, bool)
    ice[tuple(q[ty == 1].T)] = True; wat[tuple(q[ty == -1].T)] = True
    return p, ty, ice, wat
frames = []
for f in fs:
    with open(f) as fh: ty = np.fromiter((int(l.split()[4]) for l in fh), int)
    frames.append((f, (ty == -1).mean()))
print("%s: %d frames, actual liquid fractions: %s" % (label, len(fs), " ".join("%.3f" % a for _, a in frames)))
p0 = None
for w in want + [1.0]:
    f, a = min(frames, key=lambda fa: abs(fa[1] - w))
    if abs(a - w) > 0.03: print("  target %.2f: nearest frame is %.3f (skipped)" % (w, a)); continue
    p, ty, ice, wat = grids(f)
    if w == 1.0:
        c = p.mean(0); ev = np.sqrt(np.sort(np.linalg.eigvalsh(np.cov((p - c).T))))
        i0 = np.loadtxt(fs[0], dtype=int)[:, 1:4].mean(0); lab, n = ndimage.label(wat, structure=K); sz = np.bincount(lab.ravel())[1:]
        print("  fully melted: axis ratios %.3f %.3f | centroid shift %s | water bodies %d (largest holds %.1f%%)" % (
            ev[0]/ev[2], ev[1]/ev[2], np.round(c - i0, 1), n, 100*sz.max()/sz.sum())); continue
    icesurf = ice & ~ndimage.binary_erosion(ice, structure=K)
    nearw = ndimage.binary_dilation(wat, structure=K); neari = ndimage.binary_dilation(ice, structure=K)
    nwn = ndimage.convolve(wat.astype(np.int16), K.astype(np.int16), mode="constant") - wat
    lab, n = ndimage.label(wat, structure=K); sz = np.bincount(lab.ravel())[1:]
    print("  liq %.3f: ice surface wetted %4.1f%% | water touching ice %4.1f%% | mean water neighbours %4.1f | water bodies %5d, largest %4.1f%% of water" % (
        a, 100*(icesurf & nearw).sum()/icesurf.sum(), 100*(wat & neari).sum()/wat.sum(), nwn[wat].mean(), n, 100*sz.max()/sz.sum()))
