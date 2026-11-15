"""usage: analyze.py blob|slab|film INIT FINAL [slabT | r t]"""
import sys, numpy as np
from scipy import ndimage
kind, fi, ff = sys.argv[1:4]; extra = [float(x) for x in sys.argv[4:]]
def load(f):
    t = np.loadtxt(f, dtype=int); return t[:, 1:4], t[:, 4]
def axes(p):
    c = p.mean(0); ev = np.sqrt(np.sort(np.linalg.eigvalsh(np.cov((p - c).T)))); return ev[0]/ev[2], ev[1]/ev[2]
def ncomp(p):
    q = p - p.min(0) + 1; g = np.zeros(q.max(0) + 2, bool); g[q[:, 0], q[:, 1], q[:, 2]] = True
    lab, n = ndimage.label(g, structure=np.ones((3, 3, 3))); sz = np.bincount(lab.ravel())[1:]
    return n, sz.max()
p0, t0 = load(fi); p1, t1 = load(ff); w0 = p0[t0 == -1]; w1 = p1[t1 == -1]; ice = p1[t1 == 1]
if kind == "blob":
    n, big = ncomp(w1)
    print("axis ratios init %.3f %.3f -> final %.3f %.3f | centroid shift %s | components %d (largest %d of %d)" % (
        *axes(w0), *axes(w1), np.round(w1.mean(0) - w0.mean(0), 2), n, big, len(w1)))
elif kind == "slab":
    T = int(extra[0]); z = w1[:, 2] - T               # layer index above slab top (1 = first layer)
    h = z.max(); base = (z == 1).sum(); rb = np.sqrt(base/np.pi)
    top = ice[:, 2] == T; L = ice[:, 0].max() - ice[:, 0].min() + 1
    below = (z < 1).sum()
    theta = np.degrees(2*np.arctan(h/rb)) if base else float("nan")
    cxy = w1[z >= 1][:, :2].mean(0) - ice[:, :2].mean(0)
    print("height %d base radius %.1f -> contact angle %.0f deg | base covers %.0f%% of slab top | water below top %d | drop centre offset %s | max lateral reach %d of %d" % (
        h, rb, theta, 100*base/top.sum(), below, np.round(cxy, 1), max(np.ptp(w1[:, 0]), np.ptp(w1[:, 1])) + 1, L))
elif kind == "film":
    r, t = extra[:2]; c = ice.mean(0)
    g = {}
    for q in map(tuple, w1): g[q] = 1
    iceset = set(map(tuple, ice)); nb = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0)]
    surf = cov = 0
    for q in iceset:
        ns = [(q[0]+d[0], q[1]+d[1], q[2]+d[2]) for d in nb]
        if any(x not in iceset for x in ns):
            surf += 1; cov += any(x in g for x in ns)
    d = np.sqrt(((w1 - c)**2).sum(1)); n, big = ncomp(w1)
    print("ice surface covered %.1f%% | water in film zone (r<=%.1f) %.1f%% | max water radius %.1f | water centroid offset %.2f | water bodies %d" % (
        100*cov/surf, r+t+0.5, 100*(d <= r+t+0.5).mean(), d.max(), np.linalg.norm(w1.mean(0) - c), n))
