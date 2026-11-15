import sys, numpy as np
from scipy import ndimage
for f in sys.argv[1:]:
    t = np.loadtxt(f, dtype=int); p = t[:, 1:4]
    q = p - p.min(0) + 1; g = np.zeros(q.max(0) + 2, bool); g[tuple(q.T)] = True
    lab, n = ndimage.label(g, structure=np.ones((3, 3, 3))); sz = np.bincount(lab.ravel())[1:]
    keep = lab[tuple(q.T)] == (sz.argmax() + 1); allp = p.astype(float); p = allp[keep]
    def ax(x):
        ev = np.sqrt(np.sort(np.linalg.eigvalsh(np.cov((x - x.mean(0)).T)))); return ev[0]/ev[2], ev[1]/ev[2]
    c = p.mean(0); r = np.linalg.norm(p - c, axis=1); aeq = (3*len(p)/(4*np.pi))**(1/3)
    print("%-34s main body %.3f %.3f 1 (all dipoles %.3f %.3f) | rmax/aeq %.3f | bodies %d, detached %d | centre %s" % (
        f[-34:], *ax(p), *ax(allp), r.max()/aeq, n, len(allp) - len(p), np.round(c, 1)))
