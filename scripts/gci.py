#!/usr/bin/env python3
"""Three-mesh Richardson / GCI procedure (Celik et al. 2008).

Convention: index 1 = fine, 2 = medium, 3 = coarse.
h is the representative cell size (V_fluid/N)^(1/3); r21 = h2/h1, r32 = h3/h2.
Observed order p solves the transcendental equation for unequal r by
fixed-point iteration; Fs = 1.25 (three-mesh study).

Usage:  gci.py <N1> <V1> <f1>  <N2> <V2> <f2>  <N3> <V3> <f3>
   or:  import gci; gci.gci_study([(h1,f1),(h2,f2),(h3,f3)])
"""
import json
import math
import sys


def gci_study(hf):
    (h1, f1), (h2, f2), (h3, f3) = hf   # fine, medium, coarse
    r21, r32 = h2 / h1, h3 / h2
    eps21, eps32 = f2 - f1, f3 - f2
    R = eps21 / eps32 if eps32 != 0 else float('inf')
    if 0 < R < 1:
        ctype = 'monotonic convergence'
    elif -1 < R < 0:
        ctype = 'oscillatory convergence'
    else:
        ctype = 'divergent (|R| >= 1)'
    if not abs(R) < 1:
        # differences grow with refinement: outside the asymptotic range,
        # observed order / Richardson extrapolation / GCI are undefined
        return {
            'r21': r21, 'r32': r32, 'eps21': eps21, 'eps32': eps32,
            'R': R, 'convergence': ctype,
            'p': None, 'f_ext21': None, 'f_ext32': None,
            'GCI_fine_pct': None, 'GCI_medium_pct': None,
            'medium_vs_ext_pct': None,
        }
    s = math.copysign(1.0, eps32 / eps21) if eps21 != 0 else 1.0
    p = 2.0
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(abs(math.log(abs(eps32 / eps21))) + q) / math.log(r21)
        if abs(p_new - p) < 1e-12:
            p = p_new
            break
        p = 0.5 * (p + p_new)          # damped iteration
    fext21 = (r21 ** p * f1 - f2) / (r21 ** p - 1.0)
    fext32 = (r32 ** p * f2 - f3) / (r32 ** p - 1.0)
    ea21 = abs((f1 - f2) / f1)
    ea32 = abs((f2 - f3) / f2)
    return {
        'r21': r21, 'r32': r32,
        'eps21': eps21, 'eps32': eps32,
        'R': R, 'convergence': ctype,
        'p': p,
        'f_ext21': fext21, 'f_ext32': fext32,
        'GCI_fine_pct': 100 * 1.25 * ea21 / (r21 ** p - 1.0),
        'GCI_medium_pct': 100 * 1.25 * ea32 / (r32 ** p - 1.0),
        'medium_vs_ext_pct': 100 * abs((f2 - fext21) / fext21),
    }


if __name__ == '__main__':
    v = [float(x) for x in sys.argv[1:10]]
    (N1, V1, f1), (N2, V2, f2), (N3, V3, f3) = v[0:3], v[3:6], v[6:9]
    hs = [(V / N) ** (1.0 / 3.0) for N, V in ((N1, V1), (N2, V2), (N3, V3))]
    out = gci_study(list(zip(hs, (f1, f2, f3))))
    out['h_mm'] = [1000 * h for h in hs]
    out['N'] = [N1, N2, N3]
    print(json.dumps(out, indent=1))
