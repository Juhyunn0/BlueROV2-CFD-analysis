#!/usr/bin/env python3
"""KC -> 0 extrapolation of added mass per direction (Step 2b convention).

Fits BOTH  m_a = a + b*KC  and  m_a = a + b*KC^(2/3)  per direction.
Official coefficient = midpoint of the two intercepts +/- half their spread
(extrapolation-model uncertainty), reported alongside each fit's statistical
intercept SE. Sanity brackets (audit if midpoint outside): surge [3, 13] kg,
sway [5, 18] kg (BEM sealed-envelope refs 9.75 / 14.14 kg).

Inputs: heave_kc_results.json (Step 2a) + rovOsc_fit_*.json (Step 2b).
Outputs: step2b_results.json, step2b_ma_vs_KC.png.
"""
import glob
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = '/home/bdml/Desktop/OpenFoam12'
BRACKETS = {'x': (3.0, 13.0), 'y': (5.0, 18.0), 'z': None}
NAMES = {'x': 'surge', 'y': 'sway', 'z': 'heave'}
COLORS = {'x': '#2a78d6', 'y': '#1baf7a', 'z': '#eda100'}
BEM_REF = {'x': 9.75, 'y': 14.14, 'z': 35.03}


def load_points():
    pts = {'x': [], 'y': [], 'z': []}
    for r in json.load(open(f'{ROOT}/heave_kc_results.json')):
        pts['z'].append(r)
    for f in sorted(glob.glob(f'{ROOT}/rovOsc_fit_*.json')):
        r = json.load(open(f))
        if 'KC' not in r:      # Step-2a-era raw fit files (different schema)
            continue
        d = r.get('direction', 'z')
        # avoid double-counting the 2a heave fits already in the results file
        if d == 'z' and any(abs(p['KC'] - r['KC']) < 1e-6 for p in pts['z']):
            continue
        pts[d].append(r)
    for d in pts:
        pts[d].sort(key=lambda r: r['KC'])
    return pts


def fit_model(kc, ma, transform):
    x = transform(kc)
    X = np.column_stack((np.ones_like(x), x))
    coef, *_ = np.linalg.lstsq(X, ma, rcond=None)
    res = ma - X @ coef
    dof = max(len(kc) - 2, 1)
    sigma2 = float(res @ res) / dof
    cov = sigma2 * np.linalg.inv(X.T @ X)
    return {'a': float(coef[0]), 'b': float(coef[1]),
            'se_a': float(np.sqrt(cov[0, 0])),
            'max_resid': float(np.max(np.abs(res)))}


def main():
    pts = load_points()
    out = {}
    fig, ax = plt.subplots(figsize=(9, 6.5))
    for d in ('x', 'y', 'z'):
        rows = pts[d]
        if len(rows) < 2:
            continue
        kc = np.array([r['KC'] for r in rows])
        ma = np.array([r['ma'] for r in rows])
        lin = fit_model(kc, ma, lambda v: v)
        gra = fit_model(kc, ma, lambda v: v ** (2.0 / 3.0))
        mid = 0.5 * (lin['a'] + gra['a'])
        spread = 0.5 * abs(lin['a'] - gra['a'])
        bracket = BRACKETS[d]
        status = 'ok'
        if bracket and not (bracket[0] <= mid <= bracket[1]):
            status = 'HARD_STOP_AUDIT'
        out[NAMES[d]] = {
            'points': [{'case': r['case'], 'T': r['T'], 'U0': r['U0'],
                        'KC': r['KC'], 'alpha': r['alpha'], 'rhoV': r['rhoV'],
                        'ma': r['ma'], 'Cd': r['Cd'], 'r2': r['r2'],
                        'resid': r['resid']} for r in rows],
            'fit_linear': lin, 'fit_graham23': gra,
            'official_ma_kg': mid,
            'model_uncertainty_kg': spread,
            'stat_se_kg': {'linear': lin['se_a'], 'graham23': gra['se_a']},
            'bracket': bracket, 'bracket_status': status,
            'bem_reference_kg': BEM_REF[d],
        }
        c = COLORS[d]
        ax.plot(kc, ma, 'o', color=c, ms=8, zorder=5)
        kk = np.linspace(0, kc.max() * 1.05, 100)
        ax.plot(kk, lin['a'] + lin['b'] * kk, '-', color=c, lw=2, alpha=0.9)
        ax.plot(kk, gra['a'] + gra['b'] * kk ** (2 / 3), '--', color=c,
                lw=2, alpha=0.9)
        ax.plot([0, 0], [lin['a'], gra['a']], 'o', mfc='white', mec=c,
                ms=8, zorder=6)
        ax.annotate(f"{NAMES[d]}  {mid:.1f} ± {spread:.1f} kg",
                    xy=(kc.max(), ma.max()),
                    xytext=(kc.max() + 0.03, ma.max() + 0.6),
                    color='#333333', fontsize=11)
    ax.set_xlabel('KC = U0·T / L  (as-meshed L per direction)')
    ax.set_ylabel('added mass  m_a  [kg]')
    ax.set_title('BROV2 Heavy forced-oscillation CFD: m_a(KC), '
                 'solid = linear fit, dashed = KC$^{2/3}$ fit')
    ax.set_xlim(-0.05, 1.45)
    ax.set_ylim(0, 40)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(f'{ROOT}/step2b_ma_vs_KC.png', dpi=170)
    with open(f'{ROOT}/step2b_results.json', 'w') as f:
        json.dump(out, f, indent=1)
    for k, v in out.items():
        print(f"{k:6s} official m_a = {v['official_ma_kg']:.2f} "
              f"± {v['model_uncertainty_kg']:.2f} kg (model) "
              f"[lin {v['fit_linear']['a']:.2f} ± {v['stat_se_kg']['linear']:.2f}, "
              f"graham {v['fit_graham23']['a']:.2f} ± {v['stat_se_kg']['graham23']:.2f}] "
              f"-> {v['bracket_status']}")


if __name__ == '__main__':
    main()
