#!/usr/bin/env python3
"""Morison fit for ROV forced-oscillation (oscillating-inflow) cases.

Model:  F_i(t) = c0 + alpha * dU/dt + beta * U|U|,   U(t) = U0 sin(2 pi t / T)
fitted by least squares on the last 3 of 5 periods (t >= 2T).

alpha = m_a + rho*V_solid  (Froude-Krylov of the oscillating ambient flow),
so m_a = alpha - RHO_V.  beta = 0.5*rho*Cd*A with A the bbox projected area
normal to the oscillation direction.

As-meshed bbox extents of BROV2v.stl (vehicle frame, verified 2026-07-14):
x (surge) 0.457200, y (sway) 0.574845, z (heave) 0.253852 m.

Usage: morison_fit.py <caseDir> <x|y|z> <T> <U0> [--json OUT] [--rhoV KG]
--rhoV overrides the Froude-Krylov mass (default 14.40 kg, the production
medium mesh); mesh-convergence fits must pass each mesh's own rho*V_solid.
"""
import json
import sys

import numpy as np

RHO = 1000.0
RHO_V = 14.40           # rho * V_solid of the as-meshed body [kg]
EXT = {'x': 0.457200, 'y': 0.574845, 'z': 0.253852}   # bbox extents [m]
COMP = {'x': 0, 'y': 1, 'z': 2}


def area(direction):
    others = [v for k, v in EXT.items() if k != direction]
    return others[0] * others[1]


def read_forces(case):
    """forces.dat: t ((px py pz) (vx vy vz)) ((...) (...)) -> t, F[3] total."""
    path = f'{case}/postProcessing/forces/0/forces.dat'
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            vals = [float(v) for v in
                    line.replace('(', ' ').replace(')', ' ').split()]
            rows.append(vals)
    a = np.array(rows)
    t = a[:, 0]
    F = a[:, 1:4] + a[:, 4:7]          # pressure + viscous
    return t, F


def fit(case, direction, T, U0, rho_v=RHO_V):
    t, Fall = read_forces(case)
    F = Fall[:, COMP[direction]]
    w = t >= 2.0 * T                   # last 3 of 5 periods
    t, F = t[w], F[w]
    om = 2.0 * np.pi / T
    U = U0 * np.sin(om * t)
    dU = U0 * om * np.cos(om * t)
    X = np.column_stack((np.ones_like(t), dU, U * np.abs(U)))
    coef, *_ = np.linalg.lstsq(X, F, rcond=None)
    c0, alpha, beta = coef
    res = F - X @ coef
    ss_res = float(np.sum(res ** 2))
    ss_tot = float(np.sum((F - F.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot
    resid = float(np.sqrt(np.mean(res ** 2)) / np.max(np.abs(F)))
    A = area(direction)
    out = {
        'case': case.rstrip('/').split('/')[-1],
        'direction': direction,
        'T': T,
        'U0': U0,
        'KC': U0 * T / EXT[direction],
        'alpha': float(alpha),
        'rhoV': rho_v,
        'ma': float(alpha) - rho_v,
        'beta': float(beta),
        'Cd': float(beta) / (0.5 * RHO * A),
        'A_ref': A,
        'c0': float(c0),
        'r2': r2,
        'resid': resid,
    }
    return out


if __name__ == '__main__':
    case, direction, T, U0 = sys.argv[1], sys.argv[2], float(sys.argv[3]), \
        float(sys.argv[4])
    rho_v = float(sys.argv[sys.argv.index('--rhoV') + 1]) \
        if '--rhoV' in sys.argv else RHO_V
    out = fit(case, direction, T, U0, rho_v)
    print(json.dumps(out, indent=1))
    if '--json' in sys.argv:
        path = sys.argv[sys.argv.index('--json') + 1]
        with open(path, 'w') as f:
            json.dump(out, f, indent=1)
