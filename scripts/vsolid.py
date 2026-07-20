#!/usr/bin/env python3
"""Solid volume from a closed boundary-patch VTK: V = |(1/3) sum Cf . Sf|.

Reads the binary legacy VTK polydata that `foamToVTK -constant -noInternal
-excludePatches '(xmin xmax ymin ymax zmin zmax)'` writes for the BROV2
patch (VTK/BROV2/BROV2_*.vtk). Face normals of a boundary patch all point
out of the fluid (into the body), so the divergence-theorem sum gives the
enclosed body volume up to sign.

Usage: vsolid.py <patch.vtk>
"""
import struct
import sys

import numpy as np


def read_binary_polydata(path):
    with open(path, 'rb') as f:
        data = f.read()
    # locate POINTS
    i = data.index(b'POINTS')
    j = data.index(b'\n', i)
    n = int(data[i:j].split()[1])
    pts = np.frombuffer(data, dtype='>f4', count=3 * n,
                        offset=j + 1).reshape(-1, 3).astype(np.float64)
    # locate POLYGONS
    i = data.index(b'POLYGONS')
    j = data.index(b'\n', i)
    parts = data[i:j].split()
    total = int(parts[2])
    raw = np.frombuffer(data, dtype='>i4', count=total, offset=j + 1)
    polys = []
    k = 0
    while k < total:
        c = raw[k]
        polys.append(raw[k + 1:k + 1 + c])
        k += c + 1
    return pts, polys


def volume(path):
    pts, polys = read_binary_polydata(path)
    vol6 = 0.0
    for p in polys:
        # fan-triangulate each polygon; signed tetra volumes about origin
        v0 = pts[p[0]]
        for a, b in zip(p[1:-1], p[2:]):
            vol6 += np.dot(v0, np.cross(pts[a], pts[b]))
    return abs(vol6) / 6.0, len(polys)


if __name__ == '__main__':
    v, nf = volume(sys.argv[1])
    print(f'V_solid = {v:.6f} m^3   ({nf} patch faces)')
