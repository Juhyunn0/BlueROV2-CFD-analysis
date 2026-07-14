"""Re-run remesh keeping the MC mesh, decimate with watertightness repair."""
import json

import numpy as np
import trimesh
from scipy import ndimage
from skimage import measure

SRC = '/home/bdml/Desktop/OpenFoam12/rovCase/constant/geometry/BROV2.stl'
PITCH = 0.005

mesh = trimesh.load_mesh(SRC)
lo = mesh.bounds[0] - 3 * PITCH
hi = mesh.bounds[1] + 3 * PITCH
dims = np.ceil((hi - lo) / PITCH).astype(int) + 1
occ = np.zeros(dims, dtype=bool)


def mark(points):
    idx = np.floor((points - lo) / PITCH).astype(int)
    idx = np.clip(idx, 0, np.array(dims) - 1)
    occ[idx[:, 0], idx[:, 1], idx[:, 2]] = True


mark(mesh.vertices)
n_samp = int(mesh.area / PITCH ** 2 * 8)
for _ in range(4):
    pts, _ = trimesh.sample.sample_surface(mesh, n_samp // 4)
    mark(pts)

occ_closed = ndimage.binary_closing(occ, structure=np.ones((3, 3, 3)))
lbl, _ = ndimage.label(~occ_closed)
solid = ~(lbl == lbl[0, 0, 0])

verts, faces, _, _ = measure.marching_cubes(solid.astype(np.uint8), level=0.5)
wt = trimesh.Trimesh(vertices=verts * PITCH + lo, faces=faces)
trimesh.repair.fix_normals(wt)
if wt.volume < 0:
    wt.invert()
print(f'MC mesh: {len(wt.faces)} faces, watertight={wt.is_watertight}, '
      f'V={wt.volume:.5f}')

for target in (12000, 16000, 20000):
    dec = wt.simplify_quadric_decimation(face_count=target)
    dec.update_faces(dec.nondegenerate_faces())
    dec.merge_vertices()
    trimesh.repair.fill_holes(dec)
    trimesh.repair.fix_normals(dec)
    if dec.volume < 0:
        dec.invert()
    print(f'target {target}: {len(dec.faces)} faces, '
          f'watertight={dec.is_watertight}, V={dec.volume:.5f}')
    if dec.is_watertight:
        dec.export('/home/bdml/Desktop/OpenFoam12/bem/BROV2_watertight.stl')
        json.dump(dict(faces=int(len(dec.faces)),
                       watertight=True,
                       volume=float(dec.volume),
                       centroid=[float(c) for c in dec.center_mass],
                       volume_mc=float(wt.volume), pitch=PITCH),
                  open('/home/bdml/Desktop/OpenFoam12/bem/mesh_meta.json', 'w'),
                  indent=1)
        print('saved watertight decimated mesh; centroid =', dec.center_mass)
        break
else:
    raise SystemExit('decimation never watertight — need pymeshlab fallback')
