"""Final BEM mesh: decimate, remove non-manifold pinches, verify closedness."""
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
occ = ndimage.binary_closing(occ, structure=np.ones((3, 3, 3)))
lbl, _ = ndimage.label(~occ)
solid = ~(lbl == lbl[0, 0, 0])
verts, faces, _, _ = measure.marching_cubes(solid.astype(np.uint8), level=0.5)
wt = trimesh.Trimesh(vertices=verts * PITCH + lo, faces=faces)
trimesh.repair.fix_normals(wt)
if wt.volume < 0:
    wt.invert()
trimesh.smoothing.filter_taubin(wt, lamb=0.5, nu=-0.53, iterations=10)


def edge_multiplicity(m):
    unique, counts = np.unique(m.edges_sorted, axis=0, return_counts=True)
    return unique, counts


dec = wt.simplify_quadric_decimation(face_count=14000)
dec.update_faces(dec.nondegenerate_faces())
dec.merge_vertices()

for attempt in range(6):
    unique, counts = edge_multiplicity(dec)
    n_open = int((counts == 1).sum())
    n_pinch = int((counts > 2).sum())
    print(f'attempt {attempt}: faces={len(dec.faces)} open={n_open} '
          f'pinch={n_pinch} watertight={dec.is_watertight}')
    if n_open == 0 and n_pinch == 0:
        break
    # drop faces touching pinched edges, then refill the resulting holes
    bad_edges = set(map(tuple, unique[counts > 2]))
    mask = np.ones(len(dec.faces), dtype=bool)
    es = np.sort(dec.faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
    es = es.reshape(len(dec.faces), 3, 2)
    for i, ftris in enumerate(es):
        if any(tuple(e) in bad_edges for e in ftris):
            mask[i] = False
    dec.update_faces(mask)
    trimesh.repair.fill_holes(dec)
    dec.merge_vertices()
    dec.update_faces(dec.nondegenerate_faces())

trimesh.repair.fix_normals(dec)
if dec.volume < 0:
    dec.invert()
unique, counts = edge_multiplicity(dec)
n_open = int((counts == 1).sum())
n_pinch = int((counts > 2).sum())
ok = dec.is_watertight
print(f'FINAL: faces={len(dec.faces)}, watertight={ok}, open={n_open}, '
      f'pinch={n_pinch}, winding_consistent={dec.is_winding_consistent}, '
      f'V={dec.volume:.5f}, centroid={dec.center_mass}')
dec.export('/home/bdml/Desktop/OpenFoam12/bem/BROV2_watertight.stl')
json.dump(dict(faces=int(len(dec.faces)), watertight=bool(ok),
               open_edges=n_open, pinch_edges=n_pinch,
               volume=float(dec.volume),
               centroid=[float(c) for c in dec.center_mass],
               volume_mc=float(wt.volume), pitch=PITCH),
          open('/home/bdml/Desktop/OpenFoam12/bem/mesh_meta.json', 'w'),
          indent=1)
print('saved')
