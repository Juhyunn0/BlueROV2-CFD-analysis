"""Taubin smooth the MC mesh, decimate, and stitch any residual small holes."""
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
print(f'MC: {len(wt.faces)} faces, watertight={wt.is_watertight}, V={wt.volume:.5f}')

trimesh.smoothing.filter_taubin(wt, lamb=0.5, nu=-0.53, iterations=10)
print(f'after taubin: watertight={wt.is_watertight}, V={wt.volume:.5f}')


def boundary_stats(m):
    unique, counts = np.unique(m.edges_sorted, axis=0, return_counts=True)
    open_edges = unique[counts == 1]
    return len(open_edges)


best = None
for target in (12000, 16000):
    dec = wt.simplify_quadric_decimation(face_count=target)
    dec.update_faces(dec.nondegenerate_faces())
    dec.merge_vertices(merge_tex=True, merge_norm=True)
    trimesh.repair.fill_holes(dec)
    trimesh.repair.fix_normals(dec)
    if dec.volume < 0:
        dec.invert()
    nb = boundary_stats(dec)
    print(f'target {target}: {len(dec.faces)} faces, watertight={dec.is_watertight}, '
          f'open edges={nb}, V={dec.volume:.5f}')
    if dec.is_watertight:
        best = dec
        break
    if nb and nb < 60:
        # stitch remaining small boundary loops by fan triangulation
        for _ in range(4):
            unique, counts = np.unique(dec.edges_sorted, axis=0,
                                       return_counts=True)
            open_e = unique[counts == 1]
            if not len(open_e):
                break
            g = trimesh.graph.nx.Graph()
            g.add_edges_from(open_e)
            new_faces = []
            for loop in trimesh.graph.nx.cycle_basis(g):
                if len(loop) >= 3:
                    for i in range(1, len(loop) - 1):
                        new_faces.append([loop[0], loop[i], loop[i + 1]])
            if not new_faces:
                break
            dec = trimesh.Trimesh(
                vertices=dec.vertices,
                faces=np.vstack([dec.faces, np.array(new_faces)]))
            dec.merge_vertices()
            dec.update_faces(dec.nondegenerate_faces())
        trimesh.repair.fix_normals(dec)
        if dec.volume < 0:
            dec.invert()
        print(f'  after stitching: watertight={dec.is_watertight}, '
              f'open edges={boundary_stats(dec)}, V={dec.volume:.5f}')
        if dec.is_watertight:
            best = dec
            break

if best is None:
    raise SystemExit('still not watertight')

best.export('/home/bdml/Desktop/OpenFoam12/bem/BROV2_watertight.stl')
json.dump(dict(faces=int(len(best.faces)), watertight=True,
               volume=float(best.volume),
               centroid=[float(c) for c in best.center_mass],
               volume_mc=float(wt.volume), pitch=PITCH),
          open('/home/bdml/Desktop/OpenFoam12/bem/mesh_meta.json', 'w'),
          indent=1)
print(f'SAVED: {len(best.faces)} faces, V={best.volume:.5f} m^3, '
      f'centroid={best.center_mass}')
