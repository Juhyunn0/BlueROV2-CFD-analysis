"""Build a watertight ROV surface for BEM from the raw CAD assembly STL.

Method: surface voxelization at 5 mm (vertex binning + dense area-weighted
surface sampling), 1-voxel binary closing to seal pinholes, exterior flood
fill -> solid mask, marching cubes, quadric decimation to ~12k faces.
"""
import json

import numpy as np
import trimesh
from scipy import ndimage
from skimage import measure

SRC = '/home/bdml/Desktop/OpenFoam12/rovCase/constant/geometry/BROV2.stl'
OUT = '/home/bdml/Desktop/OpenFoam12/bem/BROV2_watertight.stl'
PITCH = 0.005

mesh = trimesh.load_mesh(SRC)
print(f'raw STL: {len(mesh.faces)} faces, {len(mesh.vertices)} vertices')
print(f'watertight: {mesh.is_watertight}')
print(f'bbox: {mesh.bounds[0]} .. {mesh.bounds[1]}')

# --- surface voxelization ---
lo = mesh.bounds[0] - 3 * PITCH
hi = mesh.bounds[1] + 3 * PITCH
dims = np.ceil((hi - lo) / PITCH).astype(int) + 1
print(f'voxel grid: {dims} ({dims.prod()/1e6:.1f}M voxels)')

occ = np.zeros(dims, dtype=bool)


def mark(points):
    idx = np.floor((points - lo) / PITCH).astype(int)
    idx = np.clip(idx, 0, np.array(dims) - 1)
    occ[idx[:, 0], idx[:, 1], idx[:, 2]] = True


mark(mesh.vertices)
# dense area-weighted samples to cover large flat panels
n_samp = int(mesh.area / PITCH ** 2 * 8)
print(f'sampling {n_samp/1e6:.1f}M surface points (area {mesh.area:.3f} m^2)')
for _ in range(4):
    pts, _ = trimesh.sample.sample_surface(mesh, n_samp // 4)
    mark(pts)

# seal 1-voxel pinholes so the exterior flood cannot leak into sealed bodies
occ_closed = ndimage.binary_closing(occ, structure=np.ones((3, 3, 3)))

# exterior flood fill from the corner
ext = ~occ_closed
lbl, _ = ndimage.label(ext)
outside = lbl == lbl[0, 0, 0]
solid = ~outside
print(f'solid voxels: {solid.sum()} ({solid.sum()*PITCH**3*1000:.2f} L)')

# --- marching cubes -> watertight surface ---
verts, faces, _, _ = measure.marching_cubes(solid.astype(np.uint8), level=0.5)
wt = trimesh.Trimesh(vertices=verts * PITCH + lo, faces=faces)
trimesh.repair.fix_normals(wt)
if wt.volume < 0:
    wt.invert()
print(f'marching cubes: {len(wt.faces)} faces, watertight={wt.is_watertight}, '
      f'V={wt.volume:.5f} m^3')

# --- decimate ---
dec = wt.simplify_quadric_decimation(face_count=12000)
trimesh.repair.fix_normals(dec)
if dec.volume < 0:
    dec.invert()
print(f'decimated: {len(dec.faces)} faces, watertight={dec.is_watertight}, '
      f'V={dec.volume:.5f} m^3, centroid={dec.center_mass}')

dec.export(OUT)
meta = dict(faces=int(len(dec.faces)), watertight=bool(dec.is_watertight),
            volume=float(dec.volume), centroid=[float(c) for c in dec.center_mass],
            volume_mc=float(wt.volume), pitch=PITCH)
with open('/home/bdml/Desktop/OpenFoam12/bem/mesh_meta.json', 'w') as f:
    json.dump(meta, f, indent=1)
print('saved', OUT)
