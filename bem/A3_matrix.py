"""6x6 added-mass matrix of the ROV about its volume centroid (capytaine).

Vehicle frame: x = surge (CAD +z), y = sway (CAD +x), z = heave (CAD +y, up).
Body deeply submerged (centroid 10 m below the free surface), rho = 1000.
"""
import json

import numpy as np
import trimesh
import capytaine as cpt

RHO = 1000.0
DEPTH = 10.0
DOFS = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']

wt = trimesh.load_mesh('/home/bdml/Desktop/OpenFoam12/bem/BROV2_watertight.stl')
assert wt.is_watertight
centroid_cad = wt.center_mass.copy()

# permute CAD (x,y,z) -> vehicle (x_v, y_v, z_v) = (z, x, y): right-handed
verts = wt.vertices[:, [2, 0, 1]].copy()
faces = wt.faces.copy()
centroid_v = centroid_cad[[2, 0, 1]]
verts -= centroid_v          # centroid to origin
verts[:, 2] -= DEPTH         # submerge

F4 = np.column_stack([faces, faces[:, 2]])   # tri -> capytaine quad format
mesh = cpt.Mesh(vertices=verts, faces=F4, name='BROV2')
mesh.heal_mesh()
body = cpt.FloatingBody(mesh=mesh, name='BROV2')
body.rotation_center = np.array([0.0, 0.0, -DEPTH])
body.add_all_rigid_body_dofs()

solver = cpt.BEMSolver()
mats = {}
for omega in (1.0, 5.0):
    probs = [cpt.RadiationProblem(body=body, radiating_dof=d, omega=omega,
                                  rho=RHO, water_depth=np.inf) for d in DOFS]
    results = solver.solve_all(probs, progress_bar=False)
    A = np.zeros((6, 6))
    for res in results:
        j = DOFS.index(res.radiating_dof)
        for i, d in enumerate(DOFS):
            A[i, j] = res.added_masses[d]
    mats[omega] = A
    print(f'\n=== added-mass matrix, omega = {omega} rad/s ===')
    hdr = '        ' + ''.join(f'{d:>10}' for d in DOFS)
    print(hdr)
    for i, d in enumerate(DOFS):
        print(f'{d:>7} ' + ''.join(f'{A[i, j]:10.3f}' for j in range(6)))

rel = np.max(np.abs(mats[5.0] - mats[1.0]) / (np.abs(mats[1.0]) + 1e-6))
print(f'\nmax relative difference omega=1 vs 5: {rel*100:.3f}%')
print(f'volume V = {wt.volume:.5f} m^3, rho*V = {RHO*wt.volume:.2f} kg')
print(f'centroid (CAD frame) = {centroid_cad}')

json.dump(dict(A1=mats[1.0].tolist(), A5=mats[5.0].tolist(),
               V=float(wt.volume), centroid_cad=centroid_cad.tolist(),
               dofs=DOFS),
          open('/home/bdml/Desktop/OpenFoam12/bem/added_mass.json', 'w'),
          indent=1)
print('saved added_mass.json')
