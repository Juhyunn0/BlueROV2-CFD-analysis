"""Capytaine sanity check: heave added mass of a deeply submerged sphere.

Theory: m_a = 0.5*rho*V for a sphere far from the free surface.
"""
import numpy as np
import capytaine as cpt

RHO = 1000.0
R = 0.2
V = 4 / 3 * np.pi * R ** 3
theory = 0.5 * RHO * V

mesh = cpt.mesh_sphere(radius=R, center=(0, 0, -10.0), resolution=(40, 40))
body = cpt.FloatingBody(mesh=mesh, name='sphere')
body.add_all_rigid_body_dofs()

solver = cpt.BEMSolver()
print(f'panels: {mesh.nb_faces}')
for omega in (1.0, 5.0):
    pb = cpt.RadiationProblem(body=body, radiating_dof='Heave', omega=omega,
                              rho=RHO, water_depth=np.inf)
    res = solver.solve(pb)
    ma = res.added_masses['Heave']
    err = (ma - theory) / theory * 100
    print(f'omega={omega}: m_a(heave) = {ma:.4f} kg | theory 0.5*rho*V = '
          f'{theory:.4f} kg | error = {err:+.2f}%')
