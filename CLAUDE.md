
# ROV Hydrodynamics Project — Meshing & CFD

## Environment (must follow)

- Ubuntu 22.04, OpenFOAM 12 (Foundation fork), installed at /opt/openfoam12
- Non-interactive shells do NOT source ~/.bashrc.
    Every OpenFOAM command must be run as: `source /opt/openfoam12/etc/bashrc && <command>`
- Foundation-fork caveats:
    - No simpleFoam/pimpleFoam → use `foamRun -solver incompressibleFluid`
    - Feature extraction: `surfaceFeatures` + system/surfaceFeaturesDict
    - Turbulence setup lives in constant/momentumTransport (NOT the older turbulenceProperties)
    - Most web tutorials use older or .com-fork syntax — always cross-check
      against the actual files in $FOAM_TUTORIALS for this version's syntax
- Best template: $FOAM_TUTORIALS/incompressibleFluid/motorBike
    (the canonical STL external-flow + snappyHexMesh + kOmegaSST example.
     Default approach: copy its dicts and modify them.)

## Working rules

- Save the log of every step as log.<command></command> and report only the key findings
- Check the STL with surfaceCheck before meshing (watertight, bounding box, open edges)
- If the bounding box is in the hundreds, assume mm and convert to m,
    but ask for my confirmation before converting
- checkMesh pass criteria: "Mesh OK", maxNonOrtho < 65, maxSkewness < 4
- If the same step fails 3 times, stop and report the situation
- Before any run expected to take >10 min, state the estimated time
- Parallel runs: check cores with nproc → decomposePar → mpirun -np N <tool></tool> -parallel → reconstructPar
