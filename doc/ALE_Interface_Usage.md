# ALE Interface Usage Guide

## Overview

The ALE (Arbitrary Lagrangian-Eulerian) interface extends the gsIncompressibleFlow module to handle problems with moving boundaries and deforming domains. This implementation provides:

- `gsINSSolverALE`: ALE solver extending the unsteady Navier-Stokes solver
- `gsINSAssemblerALE`: ALE-specific assembler for modified NS equations
- ALE-specific visitors for convection and GCL terms

## Key Classes

### gsINSSolverALE

Main ALE solver class that manages:
- Mesh deformation tracking
- Mesh velocity computation
- Geometry updates
- GCL enforcement

### gsINSAssemblerALE

Assembler that handles:
- ALE convection terms ((u-w)·∇)u
- Time derivatives on moving mesh
- GCL enforcement terms

## Usage Example

```cpp
// Create flow problem
gsNavStokesPde<real_t> NSpde(patches, bcInfo, &sourceFunc, viscosity);
gsFlowSolverParams<real_t> params(NSpde, discreteBases);

// Create ALE solver
gsINSSolverALE<real_t, ColMajor> aleSolver(params);

// Initialize ALE framework
aleSolver.initializeALE();
aleSolver.setEnforceGCL(true);  // Optional: enforce GCL

// Set mesh deformation boundary conditions
gsBoundaryConditions<> meshBC;
meshBC.addCondition(boundary::moving, condition_type::dirichlet, &displacement);
aleSolver.setMeshDeformationBC(meshBC);

// Initialize solver
aleSolver.initialize();

// Time loop
for (int step = 0; step < nSteps; ++step)
{
    // Get mesh displacement (from external source or computation)
    gsMultiPatch<> meshDisplacement = computeMeshDisplacement(time);
    
    // Update mesh
    aleSolver.updateMeshDeformation(meshDisplacement);
    
    // Solve flow on deformed mesh
    aleSolver.nextIteration();
    
    // Check GCL error (optional)
    real_t gclError = aleSolver.checkGCL();
}
```

## Integration with Mesh Deformation

The ALE interface is designed to work with external mesh deformation solvers. The mesh displacement should be computed separately and provided to the ALE solver at each time step.

For integration with gsElasticity's gsALE:

```cpp
// Create gsALE solver for mesh deformation
gsALE<real_t> meshSolver(geometry, displacement, interface, ale_method::HE);

// Solve for mesh displacement
gsMultiPatch<real_t> meshDisplacement;
meshSolver.constructSolution(meshDisplacement);

// Pass to ALE flow solver
aleSolver.updateMeshDeformation(meshDisplacement);
```

## Important Notes

1. **Mesh Quality**: Monitor mesh quality during deformation. Large deformations may require remeshing.

2. **Time Step**: ALE simulations may require smaller time steps for stability, especially with large mesh velocities.

3. **GCL Enforcement**: Enabling GCL enforcement improves conservation properties but adds computational cost.

4. **Boundary Conditions**: Both flow and mesh boundary conditions must be carefully specified for consistency.

## Example: Oscillating Channel

See `examples/gsIncompressibleFlow_ALE.cpp` for a complete example of flow in a channel with oscillating walls.

To run:
```bash
./gsIncompressibleFlow_ALE --ale --geo=4 --plot --animStep=5
```

This simulates flow around a flapping beam using the ALE formulation.