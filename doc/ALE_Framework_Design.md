# ALE Framework for Incompressible Flow in G+Smo

## Overview

The Arbitrary Lagrangian-Eulerian (ALE) framework extends the incompressible flow solver in G+Smo to handle problems with moving boundaries and deforming domains. This is essential for fluid-structure interaction (FSI) problems and other applications involving time-varying geometries.

## Mathematical Formulation

### ALE Navier-Stokes Equations

In the ALE formulation, the incompressible Navier-Stokes equations are written as:

```
∂u/∂t|_χ + ((u - w)·∇)u - ν∇²u + ∇p = f    (momentum equation)
∇·u = 0                                       (continuity equation)
```

where:
- `u` is the fluid velocity
- `w` is the mesh velocity
- `p` is the pressure
- `ν` is the kinematic viscosity
- `f` is the body force
- `∂/∂t|_χ` denotes the time derivative at fixed reference coordinates

### Geometric Conservation Law (GCL)

To maintain conservation properties on moving meshes, the GCL must be satisfied:

```
dJ/dt + J∇·w = 0
```

where `J` is the Jacobian determinant of the mapping from reference to current configuration.

## Architecture

### Class Hierarchy

```
gsINSSolver (base)
    └── gsINSSolverUnsteady
            └── gsINSSolverALE (new)

gsINSAssembler (base)
    └── gsINSAssemblerUnsteady
            └── gsINSAssemblerALE (new)
```

### Key Components

1. **gsINSSolverALE**: Main ALE solver class
   - Manages mesh deformation
   - Computes mesh velocity
   - Enforces GCL
   - Coordinates ALE time stepping

2. **gsINSAssemblerALE**: ALE-specific assembler
   - Assembles ALE convective terms `((u-w)·∇)u`
   - Handles time derivative on moving mesh
   - Implements GCL enforcement

3. **gsINSVisitorALEConvection**: Visitor for ALE convection
   - Evaluates relative velocity `u - w`
   - Assembles convective terms

4. **gsINSVisitorGCL**: Visitor for GCL
   - Computes GCL residual
   - Adds GCL constraints

## Implementation Details

### Mesh Deformation

The framework uses the existing `gsALE` class from `gsElasticity` module for mesh deformation:

```cpp
// Initialize ALE mesh solver
aleSolver.initializeALE(ale_method::HE); // Harmonic Extension

// Set boundary conditions for mesh deformation
gsBoundaryConditions<> meshBC;
meshBC.addCondition(boundary::moving, condition_type::dirichlet, &displacement);
aleSolver.setMeshDeformationBC(meshBC);

// Update mesh deformation
aleSolver.updateMeshDeformation(boundaryDisp);
```

### Mesh Velocity Computation

Mesh velocity is computed using backward finite differences:

```cpp
w = (d^n - d^{n-1}) / Δt
```

where `d` is the mesh displacement.

### Time Stepping Algorithm

1. Update boundary displacement
2. Solve mesh deformation equation
3. Compute mesh velocity
4. Update geometry and Jacobian
5. Assemble and solve ALE Navier-Stokes system
6. Check GCL satisfaction

## Usage Example

```cpp
// Create ALE solver
gsINSSolverALE<real_t, ColMajor> aleSolver(params);

// Initialize ALE with Harmonic Extension
aleSolver.initializeALE(ale_method::HE);

// Enable GCL enforcement
aleSolver.setEnforceGCL(true);

// Time loop
for (int step = 0; step < nSteps; ++step)
{
    // Update boundary displacement
    gsMultiPatch<> boundaryDisp = computeBoundaryMotion(time);
    
    // Update mesh
    aleSolver.updateMeshDeformation(boundaryDisp);
    
    // Solve flow on deformed mesh
    aleSolver.nextIteration();
    
    // Check GCL
    real_t gclError = aleSolver.checkGCL();
}
```

## Mesh Deformation Methods

The framework supports various mesh deformation methods through `ale_method`:

- `HE`: Harmonic Extension (Laplace equation)
- `IHE`: Incremental Harmonic Extension
- `BHE`: Biharmonic Extension
- `TINE`: Tangential-Incremental-Nonlinear-Elasticity
- `TINE_StVK`: TINE with St. Venant-Kirchhoff material

## Boundary Conditions

### Flow Boundary Conditions

Standard incompressible flow BCs are supported:
- Dirichlet (velocity)
- Neumann (traction)
- Do-nothing (natural outlet)

### Mesh Deformation Boundary Conditions

- Dirichlet: prescribed displacement
- Neumann: prescribed traction (for elastic mesh deformation)
- Sliding: tangential motion only

## Performance Considerations

1. **Mesh Quality**: Monitor mesh quality during deformation
2. **Time Step**: May need smaller time steps for large deformations
3. **GCL Enforcement**: Adds computational cost but improves conservation
4. **Mesh Update Frequency**: Can update mesh every few time steps for efficiency

## Future Extensions

1. **Adaptive Mesh Refinement**: Refine mesh in regions of high deformation
2. **Mesh Quality Indicators**: Automatic detection of mesh degradation
3. **Parallel Implementation**: Domain decomposition for large-scale problems
4. **Higher-Order Time Integration**: BDF2, Crank-Nicolson schemes
5. **FSI Coupling**: Direct integration with structural solvers

## References

1. Donea, J., Huerta, A., Ponthot, J. P., & Rodríguez‐Ferran, A. (2004). "Arbitrary Lagrangian–Eulerian Methods"
2. Formaggia, L., & Nobile, F. (1999). "A stability analysis for the arbitrary Lagrangian Eulerian formulation with finite elements"
3. Thomas, P. D., & Lombard, C. K. (1979). "Geometric conservation law and its application to flow computations on moving grids"