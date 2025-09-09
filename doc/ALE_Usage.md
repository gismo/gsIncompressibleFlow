# ALE Formulation for Incompressible Flow

## Overview

The ALE (Arbitrary Lagrangian-Eulerian) formulation has been added to gsIncompressibleFlow to enable fluid-structure interaction (FSI) simulations. The implementation modifies the convection terms to account for mesh motion.

## Key Components

### 1. ALE Terms (`gsINSTermsALE.h`)
- `gsINSTerm_ALEConvection`: Implements convection with mesh velocity: `((u - u_mesh) · ∇)φ`
- `gsINSTerm_ALEUnsteady`: Implements ALE time derivative with mesh convection correction

### 2. ALE Visitors (`gsINSVisitorsALE.h`)
- `gsINSVisitorUUnonlinALE`: Assembles nonlinear terms with ALE convection
- `gsINSVisitorUnsteadyALE`: Assembles unsteady terms with ALE correction

### 3. ALE Assemblers (`gsINSAssemblerALE.h`)
- `gsINSAssemblerALE`: Base ALE assembler
- `gsINSAssemblerUnsteadyALE`: Unsteady ALE assembler with mesh field management

### 4. ALE Solver (`gsINSSolverALE.h`)
- `gsINSSolverUnsteadyALE`: Main solver class for ALE simulations
- `gsFSIHelper`: Utility class for FSI coupling operations

## Usage Example

```cpp
// Create solver parameters
gsFlowSolverParams<> params;
params.setPde(pde);
params.addBasis(basis);  // Velocity
params.addBasis(basis);  // Pressure

// Create ALE solver
gsINSSolverUnsteadyALE<> solver(memory::make_shared(params));

// Activate ALE
solver.setALEActive(true);

// Define mesh motion function
auto meshMotion = [](real_t t) -> gsMatrix<> {
    // Return mesh displacement at time t
    // This typically comes from structure solver
};
solver.setMeshUpdateFunction(meshMotion);

// Time stepping
for (int step = 0; step < nSteps; ++step) {
    solver.nextIteration();
}
```

## FSI Coupling

For FSI simulations:

1. **Mesh Update**: 
   - Get structure displacement at interface
   - Extend displacement into fluid domain (e.g., harmonic extension)
   - Update fluid mesh using `solver.updateMesh(displacement)`

2. **Force Transfer**:
   - Compute fluid stress at interface
   - Transfer traction to structure solver

3. **Coupling Schemes**:
   - Explicit: Update mesh once per time step
   - Implicit: Iterate between fluid and structure until convergence

## Implementation Notes

- Mesh velocity is computed as: `u_mesh = (d_new - d_old) / dt`
- ALE convection term: `(u - u_mesh) · ∇u` instead of `u · ∇u`
- Geometric conservation law should be satisfied for consistency

## Future Work

- Implement geometric conservation law
- Add mesh quality checks and remeshing
- Implement strong FSI coupling schemes
- Add parallel support for large-scale FSI