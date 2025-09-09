# ALE Implementation for gsIncompressibleFlow

## Summary

This document summarizes the implementation of ALE (Arbitrary Lagrangian-Eulerian) formulation in the gsIncompressibleFlow module for fluid-structure interaction (FSI) simulations.

## Implementation Overview

### Key Files Added:

1. **gsINSTermsALE.h** - ALE-specific terms
   - `gsINSTerm_ALEConvection`: Implements the ALE convection term `((u - u_mesh) · ∇)φ`

2. **gsINSVisitorsALE.h** - ALE-specific visitors
   - `gsINSVisitorUUnonlinALE`: Visitor for ALE nonlinear terms

3. **gsINSAssemblerALE.h** - ALE-specific assemblers
   - `gsINSAssemblerUnsteadyALE`: Extends unsteady assembler with ALE capabilities

4. **gsINSSolverALE.h** - ALE-specific solver
   - `gsINSSolverUnsteadyALE`: Main solver class for ALE simulations
   - `gsFSIHelper`: Helper class for FSI coupling (framework only)

5. **examples/ale_fsi_example.cpp** - Example usage

## Key Modifications

### 1. ALE Convection Term
The standard convection term `(u · ∇)u` is replaced with `((u - u_mesh) · ∇)u` where `u_mesh` is the mesh velocity.

### 2. Mesh Velocity Management
- Mesh velocity is computed as: `u_mesh = (d_new - d_old) / dt`
- Mesh displacement and velocity fields are stored and updated each time step

### 3. Assembly Process
The ALE assembler:
- Maintains mesh velocity and displacement fields
- Uses ALE visitor when ALE is active
- Falls back to standard formulation when ALE is inactive

## Usage Example

```cpp
// Create parameters and solver
gsFlowSolverParams<real_t> params(nsPde, discreteBases);
gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));

// Activate ALE
solver.setALEActive(true);

// Define mesh motion (e.g., from structure solver)
auto meshMotion = [](real_t t) -> gsMatrix<> {
    // Return mesh displacement at time t
};
solver.setMeshUpdateFunction(meshMotion);

// Time stepping
for (int step = 0; step < nSteps; ++step) {
    solver.nextIteration();
}
```

## Technical Details

### ALE Convection Assembly
The ALE convection term is assembled by:
1. Computing relative velocity: `u_rel = u - u_mesh`
2. Transforming gradients to physical space
3. Computing convection: `(u_rel · ∇)φ_trial`
4. Assembling: `∫ φ_test * (u_rel · ∇)φ_trial dΩ`

### Mesh Update Process
1. Get new mesh displacement from FSI coupling
2. Compute mesh velocity: `(d_new - d_old) / dt`
3. Update mesh fields in assembler
4. Use updated fields in next assembly

## Limitations and Future Work

1. **Geometric Conservation Law**: Not yet implemented
2. **Mesh Quality**: No automatic remeshing or quality checks
3. **FSI Coupling**: Only framework provided, full implementation needed
4. **Time Integration**: Currently uses same scheme as fluid solver

## Compilation

To compile the example:
```bash
cd build
make ale_fsi_example
```

## Testing

Run the example:
```bash
./bin/ale_fsi_example
```

This creates output files for visualization in ParaView:
- `fluid_velocity_*.vts` - Fluid velocity field
- `fluid_pressure_*.vts` - Pressure field  
- `mesh_displacement_*.vts` - Mesh displacement field