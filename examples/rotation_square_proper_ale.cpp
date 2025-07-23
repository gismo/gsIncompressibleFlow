/** @file rotation_square_proper_ale.cpp
    @brief Example of rotating square using proper ALE with elastic mesh deformation
    
    This implementation uses the elastic equation approach for mesh deformation
    instead of decay functions. The inner boundaries follow the rotating square
    exactly while outer boundaries remain fixed.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>
#include <gsHLBFGS/gsHLBFGS.h>
// #include <gsElasticity/src/gsALE.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    // Problem parameters
    real_t Re = 100;          // Reynolds number
    real_t rho = 1000.0;      // Fluid density
    real_t mu = 0.001;        // Dynamic viscosity
    real_t meanVelocity = 1.0;// Mean inlet velocity
    
    // Time parameters
    real_t timeSpan = 0.5;    // Total simulation time
    real_t timeStep = 0.01;   // Time step size
    real_t maxAngle = 20;     // Maximum rotation angle
    index_t nTimeSteps = timeSpan / timeStep;
    
    // Domain and discretization
    index_t numRefine = 3;
    index_t degree = 2;
    index_t outputInterval = 10;  // Output every N steps
    
    // Command line parsing
    gsCmdLine cmd("Testing rotating square with proper ALE elastic mesh deformation.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numRefine);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("a","angle","Maximum rotation angle",maxAngle);
    cmd.addInt("o","output","Output interval (save every N steps)",outputInterval);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    gsInfo << "=== Rotating Square Proper ALE Example ===\n";
    gsInfo << "Reynolds number: " << Re << "\n";
    gsInfo << "Time step: " << timeStep << "\n";
    gsInfo << "Time span: " << timeSpan << "\n";
    gsInfo << "Max rotation angle: " << maxAngle << "\n\n";
    
    // Create fluid domain - cross shape with 8 patches
    /*
     * ______________________
     * |      |      |      |
     * |  2   |  4   |  7   |
     * |______|______|______|
     * |      |      |      |
     * |  1   |  sq  |  6   |
     * |______|______|______|
     * |      |      |      |
     * |  0   |  3   |  5   |
     * |______|______|______|
     */
    gsMultiPatch<> fluidDomain;
    real_t L = 1.0; // Size of each patch
    
    // Create patches
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // patch 0
    fluidDomain.patch(0).coefs() << -1,-1, 0,-1, -1,0, 0,0;
    
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // patch 1
    fluidDomain.patch(1).coefs() << -1,0, 0,0, -1,1, 0,1;
    
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // patch 2  
    fluidDomain.patch(2).coefs() << -1,1, 0,1, -1,2, 0,2;
    
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // patch 3
    fluidDomain.patch(3).coefs() << 0,-1, 1,-1, 0,0, 1,0;
    
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // patch 4
    fluidDomain.patch(4).coefs() << 0,1, 1,1, 0,2, 1,2;
    
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // patch 5
    fluidDomain.patch(5).coefs() << 1,-1, 2,-1, 1,0, 2,0;
    
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // patch 6
    fluidDomain.patch(6).coefs() << 1,0, 2,0, 1,1, 2,1;
    
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // patch 7
    fluidDomain.patch(7).coefs() << 1,1, 2,1, 1,2, 2,2;
    
    fluidDomain.computeTopology();
    
    // Create central square for visualization (not part of fluid domain)
    gsMultiPatch<> centralSquare;
    centralSquare.addPatch(gsNurbsCreator<>::BSplineSquare(1));
    centralSquare.patch(0).coefs() << 0,0, 1,0, 0,1, 1,1;

    
    // Refine mesh
    for (index_t i = 0; i < numRefine; ++i)
    {
        fluidDomain.uniformRefine();
        centralSquare.uniformRefine();
    }
    
    // Degree elevate for Taylor-Hood elements
    fluidDomain.degreeElevate();
    centralSquare.degreeElevate();
    
    // Create basis
    gsMultiBasis<> basis(fluidDomain);
    gsMultiBasis<> basisPressure(fluidDomain);
    
    // Setup boundary conditions for flow
    gsBoundaryConditions<> bcInfo;
    
    // Velocity boundary conditions
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    gsConstantFunction<> inletVel(meanVelocity, 0.0, 2);
    
    // Inlet (left side)
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &inletVel, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, &inletVel, 0);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, &inletVel, 0);
    
    // External walls (top and bottom) - fixed walls
    for (index_t d = 0; d < 2; ++d)
    {
        bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(3, boundary::south, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(4, boundary::north, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(5, boundary::south, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(7, boundary::north, condition_type::dirichlet, &zeroVel, 0, d);
    }
    
    // Note: Inner boundaries (interfaces with rotating square) have NO explicit 
    // velocity BCs in ALE formulation. The no-slip condition is automatically 
    // satisfied through the relative velocity (u-w).
    
    // Outlet (right side) - natural boundary condition (do-nothing)
    
    // Create PDE
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &zeroVel, mu);
    
    // Setup solver parameters
    std::vector<gsMultiBasis<>> discreteBases;
    discreteBases.push_back(basis);  // Velocity basis
    discreteBases.push_back(basisPressure);  // Pressure basis
    
    gsFlowSolverParams<real_t> params(nsPde, discreteBases);
    
    // Solver options
    params.options().setInt("nonlin.maxIt", 100);     // Max Picard iterations
    params.options().setReal("nonlin.tol", 1e-6);     // Picard tolerance
    params.options().setInt("lin.maxIt", 200);        // Max linear solver iterations
    params.options().setReal("lin.tol", 1e-8);        // Linear solver tolerance
    params.options().setReal("timeStep", timeStep);   // Time step size
    params.options().setSwitch("quiet", false);       // Enable verbose output
    
    // Create ALE solver
    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));
    
    // Initialize solver
    solver.initialize();
    
    // Activate ALE after initialization
    solver.setALEActive(true);
    
    // Enable mesh optimization using gsBarrierPatch
    solver.setMeshOptimization(true);
    solver.getMeshOptOptions().setInt("Verbose", 1);  // Enable verbose output for mesh optimization
    solver.getMeshOptOptions().setInt("ParamMethod", 1);
    solver.getMeshOptOptions().setInt("AAPreconditionType", 0);
    
    // Enable dynamic boundary mapping for rotating domains
    solver.setDynamicBoundaryMapping(false); // This enables special handling for rotating boundaries
    
    // Set rotation parameters for dynamic boundary mapping
    gsVector<real_t> rotationCenter(2);
    rotationCenter << 0.5, 0.5; // Center of the domain where square rotates
    solver.setRotationParameters(10.0, rotationCenter); // 10s period
    
    gsInfo << "Mesh optimization: " << (solver.isMeshOptimizationEnabled() ? "enabled" : "disabled") << "\n";
    gsInfo << "Dynamic boundary mapping: " << (solver.isDynamicBoundaryMappingEnabled() ? "enabled" : "disabled") << "\n";
    
    // Store original mesh coordinates
    gsMultiPatch<> originalFluidDomain = fluidDomain;
    gsMultiPatch<> centralSquareOld = centralSquare;
    gsMultiPatch<> centralSquareNew = centralSquare;
    
    // Store the actual fluid domain that will be updated each step
    gsMultiPatch<> currentFluidDomain = fluidDomain;
    
    // Setup boundary conditions for mesh motion
    // We need to specify which boundaries follow the rotating square
    // and which remain fixed
    gsBoundaryConditions<> meshBC;
    
    // Fixed outer boundaries - zero displacement
    gsConstantFunction<> zeroDisp(0.0, 0.0, 2);
    
    // External boundaries remain fixed
    for (index_t d = 0; d < 2; ++d)
    {
        // Outer boundaries
        meshBC.addCondition(0, boundary::west, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(0, boundary::south, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(1, boundary::west, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(2, boundary::west, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(2, boundary::north, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(3, boundary::south, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(4, boundary::north, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(5, boundary::south, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(5, boundary::east, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(6, boundary::east, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(7, boundary::north, condition_type::dirichlet, &zeroDisp, 0, d);
        meshBC.addCondition(7, boundary::east, condition_type::dirichlet, &zeroDisp, 0, d);
    }
    
    // Define mesh motion function for proper ALE
    // This function returns the CUMULATIVE displacement from the original mesh
    // We'll use a simplified approach: directly compute displacement for ALL boundary points
    auto meshMotion = [&](real_t t)
    {
        // Calculate rotation angle with sinusoidal motion
        const real_t angle = 2.2 * maxAngle * std::sin(2 * M_PI / timeSpan * t) * M_PI / 180.0; // Convert to radians
        const real_t cos_angle = std::cos(angle);
        const real_t sin_angle = std::sin(angle);
        const real_t cx = 0.5;  // rotation center x
        const real_t cy = 0.5;  // rotation center y

        const index_t udofs = solver.getAssembler()->getUdofs();
        gsMatrix<> disp(2*udofs,1); 
        disp.setZero();

        const gsDofMapper& mapper = solver.getAssembler()->getMappers()[0];
        
        // We'll set displacement for ALL free DOFs that are on inner boundaries
        // The elastic solver will handle the propagation to interior points
        
        // Set displacement for inner boundaries (interfaces with rotating square)
        // These boundaries must follow the rotating square exactly
        
        // Inner boundaries of the 8 patches that touch the central square:
        // Patch 0: east boundary touches square bottom
        // Patch 1: east boundary touches square left  
        // Patch 3: north boundary touches square bottom
        // Patch 3: west boundary touches square right
        // Patch 4: south boundary touches square top
        // Patch 6: west boundary touches square right
        // Patch 7: west boundary touches square top
        
        // Helper function to set displacement for a boundary
        auto setBoundaryDisplacement = [&](index_t patch, boxSide side)
        {
            const gsMatrix<>& originalCoefs = originalFluidDomain.patch(patch).coefs();
            const gsBasis<>& patchBasis = basis.basis(patch);
            
            // Get boundary indices
            gsMatrix<index_t> bnd = patchBasis.boundary(side);
            
            for (index_t k = 0; k < bnd.size(); ++k)
            {
                index_t i = bnd(k);
                
                if (!mapper.is_free(i, patch))
                    continue; // Skip constrained DOFs
                
                real_t x = originalCoefs(i,0);
                real_t y = originalCoefs(i,1);
                
                // Check if this point is on the inner boundary (touches the rotating square)
                // The inner square has corners at (0, 0) and (1, 1)
                const real_t tol = 1e-10;
                bool isOnInnerBoundary = false;
                
                // Check which side of the inner square this boundary touches
                if (patch == 1 && side == boundary::east && std::abs(x - 0.0) < tol)
                    isOnInnerBoundary = (y >= 0.0 - tol && y <= 1.0 + tol);
                else if (patch == 6 && side == boundary::west && std::abs(x - 1.0) < tol)
                    isOnInnerBoundary = (y >= 0.0 - tol && y <= 1.0 + tol);
                else if (patch == 3 && side == boundary::north && std::abs(y - 0.0) < tol)
                    isOnInnerBoundary = (x >= 0.0 - tol && x <= 1.0 + tol);
                else if (patch == 4 && side == boundary::south && std::abs(y - 1.0) < tol)
                    isOnInnerBoundary = (x >= 0.0 - tol && x <= 1.0 + tol);
                
                if (isOnInnerBoundary)
                {
                    // Apply rotation transformation only to points on the inner boundary
                    real_t xr = cos_angle * (x - cx) - sin_angle * (y - cy) + cx;
                    real_t yr = sin_angle * (x - cx) + cos_angle * (y - cy) + cy;
                    
                    // Displacement
                    real_t dx = xr - x;
                    real_t dy = yr - y;
                    
                    index_t base = mapper.index(i, patch);
                    disp(base) = dx;
                    disp(base + udofs) = dy;
                    
                    // Debug: count how many points are set
                    static int count = 0;
                    if (++count <= 5) // Only print first 5
                    {
                        gsInfo << "Setting displacement for patch " << patch << " side " << side 
                               << " at (" << x << "," << y << ") -> dx=" << dx << ", dy=" << dy << "\n";
                    }
                }
                // Points not on inner boundary remain at zero displacement (already initialized)
            }
        };
        
        // Set displacements for all inner boundaries that touch the central square
        // Based on the mesh layout:
        setBoundaryDisplacement(1, boundary::east);   // patch 1 east side touches square left
        setBoundaryDisplacement(3, boundary::north);  // patch 3 north side touches square bottom
        setBoundaryDisplacement(4, boundary::south);  // patch 4 south side touches square top
        setBoundaryDisplacement(6, boundary::west);   // patch 6 west side touches square right
        
        // Debug: also try setting for all possible boundaries
        if (false) // Enable for debugging
        {
            setBoundaryDisplacement(0, boundary::east);
            setBoundaryDisplacement(2, boundary::east);
            setBoundaryDisplacement(5, boundary::west);
            setBoundaryDisplacement(7, boundary::west);
        }
        
        // Now we need to solve an elastic problem to get the interior mesh displacement
        // The ALE module will handle this automatically using the boundary displacements
        
        return disp;
    };
    
    solver.setMeshUpdateFunction(meshMotion);
    
    // Output initial mesh configuration
    gsInfo << "Writing initial mesh configuration...\n";
    gsWriteParaview(fluidDomain, "proper_ale_initial_fluid_mesh", 1000, true);
    gsWriteParaview(centralSquare, "proper_ale_initial_central_square", 1000);
    
    // Solve Stokes problem for initial conditions
    gsInfo << "Solving initial Stokes problem...\n";
    solver.solveStokes();
    
    // Time stepping loop
    gsInfo << "Starting time integration with proper ALE...\n";
    
    // Prepare paraview collection for time series output
    gsParaviewCollection collectionFlow("proper_ale_rotating_square_flow");
    gsParaviewCollection collectionMesh("proper_ale_rotating_square_mesh");
    
    // Update recalculation of nTimeSteps after command line parsing
    nTimeSteps = timeSpan / timeStep;
    
    for (index_t step = 0; step <= nTimeSteps; ++step)
    {
        real_t time = step * timeStep;
        gsInfo << "\nTime step " << step << ", t = " << time << "\n";
        
        // Update central square position for visualization
        real_t angle = 2.2 * maxAngle * std::sin(2 * M_PI / timeSpan * time);
        centralSquareNew = centralSquare;
        gsNurbsCreator<>::rotate2D(centralSquareNew.patch(0), angle, 0.5, 0.5);
        
        // Solve fluid problem with ALE
        if (step > 0)
        {
            solver.nextIteration();
        }
        
        // Get solution
        gsMatrix<> velCoefs = solver.solutionCoefs(0);
        gsMatrix<> presCoefs = solver.solutionCoefs(1);
        
        // Output some statistics
        gsInfo << "  Max velocity: " << velCoefs.lpNorm<gsEigen::Infinity>() << "\n";
        gsInfo << "  Max pressure: " << presCoefs.lpNorm<gsEigen::Infinity>() << "\n";
        
        if (solver.isALEActive() && step > 0)
        {
            gsField<> meshDispField = solver.getMeshDisplacementField();
            gsInfo << "  ‖meshDisp‖_∞: "
                   << meshDispField.coefficientVector().template lpNorm<gsEigen::Infinity>() << "\n";
                   
            // Get the actual deformed mesh from the solver
            // The solver should already have the updated mesh after nextIteration()
            const gsMultiPatch<>& solverPatches = solver.getAssembler()->getPatches();
            currentFluidDomain = solverPatches;
            
            // Note: gsBarrierPatch optimization is now handled inside the solver
            // if mesh optimization is enabled
            
            // Output deformed mesh only every outputInterval steps
            if (step % outputInterval == 0)
            {
                std::string meshname = "proper_ale_deformed_mesh_" + std::to_string(step);
                std::string dispname = "proper_ale_mesh_disp_" + std::to_string(step);
                
                // Write mesh with all patches
                gsWriteParaview(currentFluidDomain, meshname, 1000, true);
                gsWriteParaview(meshDispField, dispname, 1000, true);
                
                // Add all patches to collection
                for (size_t p = 0; p < currentFluidDomain.nPatches(); ++p)
                {
                    std::ostringstream meshnameWithPatch, dispnameWithPatch;
                    meshnameWithPatch << meshname << "_" << p << ".vts";
                    dispnameWithPatch << dispname << "_" << p << ".vts";
                    collectionMesh.addPart(meshnameWithPatch.str(), time, "mesh", p);
                    collectionMesh.addPart(dispnameWithPatch.str(), time, "displacement", p);
                }
                
                // Also output the central square for visualization
                gsWriteParaview(centralSquareNew, "proper_ale_central_square_" + std::to_string(step), 1000);
            }
            
            // Update fluidDomain for subsequent output
            fluidDomain = currentFluidDomain;
        }
        
        // Export solution periodically
        if (step % outputInterval == 0)
        {
            // Get velocity and pressure coefficient vectors from already computed solution
            gsMatrix<> velCoefsFull = velCoefs;
            gsMatrix<> presCoefsFull = presCoefs;
            
            // Create multipatch for velocity (2D vector field)
            gsMultiPatch<> velPatches(fluidDomain);
            gsMultiPatch<> presPatches(fluidDomain);
            
            // Update coefficients for each patch
            const gsDofMapper& uMapper = solver.getAssembler()->getMappers()[0];
            const gsDofMapper& pMapper = solver.getAssembler()->getMappers()[1];
            const index_t shift = solver.getAssembler()->getUdofs();
            
            for (size_t p = 0; p < fluidDomain.nPatches(); ++p)
            {
                // Velocity field (2 components)
                index_t nCoefs = basis.basis(p).size();
                gsMatrix<> velPatchCoefs(nCoefs, 2);
                for (index_t i = 0; i < nCoefs; ++i)
                {
                    if (uMapper.is_free(i, p))
                    {
                        index_t idx = uMapper.index(i, p);
                        velPatchCoefs(i, 0) = velCoefsFull(idx);
                        velPatchCoefs(i, 1) = velCoefsFull(idx + shift);
                    }
                    else
                    {
                        velPatchCoefs(i, 0) = 0.0;
                        velPatchCoefs(i, 1) = 0.0;
                    }
                }
                velPatches.patch(p).coefs() = velPatchCoefs;
                
                // Pressure field (scalar)
                index_t nPCoefs = basisPressure.basis(p).size();
                gsMatrix<> presPatchCoefs(nPCoefs, 1);
                for (index_t i = 0; i < nPCoefs; ++i)
                {
                    if (pMapper.is_free(i, p))
                    {
                        index_t idx = pMapper.index(i, p);
                        presPatchCoefs(i, 0) = presCoefsFull(idx);
                    }
                    else
                    {
                        presPatchCoefs(i, 0) = 0.0;
                    }
                }
                presPatches.patch(p).coefs() = presPatchCoefs;
            }
            
            // Create fields
            gsField<> velocity(fluidDomain, velPatches);
            gsField<> pressure(fluidDomain, presPatches);
            
            // Output flow fields
            unsigned samp = 1000; // sampling points per direction
            std::string vname = "proper_ale_velocity_" + std::to_string(step);
            std::string pname = "proper_ale_pressure_" + std::to_string(step);
            
            // Write velocity and pressure fields with all patches
            gsWriteParaview(velocity, vname, samp, true); // true = mesh flag
            gsWriteParaview(pressure, pname, samp, true);
            
            // Add all patches to collection
            // Files are generated as name0.vts, name1.vts, etc (not name_0.vts)
            for (size_t p = 0; p < fluidDomain.nPatches(); ++p)
            {
                std::ostringstream vnameWithPatch, pnameWithPatch;
                vnameWithPatch << vname << p << ".vts";
                pnameWithPatch << pname << p << ".vts";
                collectionFlow.addPart(vnameWithPatch.str(), time, "velocity", p);
                collectionFlow.addPart(pnameWithPatch.str(), time, "pressure", p);
            }
            
            // Output mesh and central square
            if (solver.isALEActive())
            {
                // Output the central square
                std::string sqname = "proper_ale_central_square_" + std::to_string(step);
                gsWriteParaview(centralSquareNew, sqname, samp);
            }
            
            // Output the current fluid mesh
            std::string mname = "proper_ale_fluid_mesh_" + std::to_string(step);
            gsWriteParaview(fluidDomain, mname, samp);
        }
    }
    
    // Save paraview collections
    collectionFlow.save();
    collectionMesh.save();
    
    gsInfo << "\nProper ALE simulation completed!\n";
    gsInfo << "Key differences from decay function approach:\n";
    gsInfo << "- Inner boundaries strictly follow rotating square motion\n";
    gsInfo << "- Mesh deformation computed by solving elastic equations\n";
    gsInfo << "- No artificial decay functions used\n";
    
    return 0;
}