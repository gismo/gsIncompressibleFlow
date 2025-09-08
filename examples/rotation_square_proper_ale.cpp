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
    real_t Re = 200;          // Reynolds number
    real_t meanVelocity = 20.0;// Mean inlet velocity
    real_t L = 1.0;           // Characteristic length
    
    // Calculate kinematic viscosity from Reynolds number
    // Re = U*L/nu => nu = U*L/Re
    real_t nu = meanVelocity * L / Re;  // nu = 1.0 * 1.0 / 100 = 0.01
    
    // Time parameters
    real_t timeSpan = 3.0;    // Total simulation time
    real_t timeStep = 0.001;   // Time step size
    real_t maxAngle = 5;      // Maximum rotation angle in degrees
    index_t nTimeSteps = timeSpan / timeStep;
    
    // Domain and discretization
    index_t numRefine = 3;
    index_t degree = 3;
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
    gsInfo << "Kinematic viscosity: " << nu << "\n";
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
    // real_t L = 1.0; // Size of each patch - already defined above as characteristic length
    
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
    
    // ================== New basis construction for Taylor-Hood elements ==================
    // We use a simpler and more robust approach to create P2/P1 Taylor-Hood elements
    
    // 1. Create a copy of the domain for pressure (keeping it at degree 1)
    //    At this point fluidDomain is still degree 1 after refinement
    gsMultiPatch<> pressureDomain = fluidDomain;
    
    // 2. Elevate the main domain's degree for velocity field
    fluidDomain.degreeElevate(1);    // fluidDomain is now degree 2
    centralSquare.degreeElevate(1);  // Keep visualization geometry consistent
    
    // 3. Create basis functions from geometries of different degrees
    gsMultiBasis<> basisVelocity(fluidDomain);     // Velocity basis: P2 (from degree 2 geometry)
    gsMultiBasis<> basisPressure(pressureDomain);  // Pressure basis: P1 (from degree 1 geometry)
    
    // ==========================================================
    
    // Setup boundary conditions for flow
    gsBoundaryConditions<> bcInfo;
    
    // Velocity boundary conditions
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    gsConstantFunction<> inletVel(meanVelocity, 0.0, 2);
    
    // Inlet (left side) - set velocity boundary condition
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
    
    // Inner boundaries (interfaces with rotating square) - ALE automatic no-slip
    // In ALE formulation, when mesh motion is properly defined, the no-slip condition u = w
    // is automatically satisfied at moving boundaries. The mesh motion function defines w,
    // and the ALE solver ensures u = w at the fluid-structure interface.
    
    // Outlet (right side) - natural boundary condition (do-nothing)
    
    // Add pressure reference point to avoid floating pressure
    gsConstantFunction<> zeroPressure(0.0, 2);
    bcInfo.addCondition(5, boundary::east, condition_type::dirichlet, &zeroPressure, 1); // patch 5, east side, pressure
    
    // Create PDE
    // Use the kinematic viscosity calculated from Reynolds number
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &zeroVel, nu);
    
    // Setup solver parameters
    std::vector<gsMultiBasis<>> discreteBases;
    discreteBases.push_back(basisVelocity);  // Velocity basis
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
    
    // Artificial mesh motion
    auto meshMotion = [&](real_t t) -> gsMatrix<real_t>
    {
        // Calculate rotation angle with your specified motion
        // Continuously rotate to the right (clockwise) with constant angular velocity
        const real_t angle = -maxAngle * t / timeSpan * M_PI / 180.0; // Convert to radians, negative for clockwise
        const real_t cos_a = std::cos(angle);
        const real_t sin_a = std::sin(angle);
        const real_t cx = 0.5;  // rotation center x (center of the central square)
        const real_t cy = 0.5;  // rotation center y (center of the central square)

        const index_t udofs = solver.getAssembler()->getUdofs();
        gsMatrix<> disp(2 * udofs, 1); 
        disp.setZero();

        // Get the velocity DOF mapper to iterate through control points
        const gsDofMapper& mapper = solver.getAssembler()->getMappers()[0];
        
        // Debug counters and tracking
        static int call_count = 0;
        call_count++;
        int boundary_points = 0;
        real_t max_disp = 0.0;
        
        // Traverse all patches and all control points (DOFs)
        for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
        {
            const gsMatrix<>& originalCoefs = originalFluidDomain.patch(p).coefs();
            for (index_t i = 0; i < originalCoefs.rows(); ++i)
            {
                // Only process free DOFs (not constrained by boundary conditions)
                if (!mapper.is_free(i, p))
                    continue;

                const real_t x = originalCoefs(i, 0);
                const real_t y = originalCoefs(i, 1);
                const real_t tol = 1e-9;

                // Here set the boundary condition for boundary points
                bool isOnInnerBoundary = 
                    ( (std::abs(x - 0.0) < tol || std::abs(x - 1.0) < tol) && (y >= 0.0 - tol && y <= 1.0 + tol) ) ||
                    ( (std::abs(y - 0.0) < tol || std::abs(y - 1.0) < tol) && (x >= 0.0 - tol && x <= 1.0 + tol) );

                if (isOnInnerBoundary)
                {
                    boundary_points++;
                    
                    // Apply rotation transformation around the center
                    real_t xr = cos_a * (x - cx) - sin_a * (y - cy) + cx;
                    real_t yr = sin_a * (x - cx) + cos_a * (y - cy) + cy;
                    
                    // Calculate displacement from original to rotated position
                    real_t dx = xr - x;
                    real_t dy = yr - y;
                    
                    real_t disp_mag = std::sqrt(dx*dx + dy*dy);
                    max_disp = std::max(max_disp, disp_mag);
                    
                    // Map to global DOF indices and set displacement
                    index_t base = mapper.index(i, p);
                    if (base < udofs && base + udofs < disp.size())
                    {
                        disp(base) = dx;           // x-displacement
                        disp(base + udofs) = dy;   // y-displacement
                        
                        // Debug first few calls
                        if (call_count <= 2 && boundary_points <= 3)
                        {
                            gsInfo << "    Point (" << x << "," << y << ") -> (" << xr << "," << yr 
                                   << ") disp=(" << dx << "," << dy << ") mag=" << disp_mag 
                                   << " idx=" << base << "\n";
                        }
                    }
                    else
                    {
                        gsInfo << "WARNING: Invalid DOF index " << base << " (udofs=" << udofs 
                               << ", disp.size=" << disp.size() << ")\n";
                    }
                }
                // For control points not on inner boundary, displacement remains zero
                // This allows the mesh optimization to handle the deformation smoothly
            }
        }
        
        // Debug output for first few calls
        if (call_count <= 3)
        {
            gsInfo << "meshMotion call " << call_count << ": t=" << t 
                   << ", angle=" << angle << " rad, boundary_points=" << boundary_points 
                   << ", max_disp=" << max_disp << "\n";
        }
        
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
    gsParaviewCollection collectionVelocity("proper_ale_velocity_animation");
    gsParaviewCollection collectionPressure("proper_ale_pressure_animation");
    gsParaviewCollection collectionMesh("proper_ale_rotating_square_mesh");
    
    // Update recalculation of nTimeSteps after command line parsing
    nTimeSteps = timeSpan / timeStep;
    
    for (index_t step = 0; step <= nTimeSteps; ++step)
    {
        real_t time = step * timeStep;
        gsInfo << "\nTime step " << step << ", t = " << time << "\n";
        
        // Update central square position for visualization
        // Same angle formula for visualization (in degrees)
        real_t angle = -maxAngle * time / timeSpan; // negative for clockwise
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
        gsInfo << "  solver time = " << solver.getSimulationTime() << "\n";
        
        if (solver.isALEActive() && step > 0)
        {
            gsField<> meshDispField = solver.getMeshDisplacementField();
            gsField<> meshVelField = solver.getMeshVelocityField();
            
            real_t meshDispNorm = meshDispField.coefficientVector().template lpNorm<gsEigen::Infinity>();
            real_t meshVelNorm = meshVelField.coefficientVector().template lpNorm<gsEigen::Infinity>();
            
            gsInfo << "  ||meshDisp||_inf: " << meshDispNorm << "\n";
            gsInfo << "  ||meshVel||_inf: " << meshVelNorm << "\n";
            
            // =================== VERIFICATION: u vs w independence ===================
            if (step <= 3 || step % 50 == 0) // Verify at early steps and periodically
            {
                gsInfo << "\n  === VERIFYING u vs w INDEPENDENCE ===\n";
                
                // Get fluid velocity coefficients  
                gsMatrix<> fluidVelCoefs = solver.solutionCoefs(0);
                gsMatrix<> meshVelCoefs = meshVelField.coefficientVector();
                
                // Get DOF mappers
                const gsDofMapper& uMapper = solver.getAssembler()->getMappers()[0];
                const index_t udofs = solver.getAssembler()->getUdofs();
                
                gsInfo << "  Fluid velocity DOFs: " << fluidVelCoefs.rows() << "\n";
                gsInfo << "  Mesh velocity DOFs: " << meshVelCoefs.rows() << "\n";
                gsInfo << "  udofs: " << udofs << "\n";
                
                // Sample some points to compare u and w
                int boundary_samples = 0;
                int interior_samples = 0;
                real_t boundary_diff_sum = 0.0;
                real_t interior_diff_sum = 0.0;
                
                // Check a few patches in the middle region
                std::vector<index_t> test_patches = {1, 3, 4, 6}; // patches around central square
                
                for (index_t p : test_patches)
                {
                    const gsMatrix<>& originalCoefs = originalFluidDomain.patch(p).coefs();
                    gsInfo << "  Checking patch " << p << " with " << originalCoefs.rows() << " control points\n";
                    
                    for (index_t i = 0; i < std::min(originalCoefs.rows(), (index_t)10); ++i) // Check first 10 points per patch
                    {
                        if (!uMapper.is_free(i, p)) continue;
                        
                        const real_t x = originalCoefs(i, 0);
                        const real_t y = originalCoefs(i, 1);
                        const real_t tol = 1e-9;
                        
                        // Check if this point is on boundary or interior
                        bool isOnInnerBoundary = 
                            ( (std::abs(x - 0.0) < tol || std::abs(x - 1.0) < tol) && (y >= 0.0 - tol && y <= 1.0 + tol) ) ||
                            ( (std::abs(y - 0.0) < tol || std::abs(y - 1.0) < tol) && (x >= 0.0 - tol && x <= 1.0 + tol) );
                        
                        bool isOnOuterBoundary = 
                            (p == 0 && (std::abs(x + 1.0) < tol || std::abs(y + 1.0) < tol)) ||
                            (p == 2 && (std::abs(x + 1.0) < tol || std::abs(y - 2.0) < tol)) ||
                            (p == 5 && (std::abs(x - 2.0) < tol || std::abs(y + 1.0) < tol)) ||
                            (p == 7 && (std::abs(x - 2.0) < tol || std::abs(y - 2.0) < tol));
                        
                        bool isBoundary = isOnInnerBoundary || isOnOuterBoundary;
                        
                        // Get DOF index
                        index_t base = uMapper.index(i, p);
                        if (base >= udofs || base + udofs >= fluidVelCoefs.rows()) continue;
                        if (base >= meshVelCoefs.rows()/2) continue;
                        
                        // Get velocity components
                        real_t u_x = fluidVelCoefs(base);
                        real_t u_y = fluidVelCoefs(base + udofs);
                        real_t w_x = meshVelCoefs(2*base);
                        real_t w_y = meshVelCoefs(2*base + 1);
                        
                        real_t diff_mag = std::sqrt((u_x - w_x)*(u_x - w_x) + (u_y - w_y)*(u_y - w_y));
                        
                        if (isBoundary)
                        {
                            boundary_samples++;
                            boundary_diff_sum += diff_mag;
                            if (boundary_samples <= 5) // Show first few boundary points
                            {
                                gsInfo << "    BOUNDARY Point (" << x << "," << y << "): u=(" << u_x << "," << u_y 
                                       << ") w=(" << w_x << "," << w_y << ") |u-w|=" << diff_mag << "\n";
                            }
                        }
                        else
                        {
                            interior_samples++;
                            interior_diff_sum += diff_mag;
                            if (interior_samples <= 5) // Show first few interior points
                            {
                                gsInfo << "    INTERIOR Point (" << x << "," << y << "): u=(" << u_x << "," << u_y 
                                       << ") w=(" << w_x << "," << w_y << ") |u-w|=" << diff_mag << "\n";
                            }
                        }
                    }
                }
                
            }
            
                   
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
                gsWriteParaview(currentFluidDomain, meshname, 1000, false);
                gsWriteParaview(meshDispField, dispname, 1000, false);
                
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
            // Note: We keep the deformed geometry for visualization
            fluidDomain = currentFluidDomain;
        }
        
        // Export solution periodically
        if (step % outputInterval == 0)
        {
            // Get velocity and pressure coefficient vectors from already computed solution
            gsMatrix<> velCoefsFull = velCoefs;
            gsMatrix<> presCoefsFull = presCoefs;
            
            // Create multipatch for velocity (2D vector field)
            // Note: velocity uses P2 basis, pressure uses P1
            gsMultiPatch<> velPatches;
            gsMultiPatch<> presPatches;
            
            // Update coefficients for each patch
            const gsDofMapper& uMapper = solver.getAssembler()->getMappers()[0];
            const gsDofMapper& pMapper = solver.getAssembler()->getMappers()[1];
            const index_t shift = solver.getAssembler()->getUdofs();
            
            // Debug: print number of patches
            if (step == 0)
                gsInfo << "Processing " << fluidDomain.nPatches() << " patches for output\n";
                
            for (size_t p = 0; p < fluidDomain.nPatches(); ++p)
            {
                // Velocity field (2 components)
                index_t nCoefs = basisVelocity.basis(p).size();
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
                        // For constrained DOFs, we need to determine the actual boundary condition
                        // Get the anchor point of this basis function in the parametric domain
                        gsMatrix<> anchors = basisVelocity.basis(p).anchors();
                        
                        // Check if this is a boundary DOF
                        real_t u = anchors(0, i % anchors.cols());
                        real_t v = anchors(1, i / anchors.cols());
                        
                        // Check which boundary this DOF is on
                        const real_t tol = 1e-10;
                        bool isWest = (std::abs(u - 0.0) < tol);
                        bool isEast = (std::abs(u - 1.0) < tol);
                        bool isSouth = (std::abs(v - 0.0) < tol);
                        bool isNorth = (std::abs(v - 1.0) < tol);
                        
                        // For patches 0, 1, 2, west boundary is inlet
                        if (p <= 2 && isWest)
                        {
                            velPatchCoefs(i, 0) = meanVelocity;
                            velPatchCoefs(i, 1) = 0.0;
                        }
                        else
                        {
                            // All other boundaries are walls (zero velocity)
                            velPatchCoefs(i, 0) = 0.0;
                            velPatchCoefs(i, 1) = 0.0;
                        }
                    }
                }

                // Create geometry using the basis and coefficients
                auto velPatch = basisVelocity.basis(p).makeGeometry(give(velPatchCoefs));
                velPatches.addPatch(give(velPatch));
                
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
                        presPatchCoefs(i, 0) = 0.0; // TODO: implement proper extrapolation
                    }
                }
                
                // Create geometry using the basis and coefficients
                auto presPatch = basisPressure.basis(p).makeGeometry(give(presPatchCoefs));
                presPatches.addPatch(give(presPatch));
            }
            
            // Create fields
            if (step == 0)
            {
                gsInfo << "Created velocity patches: " << velPatches.nPatches() << "\n";
                gsInfo << "Created pressure patches: " << presPatches.nPatches() << "\n";
            }
            gsField<> velocity(fluidDomain, velPatches);
            gsField<> pressure(fluidDomain, presPatches);
            
            // Output flow fields
            unsigned samp = 1000; // sampling points per direction
            std::string vname = "proper_ale_velocity_" + std::to_string(step);
            std::string pname = "proper_ale_pressure_" + std::to_string(step);
            
            // Write velocity and pressure fields
            // For animation, we need to write without creating individual PVD files
            gsWriteParaview(velocity, vname, samp, false);
            gsWriteParaview(pressure, pname, samp, false);
            
            // Add all patches to the time series collection
            // Each patch is written as a separate VTS file
            for (size_t p = 0; p < fluidDomain.nPatches(); ++p)
            {
                // The files are named as: name0.vts, name1.vts, etc.
                std::string vfile = vname + std::to_string(p) + ".vts";
                std::string pfile = pname + std::to_string(p) + ".vts";
                
                // Add each patch as a separate part
                collectionVelocity.addPart(vfile, time, std::to_string(p));
                collectionPressure.addPart(pfile, time, std::to_string(p));
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
    collectionVelocity.save();
    collectionPressure.save();
    collectionMesh.save();
    
    
    return 0;
}