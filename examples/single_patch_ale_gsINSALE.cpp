/// Example showing how to connect inner and outer boundaries as single curves
/// and create a single patch parameterization
///
/// Author: J.Li
#include <gismo.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>
#include <gsHLBFGS/gsHLBFGS.h>

using namespace gismo;

// Create a single patch from connected boundary curves
gsMultiPatch<> createConnectedBoundaryPatch(real_t angle = 0.0)
{
    // Define outer boundary curve (square from -1 to 2)
    gsMatrix<> outerCoefs(5, 2);
    outerCoefs << -1, -1,
                    2, -1,
                    2,  2,
                   -1,  2,
                   -1, -1;
    gsKnotVector<> kv(0, 1, 3, 2);
    gsBSpline<> outerCurve(kv, outerCoefs);
    
    // Define inner boundary curve (unit square)
    gsMatrix<> innerCoefs(5, 2);
    innerCoefs << 0, 0,
                   1, 0,
                   1, 1,
                   0, 1,
                   0, 0;
    gsBSpline<> innerCurve(kv, innerCoefs);
    
    // For outer curve, always split at corner (xi = 0.25) for static division
    real_t outerXi = 0.25; // Always split at bottom-left corner
    
    // Split and reorder outer curve at fixed position
    gsBSpline<> left, right;
    outerCurve.splitAt(outerXi, left, right);
    
    // Apply rotation to inner curve
    auto innerCurve2 = innerCurve;
    if (std::abs(angle) > 1e-12)
    {
        innerCurve2.translate(gsEigen::Vector2d(-0.5, -0.5));
        innerCurve2.rotate(angle);
        innerCurve2.translate(-gsEigen::Vector2d(-0.5, -0.5));
    }
    
    // For inner curve, calculate split based on rotation angle
    real_t innerXi = 0.25; // Default at corner
    
    if (std::abs(angle) > 1e-12)
    {
        // The inner curve needs to split at a position that aligns with 
        // the fixed outer split through the rotation center
        // Since outer splits at corner (-1,-1) and rotation center is (0.5,0.5),
        // we need to find where this radial line intersects the rotated inner square
        
        // Convert rotation angle to parameter shift on inner curve
        // Inner curve parameter: 0-0.25 (bottom), 0.25-0.5 (right), 0.5-0.75 (top), 0.75-1 (left)
        innerXi = 0.25 - angle/(2.0*M_PI);
        
        // Normalize to [0,1]
        while (innerXi < 0.0) innerXi += 1.0;
        while (innerXi > 1.0) innerXi -= 1.0;
    }
    
    gsInfo << "Outer split at xi = " << outerXi << ", inner split at xi = " << innerXi << "\n";
    
    // Split inner curve at calculated position
    gsBSpline<> innerLeft, innerRight;
    innerCurve2.splitAt(innerXi, innerLeft, innerRight);
    innerRight.merge(&innerLeft);
    innerRight.knots(0).affineTransformTo(0, 1);
    innerCurve2 = innerRight;
    
    // Merge to create connected curve starting from split point
    right.merge(&left);
    right.knots(0).affineTransformTo(0, 1);
    
    // Ensure compatible knot vectors
    std::vector<real_t> diffKnots;
    right.knots(0).difference(innerCurve2.knots(0), diffKnots);
    innerCurve2.insertKnots(diffKnots.begin(), diffKnots.end());
    innerCurve2.knots(0).difference(right.knots(0), diffKnots);
    right.insertKnots(diffKnots.begin(), diffKnots.end());
    
    // Create surface with outer and inner boundary curves
    gsMatrix<> surfCoefs;
    surfCoefs.resize(innerCurve2.numCoefs() * 2, 2);
    surfCoefs.topRows(innerCurve2.numCoefs()) = right.coefs();
    surfCoefs.bottomRows(innerCurve2.numCoefs()) = innerCurve2.coefs();
    
    gsKnotVector<> uKnot = innerCurve2.knots(0);
    gsKnotVector<> vKnot(0, 1, 0, 2);
    uKnot.affineTransformTo(0, 4);
    
    gsTensorBSpline<2> surf2(uKnot, vKnot, surfCoefs);
    surf2.degreeElevate(1);
    
    // Optimize using barrier method
    gsBarrierPatch<2, real_t> optthisone(surf2, true);
    optthisone.options().setInt("Verbose", 0);
    optthisone.options().setInt("ParamMethod", 1);
    optthisone.compute();
    surf2 = dynamic_cast<gsTensorBSpline<2>&>(optthisone.result().patch(0));
    
    // Return the single optimized patch
    gsMultiPatch<> result;
    result.addPatch(surf2.clone());
    
    gsInfo << "Created single connected boundary patch\n";
    
    return result;
}

// Create multiple patches from the connected boundary
gsMultiPatch<> createMultiPatchFromConnectedBoundary(real_t angle = 0.0)
{
    // First create single patch
    gsMultiPatch<> singlePatch = createConnectedBoundaryPatch(angle);
    gsTensorBSpline<2,real_t>& surf = static_cast<gsTensorBSpline<2,real_t>&>(singlePatch.patch(0));
    
    // Define corner points for splitting (outer boundary corners)
    gsMatrix<> outerCorners(5, 2);
    outerCorners << -1, -1,
                     2, -1,
                     2,  2,
                    -1,  2,
                    -1, -1;
    
    gsMultiPatch<> result;
    gsTensorBSpline<2, real_t> currentPatch = surf;
    
    // Split at corner locations and tag boundaries
    for (index_t i = 1; i < outerCorners.rows() - 1; ++i)
    {
        gsMatrix<> targetPoint = outerCorners.row(i);
        gsMatrix<> paramPoint;
        
        // Find parameter value for corner point on the outer boundary
        currentPatch.invertPoints(targetPoint.transpose(), paramPoint);
        
        gsInfo << "Splitting at u = " << paramPoint(0, 0) << " for corner " << i << "\n";
        
        // Split the patch
        gsTensorBSpline<2, real_t> leftPatch, rightPatch;
        currentPatch.splitAt(0, paramPoint(0, 0), leftPatch, rightPatch);
        
        result.addPatch(leftPatch.clone());
        
        currentPatch = rightPatch;
    }
    
    // Add the last piece
    result.addPatch(currentPatch.clone());
    
    // Compute topology to establish patch connectivity
    result.computeTopology();
    
    // After computing topology, we need to explicitly mark the outer boundaries
    // In our parameterization, v=1 (north) is always the outer physical boundary
    for (index_t p = 0; p < result.nPatches(); ++p)
    {
        // Mark the north side as boundary (not interface)
        patchSide ps(p, boundary::north);
        result.addBoundary(ps);
    }
    
    gsInfo << "Created " << result.nPatches() << " patches\n";
    gsInfo << "Total boundaries: " << result.boundaries().size() << "\n";
    gsInfo << "Total interfaces: " << result.interfaces().size() << "\n";
    
    return result;
}

int main(int argc, char* argv[]){
    gsInfo << "Connected boundary single-patch flow with ALE.\n";

    //=====================================//
                    // Input //
    //=====================================//

    // Problem parameters
    real_t Re = 100;          // Reynolds number
    real_t meanVelocity = 1.0;// Mean inlet velocity
    real_t L = 1.0;           // Characteristic length
    
    // Calculate kinematic viscosity from Reynolds number
    real_t nu = meanVelocity * L / Re;  // nu = 1.0 * 1.0 / 100 = 0.01
    
    // Time parameters
    real_t timeSpan = 0.2;    // Total simulation time
    real_t timeStep = 0.01;   // Time step size
    real_t maxAngle = 20.0;   // Maximum rotation angle in degrees
    index_t nTimeSteps = timeSpan / timeStep;
    
    // Domain and discretization
    index_t numRefine = 3;
    index_t degree = 3;
    index_t outputInterval = 1;  // Output every N steps
    
    gsCmdLine cmd("Connected boundary ALE example.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numRefine);
    cmd.addReal("t","time","Time span",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addReal("a","angle","Maximum rotation angle in degrees",maxAngle);
    cmd.addInt("o","output","Output interval (save every N steps)",outputInterval);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "=== Connected Boundary ALE Example ===\n";
    gsInfo << "Reynolds number: " << Re << "\n";
    gsInfo << "Kinematic viscosity: " << nu << "\n";
    gsInfo << "Time step: " << timeStep << "\n";
    gsInfo << "Time span: " << timeSpan << "\n";
    gsInfo << "Max rotation angle: " << maxAngle << "\n\n";
    
    // Create initial fluid domain: start from a single connected patch,
    // then split along the four outer corners to obtain 4 patches so that
    // inlet/outlet/top/bottom BCs can be applied on distinct sides.
    gsMultiPatch<> fluidDomain = createMultiPatchFromConnectedBoundary(0.0);
    
    // Store original domain for reference
    gsMultiPatch<> originalFluidDomain = fluidDomain;
    
    // Refine mesh
    for (index_t i = 0; i < numRefine; ++i)
    {
        fluidDomain.uniformRefine();
        originalFluidDomain.uniformRefine();
    }
    
    // Create basis for Taylor-Hood elements
    // First create a copy for pressure (keeping it at lower degree)
    gsMultiPatch<> pressureDomain = fluidDomain;
    
    // Elevate degree for velocity field
    fluidDomain.degreeElevate(1);    // fluidDomain is now higher degree
    originalFluidDomain.degreeElevate(1);  // Keep original consistent
    
    // Create basis functions from geometries
    gsMultiBasis<> basisVelocity(fluidDomain);     // Velocity basis: P2
    gsMultiBasis<> basisPressure(pressureDomain);  // Pressure basis: P1
        
    
    // Setup boundary conditions for flow
    gsBoundaryConditions<> bcInfo;
    
    // Velocity and pressure boundary conditions
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    gsConstantFunction<> inletVel(meanVelocity, 0.0, 2);
    gsConstantFunction<> zeroPressure(0.0, 2);
    
    // Note: Inner boundaries (v=0) have NO explicit velocity BCs in ALE.
    // The no-slip condition at the moving inner boundary is enforced implicitly by u-w.
    // Outer boundary is always v=1 (north) per createMultiPatchFromConnectedBoundary.
    
    gsInfo << "Setting boundary conditions on 4 outer segments (per patch north side).\n";
    bool pressureFixed = false;
    const real_t tol = 1e-6;
    for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
    {
        // Determine which physical side this patch's north boundary belongs to
        // Use average of north-side control points to classify
        const gsMatrix<>& C = fluidDomain.patch(p).coefs();
        const index_t nU = basisVelocity.basis(p).component(0).size();
        const index_t nV = basisVelocity.basis(p).component(1).size();
        gsVector<> avg(2); avg.setZero();
        for (index_t i = 0; i < nU; ++i)
        {
            const index_t idx = i + (nV - 1) * nU; // north row
            avg(0) += C(idx, 0);
            avg(1) += C(idx, 1);
        }
        avg /= static_cast<real_t>(nU);
        
        if (std::abs(avg(0) - (-1.0)) < 1e-3) // Left boundary -> inlet
        {
            bcInfo.addCondition(p, boundary::north, condition_type::dirichlet, &inletVel, 0);
            gsInfo << "  Patch " << p << " north classified as LEFT (inlet).\n";
        }
        else if (std::abs(avg(0) - 2.0) < 1e-3) // Right boundary -> outlet (do-nothing)
        {
            // Do not impose velocity; add a single pressure fix here (once)
            if (!pressureFixed)
            {
                bcInfo.addCondition(p, boundary::north, condition_type::dirichlet, &zeroPressure, 1);
                pressureFixed = true;
                gsInfo << "  Patch " << p << " north classified as RIGHT (outlet), pressure fixed.\n";
            }
            else
            {
                gsInfo << "  Patch " << p << " north classified as RIGHT (outlet).\n";
            }
        }
        else if (std::abs(avg(1) - 2.0) < 1e-3) // Top wall
        {
            bcInfo.addCondition(p, boundary::north, condition_type::dirichlet, &zeroVel, 0);
            gsInfo << "  Patch " << p << " north classified as TOP wall.\n";
        }
        else if (std::abs(avg(1) - (-1.0)) < 1e-3) // Bottom wall
        {
            bcInfo.addCondition(p, boundary::north, condition_type::dirichlet, &zeroVel, 0);
            gsInfo << "  Patch " << p << " north classified as BOTTOM wall.\n";
        }
        else
        {
            gsWarn << "  Patch " << p << ": could not classify north boundary reliably (avg x="
                   << avg(0) << ", y=" << avg(1) << ")\n";
        }
    }
    if (!pressureFixed)
    {
        // As a fallback, fix pressure at the first patch's north side
        bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &zeroPressure, 1);
        gsWarn << "Pressure was not fixed on outlet; fixed at patch 0 north as fallback.\n";
    }
    
    // Create PDE
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &zeroVel, nu);
    
    // Debug: Print boundary conditions summary
    gsInfo << "\nBoundary conditions summary:\n";
    gsInfo << "Number of boundary condition containers: " << bcInfo.size() << "\n";
    
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
    
    // Debug: Print solver options to verify they are set correctly
    gsInfo << "\nSolver options:\n";
    gsInfo << "  Time step: " << params.options().getReal("timeStep") << "\n";
    gsInfo << "  Nonlinear tolerance: " << params.options().getReal("nonlin.tol") << "\n";
    gsInfo << "  Linear tolerance: " << params.options().getReal("lin.tol") << "\n";
    
    // Create ALE solver
    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));
    
    // Initialize solver
    solver.initialize();
    
    // Debug: Check solver state after initialization
    gsInfo << "\nSolver state after initialization:\n";
    gsInfo << "  Is ALE initially active: " << solver.isALEActive() << "\n";
    
    // Activate ALE after initialization
    solver.setALEActive(true);
    gsInfo << "  Is ALE active after setting: " << solver.isALEActive() << "\n";
    
    // Enable mesh optimization using gsBarrierPatch
    solver.setMeshOptimization(true );
    solver.getMeshOptOptions().setInt("Verbose", 1);
    solver.getMeshOptOptions().setInt("ParamMethod", 1);
    solver.getMeshOptOptions().setInt("AAPreconditionType", 0);
    
    // Define mesh motion function for ALE
    // This function returns the CUMULATIVE displacement from the original mesh
    auto meshMotion = [&](real_t t) -> gsMatrix<real_t>
    {
        // Calculate rotation angle with sinusoidal motion
        const real_t angle = maxAngle * std::sin(2 * M_PI * t / timeSpan) * M_PI / 180.0;
        const real_t cos_a = std::cos(angle);
        const real_t sin_a = std::sin(angle);
        const real_t cx = 0.5;  // rotation center x
        const real_t cy = 0.5;  // rotation center y

        const index_t udofs = solver.getAssembler()->getUdofs();
        gsMatrix<> disp(2 * udofs, 1); 
        disp.setZero();

        const gsDofMapper& mapper = solver.getAssembler()->getMappers()[0];
        
        // Loop all patches; inner boundary is always v=0 (south)
        for (index_t p = 0; p < originalFluidDomain.nPatches(); ++p)
        {
            const gsMatrix<>& originalCoefs = originalFluidDomain.patch(p).coefs();
            const index_t nU = basisVelocity.basis(p).component(0).size();
            const index_t nV = basisVelocity.basis(p).component(1).size();

            if (t == timeStep)
            {
                gsInfo << "\n=== Mesh Motion Debug (patch " << p << ") ===\n";
                gsInfo << "Time: " << t << ", Angle: " << angle * 180.0 / M_PI << " deg, nU=" << nU << ", nV=" << nV << "\n";
            }

            // Move inner boundary control points (south row)
            for (index_t i = 0; i < nU; ++i)
            {
                const index_t cpIndex = i; // v=0 row
                const real_t x = originalCoefs(cpIndex, 0);
                const real_t y = originalCoefs(cpIndex, 1);

                if (!mapper.is_free(cpIndex, p))
                    continue; // skip constrained velocity DOFs

                // Rotate around (cx,cy)
                const real_t xr = cos_a * (x - cx) - sin_a * (y - cy) + cx;
                const real_t yr = sin_a * (x - cx) + cos_a * (y - cy) + cy;

                const real_t dx = xr - x;
                const real_t dy = yr - y;

                const index_t base = mapper.index(cpIndex, p);
                if (base >= 0 && base < udofs)
                {
                    disp(base) = dx;
                    disp(base + udofs) = dy;
                }
            }
        }
        
        // Count non-zero displacements
        if (t == timeStep)
        {
            index_t nonZeroCount = 0;
            for (index_t i = 0; i < disp.rows(); ++i)
            {
                if (std::abs(disp(i, 0)) > 1e-10)
                    nonZeroCount++;
            }
            gsInfo << "\nTotal non-zero displacements: " << nonZeroCount << " out of " << disp.rows() << "\n";
            gsInfo << "Max displacement magnitude: " << disp.lpNorm<gsEigen::Infinity>() << "\n";
            gsInfo << "==============================\n\n";
        }
        
        return disp;
    };
    
    solver.setMeshUpdateFunction(meshMotion);
    
    // Output initial configuration
    gsInfo << "Writing initial configuration...\n";
    gsWriteParaview(fluidDomain, "ale_single_initial_mesh", 1000, true);
    
    // Solve Stokes problem for initial conditions
    gsInfo << "Solving initial Stokes problem...\n";
    solver.solveStokes();
    
    // Time stepping loop
    gsInfo << "Starting time integration with ALE...\n";
    
    // Prepare paraview collection for time series output
    gsParaviewCollection collectionVelocity("ale_single_velocity");
    gsParaviewCollection collectionPressure("ale_single_pressure");
    gsParaviewCollection collectionMesh("ale_single_mesh");
    
    // Update recalculation of nTimeSteps after command line parsing
    nTimeSteps = timeSpan / timeStep;
    
    for (index_t step = 0; step <= nTimeSteps; ++step)
    {
        real_t time = step * timeStep;
        gsInfo << "\nTime step " << step << ", t = " << time << "\n";
        
        // Solve fluid problem with ALE
        if (step > 0)
        {
            gsInfo << "  Calling nextIteration() for time integration...\n";
            solver.nextIteration();
            gsInfo << "  nextIteration() completed\n";
        }
        
        // Get solution
        gsMatrix<> velCoefs = solver.solutionCoefs(0);
        gsMatrix<> presCoefs = solver.solutionCoefs(1);
        
        // Output some statistics
        gsInfo << "  Max velocity: " << velCoefs.lpNorm<gsEigen::Infinity>() << "\n";
        
        // Proper pressure monitoring: show range instead of max absolute value
        real_t pmin = presCoefs.minCoeff();
        real_t pmax = presCoefs.maxCoeff();
        real_t prange = pmax - pmin;
        gsInfo << "  Pressure range: p_min = " << pmin << ", p_max = " << pmax << ", Δp = " << prange << "\n";
        
        if (solver.isALEActive() && step > 0)
        {
            gsField<> meshDispField = solver.getMeshDisplacementField();
            gsInfo << "  ‖meshDisp‖_∞: "
                   << meshDispField.coefficientVector().template lpNorm<gsEigen::Infinity>() << "\n";
        }
        
        // Export solution periodically
        if (step % outputInterval == 0)
        {
            // Create multipatch for velocity and pressure from coefficients
            gsMultiPatch<> velPatches;
            gsMultiPatch<> presPatches;
            
            // Get the current deformed mesh
            const gsMultiPatch<>& currentFluidDomain = solver.getAssembler()->getPatches();
            
            // Update coefficients for all patches
            const gsDofMapper& uMapper = solver.getAssembler()->getMappers()[0];
            const gsDofMapper& pMapper = solver.getAssembler()->getMappers()[1];
            const index_t shift = solver.getAssembler()->getUdofs();
            
            // Process each patch
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
                        velPatchCoefs(i, 0) = velCoefs(idx);
                        velPatchCoefs(i, 1) = velCoefs(idx + shift);
                    }
                    else
                    {
                        // Handle boundary conditions for constrained DOFs
                        // Get the boundary condition value from the assembler
                        gsMatrix<> fixedValue(2, 1);
                        fixedValue.setZero();
                        
                        // Check if this DOF is on a boundary with Dirichlet conditions
                        index_t nU = basisVelocity.basis(p).component(0).size();
                        index_t nV = basisVelocity.basis(p).component(1).size();
                        
                        // Determine which boundary is outer (similar check as in BC setup)
                        bool isOuterBoundary = false;
                        gsMatrix<> param(2, 1);
                        
                        if (i >= nU * (nV - 1)) // North boundary (v=1)
                        {
                            // Check if north is outer by comparing sizes
                            // (In practice, we know from setup that north is outer if patch extends from -1 to 2)
                            param(0, 0) = 0.5;
                            param(1, 0) = 1.0;
                            gsMatrix<> testPt = currentFluidDomain.patch(p).eval(param);
                            if (std::abs(testPt(0, 0)) > 1.5 || std::abs(testPt(1, 0)) > 1.5)
                            {
                                isOuterBoundary = true;
                                param(0, 0) = real_t(i % nU) / (nU - 1);
                                param(1, 0) = 1.0;
                            }
                        }
                        else if (i < nU) // South boundary (v=0)
                        {
                            // Check if south is outer
                            param(0, 0) = 0.5;
                            param(1, 0) = 0.0;
                            gsMatrix<> testPt = currentFluidDomain.patch(p).eval(param);
                            if (std::abs(testPt(0, 0)) > 1.5 || std::abs(testPt(1, 0)) > 1.5)
                            {
                                isOuterBoundary = true;
                                param(0, 0) = real_t(i % nU) / (nU - 1);
                                param(1, 0) = 0.0;
                            }
                        }
                        
                        if (isOuterBoundary)
                        {
                            // Get physical location to determine boundary type
                            gsMatrix<> physPt = currentFluidDomain.patch(p).eval(param);
                            real_t x = physPt(0, 0);
                            
                            if (std::abs(x - (-1.0)) < 0.1)
                            {
                                // Left boundary - inlet
                                velPatchCoefs(i, 0) = meanVelocity;
                                velPatchCoefs(i, 1) = 0.0;
                            }
                            else
                            {
                                // Other boundaries - walls or outlet
                                velPatchCoefs(i, 0) = 0.0;
                                velPatchCoefs(i, 1) = 0.0;
                            }
                        }
                        else
                        {
                            // Inner boundary (v=0) or other constrained DOFs - zero velocity
                            velPatchCoefs(i, 0) = 0.0;
                            velPatchCoefs(i, 1) = 0.0;
                        }
                    }
                }
                
                // Create velocity patch
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
                        presPatchCoefs(i, 0) = presCoefs(idx);
                    }
                    else
                    {
                        presPatchCoefs(i, 0) = 0.0;
                    }
                }
                
                // Create pressure patch
                auto presPatch = basisPressure.basis(p).makeGeometry(give(presPatchCoefs));
                presPatches.addPatch(give(presPatch));
            }
            
            // Create fields
            gsField<> velocity(currentFluidDomain, velPatches);
            gsField<> pressure(currentFluidDomain, presPatches);
            
            // Output
            std::string vname = "ale_single_velocity_" + std::to_string(step);
            std::string pname = "ale_single_pressure_" + std::to_string(step);
            std::string mname = "ale_single_mesh_" + std::to_string(step);
            
            gsWriteParaview(velocity, vname, 1000, false);
            gsWriteParaview(pressure, pname, 1000, false);
            gsWriteParaview(currentFluidDomain, mname, 1000, false);
            
            // Add all patches to collections
            for (size_t p = 0; p < currentFluidDomain.nPatches(); ++p)
            {
                collectionVelocity.addPart(vname + std::to_string(p) + ".vts", time, std::to_string(p));
                collectionPressure.addPart(pname + std::to_string(p) + ".vts", time, std::to_string(p));
                collectionMesh.addPart(mname + std::to_string(p) + ".vts", time, std::to_string(p));
            }
        }
    }
    
    // Save paraview collections
    collectionVelocity.save();
    collectionPressure.save();
    collectionMesh.save();
    
    gsInfo << "Simulation completed successfully!\n";
    
    return 0;
}
