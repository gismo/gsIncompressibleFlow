/** @file ale_simple_test.cpp
    @brief Simple test of ALE formulation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    // Simple 2D channel flow
    real_t L = 2.0;     // Channel length
    real_t H = 1.0;     // Channel height
    real_t Re = 10;     // Low Reynolds number for stability
    real_t rho = 1.0;
    real_t mu = rho * L / Re;
    real_t dt = 0.1;    // Large time step for testing
    index_t nSteps = 10;
    
    gsInfo << "=== Simple ALE Test ===\n";
    gsInfo << "Reynolds number: " << Re << "\n";
    gsInfo << "Time step: " << dt << "\n\n";
    
    // Create simple rectangular domain
    gsMultiPatch<> domain;
    domain.addPatch(gsNurbsCreator<>::BSplineRectangle(0, 0, L, H));
    
    // Simple discretization
    gsMultiBasis<> basisVel(domain);
    gsMultiBasis<> basisPres(domain);
    
    // Use Taylor-Hood elements: P2 for velocity, P1 for pressure
    basisVel.setDegree(2);
    basisPres.setDegree(1);
    
    // Minimal refinement to avoid too many DOFs
    basisVel.uniformRefine();
    basisPres.uniformRefine();
    
    // Boundary conditions
    gsBoundaryConditions<> bcInfo;
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    gsConstantFunction<> inletVel(1.0, 0.0, 2);
    
    // Velocity BCs - all boundaries have Dirichlet conditions
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &inletVel, 0);   // Inlet
    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &zeroVel, 0);    // Outlet  
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &zeroVel, 0);  // Bottom wall
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &zeroVel, 0);  // Top wall
    
    // For pressure, we need natural boundary conditions
    // The outlet will have zero pressure (natural BC)
    // Note: In incompressible flow, only velocity gradients matter, 
    // so we don't need to explicitly set pressure BCs
    
    // Create PDE
    gsNavStokesPde<real_t> nsPde(domain, bcInfo, &zeroVel, mu);
    
    // Setup bases  
    std::vector<gsMultiBasis<>> bases;
    bases.push_back(basisVel);   // Velocity
    bases.push_back(basisPres);  // Pressure
    
    // Create parameters
    gsFlowSolverParams<real_t> params(nsPde, bases);
    params.options().setReal("timeStep", dt);
    params.options().setInt("nonlin.maxIt", 3);
    params.options().setReal("nonlin.tol", 1e-4);
    params.options().setSwitch("quiet", false);
    
    // Use iterative solver instead of direct solver to handle potential singularity
    params.options().setString("lin.solver", "iter");
    params.options().setString("lin.krylov", "BiCGSTAB");
    params.options().setInt("lin.maxIt", 1000);
    params.options().setReal("lin.tol", 1e-8);
    params.options().setString("lin.precType", "MSIMPLER_FdiagEqual");
    
    // Create ALE solver
    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));
    
    // Initialize
    gsInfo << "Initializing solver...\n";
    solver.initialize();
    
    // Solve initial Stokes
    gsInfo << "Solving initial Stokes problem...\n";
    try {
        solver.solveStokes();
        gsInfo << "Stokes solve completed successfully\n";
    } catch (const std::exception& e) {
        gsInfo << "Error in solveStokes: " << e.what() << "\n";
        return 1;
    }
    
    // Test without ALE first
    gsInfo << "\n--- Testing without ALE ---\n";
    solver.setALEActive(false);
    
    for (index_t step = 0; step < 3; ++step)
    {
        gsInfo << "Time step " << step << "\n";
        solver.nextIteration();
    }
    
    // Now test with ALE
    gsInfo << "\n--- Testing with ALE ---\n";
    solver.setALEActive(true);
    
    // Get free DOFs for mesh displacement (after Dirichlet elimination)
    index_t udofs = solver.getAssembler()->getUdofs();
    index_t totalDofs = udofs * 2;  // 2D problem
    
    gsInfo << "Free velocity DOFs: " << udofs << " (total: " << totalDofs << ")\n";
    
    // Create zero mesh displacement for free DOFs only
    gsMatrix<> meshDisp(totalDofs, 1);
    meshDisp.setZero();
    
    // Small test displacement in x-direction only for free DOFs
    for (index_t i = 0; i < udofs; ++i)
    {
        meshDisp(i, 0) = 0.01;  // Small uniform displacement in x
    }
    
    // Update mesh
    solver.getALEAssembler()->updateMesh(meshDisp);
    
    // Run a few more steps with ALE
    for (index_t step = 3; step < 6; ++step)
    {
        gsInfo << "Time step " << step << " (with ALE)\n";
        solver.nextIteration();
    }
    
    gsInfo << "\nTest completed successfully!\n";
    
    return 0;
}