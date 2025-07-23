/** @file ale_fsi_example.cpp
    @brief Example of using ALE formulation for FSI simulation

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
    // Problem parameters
    real_t Re = 100;          // Reynolds number
    real_t rho = 1.0;         // Fluid density
    real_t mu = rho / Re;     // Dynamic viscosity
    real_t U = 1.0;           // Characteristic velocity
    real_t L = 1.0;           // Characteristic length
    
    // Time parameters
    real_t T = 10.0;          // Total simulation time
    real_t dt = 0.01;         // Time step size
    index_t nTimeSteps = T / dt;
    
    // Domain and discretization
    index_t numRefine = 2;
    index_t degree = 2;
    
    gsInfo << "=== ALE FSI Example ===\n";
    gsInfo << "Reynolds number: " << Re << "\n";
    gsInfo << "Time step: " << dt << "\n";
    gsInfo << "Number of time steps: " << nTimeSteps << "\n\n";
    
    // Create fluid domain (unit square)
    gsMultiPatch<> fluidDomain;
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1));
    
    // Refine mesh
    for (index_t i = 0; i < numRefine; ++i)
        fluidDomain.uniformRefine();
    
    // Create basis
    gsMultiBasis<> basis(fluidDomain);
    basis.setDegree(degree);
    
    // Setup boundary conditions
    gsBoundaryConditions<> bcInfo;
    
    // Velocity boundary conditions
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    gsConstantFunction<> inletVel(U, 0.0, 2);
    
    // Fixed walls (top and bottom)
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &zeroVel, 0);
    
    // Inlet (left) - now also zero velocity
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &zeroVel, 0);
    
    // FSI interface (right) - will be updated dynamically
    // Leave it as natural BC for now
    
    // Pressure boundary conditions
    // Natural BC at outlet/interface
    
    // Create PDE
    // Create Navier-Stokes PDE
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &zeroVel, mu);
    
    // Setup solver parameters
    std::vector<gsMultiBasis<>> discreteBases;
    discreteBases.push_back(basis);  // Velocity basis
    discreteBases.push_back(basis);  // Pressure basis
    
    // For Taylor-Hood elements, elevate velocity space
    discreteBases[0].degreeElevate(1);
    
    gsFlowSolverParams<real_t> params(nsPde, discreteBases);
    
    // Solver options
    params.options().setInt("nonlin.maxIt", 100);     // Max Picard iterations
    params.options().setReal("nonlin.tol", 1e-6);     // Picard tolerance
    params.options().setInt("lin.maxIt", 200);        // Max linear solver iterations
    params.options().setReal("lin.tol", 1e-8);        // Linear solver tolerance
    params.options().setReal("timeStep", dt);         // Time step size
    params.options().setSwitch("quiet", false);    // Enable verbose output
    
    // Create ALE solver
    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));
    
    // Initialize solver
    solver.initialize();
    
    // Activate ALE after initialization
    solver.setALEActive(true);
    
    // Define mesh motion (example: oscillating boundary)
    // 圆周运动：R=0.05, f=1 Hz
    auto meshMotion = [&](real_t t)
    {
        const real_t R = 0.05, w = 2*M_PI;
        const real_t dx = R*std::cos(w*t);
        const real_t dy = R*std::sin(w*t);

        const index_t udofs = solver.getAssembler()->getUdofs();
        gsMatrix<> disp(2*udofs,1); disp.setZero();

        // 控制点物理坐标
        const gsMatrix<>& C = fluidDomain.patch(0).coefs();
        const gsDofMapper& mapper = solver.getAssembler()->getMappers()[0];

        for (index_t i=0; i<C.rows(); ++i)
        {
            real_t x = C(i,0), y = C(i,1);

            // 边界控制点：x==0/1 或 y==0/1，保持 0
            if (x==0 || x==1 || y==0 || y==1) continue;

            if (mapper.is_free(i,0))
            {
                index_t base = mapper.index(i,0);   // x-分量索引
                disp(base      ) = dx;
                disp(base+udofs) = dy;
            }
        }
        return disp;
    };
    solver.setMeshUpdateFunction(meshMotion);
    
    // Solve Stokes problem for initial conditions
    gsInfo << "Solving initial Stokes problem...\n";
    solver.solveStokes();
    
    // Time stepping loop
    gsInfo << "Starting time integration...\n";
    
    for (index_t step = 0; step < nTimeSteps; ++step)
    {
        real_t time = step * dt;
        gsInfo << "\nTime step " << step << ", t = " << time << "\n";
        
        // Solve fluid problem with ALE
        solver.nextIteration();
        
        // Get solution
        gsMatrix<> velCoefs = solver.solutionCoefs(0);
        gsMatrix<> presCoefs = solver.solutionCoefs(1);
        gsField<> meshDispField = solver.getMeshDisplacementField();
        gsInfo << "  ‖meshDisp‖_∞: "
               << meshDispField.coefficientVector().template lpNorm<gsEigen::Infinity>() << "\n";
        
        // Output some statistics
        gsInfo << "  Max velocity: " << velCoefs.lpNorm<gsEigen::Infinity>() << "\n";
        gsInfo << "  Max pressure: " << presCoefs.lpNorm<gsEigen::Infinity>() << "\n";
        
        // Export solution periodically
        if (step % 10 == 0)
        {
            // Create fields
            // Construct solution fields from coefficients
            gsField<> velocity = solver.constructSolution(0);  // 0 for velocity
            gsField<> pressure = solver.constructSolution(1);  // 1 for pressure
            
            // 当 ALE 激活时，输出网格位移
            unsigned samp = 1000; // 每方向采样点数
            std::string vname = "fluid_velocity_" + std::to_string(step);
            std::string pname = "fluid_pressure_" + std::to_string(step);
            gsWriteParaview(velocity, vname, samp);
            gsWriteParaview(pressure, pname, samp);

            if (solver.isALEActive())
            {
                gsField<> meshDisp = solver.getMeshDisplacementField();
                std::string mname = "mesh_displacement_" + std::to_string(step);
                gsWriteParaview(meshDisp, mname, samp);
            }
        }
    }
    
    gsInfo << "\nSimulation completed!\n";
    
    // No need to delete - nsPde is a stack object
    return 0;
}