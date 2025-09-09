/** @file simple_rotation_test.cpp
    @brief Simple test for rotation with ALE - channel flow with rotating obstacle
*/

#include <gismo.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    // Parameters
    real_t Re = 100;
    real_t meanVelocity = 1.0;
    real_t nu = 0.01;
    real_t timeStep = 0.01;
    index_t nSteps = 20;
    real_t rotationSpeed = 10.0; // degrees per time unit
    
    gsInfo << "=== Simple Rotation Test with ALE ===\n";
    gsInfo << "Reynolds number: " << Re << "\n";
    gsInfo << "Time step: " << timeStep << "\n";
    gsInfo << "Rotation speed: " << rotationSpeed << " deg/s\n\n";
    
    // Create a simple channel with an obstacle in the middle
    // Domain: [0,3] x [0,1] with obstacle at [1.3,1.7] x [0.3,0.7]
    gsMultiPatch<> fluidDomain;
    
    // Create 4 patches around the obstacle
    // Left channel
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1));
    fluidDomain.patch(0).coefs() << 0,0, 1.3,0, 0,1, 1.3,1;
    
    // Bottom channel  
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1));
    fluidDomain.patch(1).coefs() << 1.3,0, 1.7,0, 1.3,0.3, 1.7,0.3;
    
    // Top channel
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1));
    fluidDomain.patch(2).coefs() << 1.3,0.7, 1.7,0.7, 1.3,1, 1.7,1;
    
    // Right channel
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1));
    fluidDomain.patch(3).coefs() << 1.7,0, 3,0, 1.7,1, 3,1;
    
    fluidDomain.computeTopology();
    
    // Refine
    for (int i = 0; i < 3; ++i)
        fluidDomain.uniformRefine();
    
    // Store original domain
    gsMultiPatch<> originalDomain = fluidDomain;
    
    // Create basis  
    gsMultiPatch<> pressureDomain = fluidDomain;
    fluidDomain.degreeElevate(1);
    
    gsMultiBasis<> basisVelocity(fluidDomain);
    gsMultiBasis<> basisPressure(pressureDomain);
    
    // Boundary conditions
    gsBoundaryConditions<> bcInfo;
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    gsConstantFunction<> inletVel(meanVelocity, 0.0, 2);
    
    // Inlet (left)
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &inletVel, 0);
    
    // Walls (top and bottom of channel)
    for (index_t d = 0; d < 2; ++d)
    {
        bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(3, boundary::south, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &zeroVel, 0, d);
    }
    
    // Obstacle boundaries (no explicit BC - will move with mesh)
    
    // Create PDE
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &zeroVel, nu);
    
    // Setup solver
    std::vector<gsMultiBasis<>> discreteBases;
    discreteBases.push_back(basisVelocity);
    discreteBases.push_back(basisPressure);
    
    gsFlowSolverParams<real_t> params(nsPde, discreteBases);
    params.options().setInt("nonlin.maxIt", 20);
    params.options().setReal("nonlin.tol", 1e-6);
    params.options().setReal("timeStep", timeStep);
    params.options().setSwitch("quiet", false);
    
    // Create ALE solver
    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));
    solver.initialize();
    solver.setALEActive(true);
    
    // Define mesh motion - rotate obstacle boundary nodes
    auto meshMotion = [&](real_t t) -> gsMatrix<real_t>
    {
        const index_t udofs = solver.getAssembler()->getUdofs();
        const gsDofMapper& mapper = solver.getAssembler()->getMappers()[0];
        gsMatrix<> disp(2 * udofs, 1); 
        disp.setZero();
        
        // Rotation parameters
        real_t angle = rotationSpeed * t * M_PI / 180.0; // Convert to radians
        real_t cos_a = std::cos(angle);
        real_t sin_a = std::sin(angle);
        real_t cx = 1.5; // Center of obstacle
        real_t cy = 0.5;
        
        // Move nodes on obstacle boundaries
        for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
        {
            const gsMatrix<>& originalCoefs = originalDomain.patch(p).coefs();
            
            for (index_t i = 0; i < originalCoefs.rows(); ++i)
            {
                if (!mapper.is_free(i, p))
                    continue;
                
                real_t x = originalCoefs(i, 0);
                real_t y = originalCoefs(i, 1);
                
                // Check if on obstacle boundary (with tolerance)
                real_t tol = 1e-6;
                bool onObstacle = false;
                
                // Check each patch's boundaries that touch the obstacle
                if (p == 0 && std::abs(x - 1.3) < tol && y >= 0.3-tol && y <= 0.7+tol) // Right boundary of left patch
                    onObstacle = true;
                else if (p == 1 && y >= 0.3-tol && y <= 0.3+tol && x >= 1.3-tol && x <= 1.7+tol) // Top boundary of bottom patch  
                    onObstacle = true;
                else if (p == 2 && y >= 0.7-tol && y <= 0.7+tol && x >= 1.3-tol && x <= 1.7+tol) // Bottom boundary of top patch
                    onObstacle = true;
                else if (p == 3 && std::abs(x - 1.7) < tol && y >= 0.3-tol && y <= 0.7+tol) // Left boundary of right patch
                    onObstacle = true;
                
                if (onObstacle)
                {
                    // Apply rotation
                    real_t xr = cos_a * (x - cx) - sin_a * (y - cy) + cx;
                    real_t yr = sin_a * (x - cx) + cos_a * (y - cy) + cy;
                    
                    index_t idx = mapper.index(i, p);
                    disp(idx) = xr - x;
                    disp(idx + udofs) = yr - y;
                }
            }
        }
        
        return disp;
    };
    
    solver.setMeshUpdateFunction(meshMotion);
    
    // Solve Stokes
    gsInfo << "Solving initial Stokes problem...\n";
    solver.solveStokes();
    
    gsMatrix<> stokesVel = solver.solutionCoefs(0);
    gsMatrix<> stokesPres = solver.solutionCoefs(1);
    gsInfo << "Initial: Max velocity = " << stokesVel.lpNorm<gsEigen::Infinity>() 
           << ", Max pressure = " << stokesPres.lpNorm<gsEigen::Infinity>() << "\n\n";
    
    // Time stepping
    for (index_t step = 1; step <= nSteps; ++step)
    {
        real_t time = step * timeStep;
        gsInfo << "Time step " << step << ", t = " << time 
               << ", rotation = " << rotationSpeed * time << " degrees\n";
        
        solver.nextIteration();
        
        gsMatrix<> vel = solver.solutionCoefs(0);
        gsMatrix<> pres = solver.solutionCoefs(1);
        
        gsField<> meshVelField = solver.getMeshVelocityField();
        real_t maxMeshVel = meshVelField.coefficientVector().template lpNorm<gsEigen::Infinity>();
        
        gsInfo << "  Max velocity = " << vel.lpNorm<gsEigen::Infinity>()
               << ", Max pressure = " << pres.lpNorm<gsEigen::Infinity>()
               << ", Max mesh vel = " << maxMeshVel << "\n";
               
        // Write output every 5 steps
        if (step % 5 == 0)
        {
            gsField<> velField = solver.constructSolution(vel);
            gsField<> presField = solver.constructPressure(pres);
            
            gsWriteParaview(velField, "simple_rotation_velocity_" + std::to_string(step), 1000, true);
            gsWriteParaview(presField, "simple_rotation_pressure_" + std::to_string(step), 1000, true);
            gsWriteParaview(solver.getAssembler()->getPatches(), "simple_rotation_mesh_" + std::to_string(step), 1000, true);
        }
    }
    
    return 0;
}