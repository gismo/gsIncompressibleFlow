/** @file test_rotation_debug.cpp
    @brief Debug test for rotation with ALE
*/

#include <gismo.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    // Simple test parameters
    real_t Re = 100;
    real_t meanVelocity = 1.0;
    real_t L = 1.0;
    real_t nu = meanVelocity * L / Re;
    real_t timeStep = 0.01;
    real_t timeSpan = 0.1;  // Short test
    
    gsInfo << "=== Debug Rotation Test ===\n";
    gsInfo << "Reynolds number: " << Re << "\n";
    gsInfo << "Kinematic viscosity: " << nu << "\n";
    gsInfo << "Time step: " << timeStep << "\n\n";
    
    // Create simple 2-patch domain
    gsMultiPatch<> fluidDomain;
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // Left patch
    fluidDomain.patch(0).coefs() << 0,0, 0.5,0, 0,1, 0.5,1;
    
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1)); // Right patch  
    fluidDomain.patch(1).coefs() << 0.5,0, 1,0, 0.5,1, 1,1;
    
    fluidDomain.computeTopology();
    fluidDomain.uniformRefine(2);
    
    // Create basis
    gsMultiPatch<> pressureDomain = fluidDomain;
    fluidDomain.degreeElevate(1);
    
    gsMultiBasis<> basisVelocity(fluidDomain);
    gsMultiBasis<> basisPressure(pressureDomain);
    
    // Boundary conditions
    gsBoundaryConditions<> bcInfo;
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    gsConstantFunction<> inletVel(meanVelocity, 0.0, 2);
    
    // Inlet
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &inletVel, 0);
    
    // Walls
    for (index_t d = 0; d < 2; ++d)
    {
        bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, &zeroVel, 0, d);
        bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &zeroVel, 0, d);
    }
    
    // Create PDE
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &zeroVel, nu);
    
    // Setup solver
    std::vector<gsMultiBasis<>> discreteBases;
    discreteBases.push_back(basisVelocity);
    discreteBases.push_back(basisPressure);
    
    gsFlowSolverParams<real_t> params(nsPde, discreteBases);
    params.options().setInt("nonlin.maxIt", 10);
    params.options().setReal("nonlin.tol", 1e-6);
    params.options().setReal("timeStep", timeStep);
    params.options().setSwitch("quiet", false);
    
    // Test 1: Non-ALE solver
    gsInfo << "\n--- Test 1: Non-ALE Solver ---\n";
    {
        gsINSSolverUnsteady<> nonALESolver(memory::make_shared_not_owned(&params));
        nonALESolver.initialize();
        nonALESolver.solveStokes();
        
        gsMatrix<> stokesVel = nonALESolver.solutionCoefs(0);
        gsMatrix<> stokesPres = nonALESolver.solutionCoefs(1);
        gsInfo << "Stokes: Max velocity = " << stokesVel.lpNorm<gsEigen::Infinity>() 
               << ", Max pressure = " << stokesPres.lpNorm<gsEigen::Infinity>() << "\n";
        
        for (int i = 0; i < 3; ++i)
        {
            nonALESolver.nextIteration();
            gsMatrix<> vel = nonALESolver.solutionCoefs(0);
            gsMatrix<> pres = nonALESolver.solutionCoefs(1);
            gsInfo << "Step " << i << ": Max velocity = " << vel.lpNorm<gsEigen::Infinity>()
                   << ", Max pressure = " << pres.lpNorm<gsEigen::Infinity>() << "\n";
        }
    }
    
    // Test 2: ALE solver with zero mesh motion
    gsInfo << "\n--- Test 2: ALE Solver (Zero Motion) ---\n";
    {
        gsINSSolverUnsteadyALE<> aleSolver(memory::make_shared_not_owned(&params));
        aleSolver.initialize();
        aleSolver.setALEActive(true);
        aleSolver.solveStokes();
        
        gsMatrix<> stokesVel = aleSolver.solutionCoefs(0);
        gsMatrix<> stokesPres = aleSolver.solutionCoefs(1);
        gsInfo << "Stokes: Max velocity = " << stokesVel.lpNorm<gsEigen::Infinity>() 
               << ", Max pressure = " << stokesPres.lpNorm<gsEigen::Infinity>() << "\n";
        
        // Zero mesh motion
        auto zeroMotion = [&](real_t t) -> gsMatrix<real_t>
        {
            const index_t udofs = aleSolver.getAssembler()->getUdofs();
            gsMatrix<> disp(2 * udofs, 1); 
            disp.setZero();
            return disp;
        };
        aleSolver.setMeshUpdateFunction(zeroMotion);
        
        for (int i = 0; i < 3; ++i)
        {
            aleSolver.nextIteration();
            gsMatrix<> vel = aleSolver.solutionCoefs(0);
            gsMatrix<> pres = aleSolver.solutionCoefs(1);
            gsInfo << "Step " << i << ": Max velocity = " << vel.lpNorm<gsEigen::Infinity>()
                   << ", Max pressure = " << pres.lpNorm<gsEigen::Infinity>() << "\n";
        }
    }
    
    // Test 3: ALE solver with small mesh motion
    gsInfo << "\n--- Test 3: ALE Solver (Small Motion) ---\n";
    {
        gsINSSolverUnsteadyALE<> aleSolver(memory::make_shared_not_owned(&params));
        aleSolver.initialize();
        aleSolver.setALEActive(true);
        aleSolver.solveStokes();
        
        gsMatrix<> stokesVel = aleSolver.solutionCoefs(0);
        gsMatrix<> stokesPres = aleSolver.solutionCoefs(1);
        gsInfo << "Stokes: Max velocity = " << stokesVel.lpNorm<gsEigen::Infinity>() 
               << ", Max pressure = " << stokesPres.lpNorm<gsEigen::Infinity>() << "\n";
        
        // Small oscillating motion
        auto smallMotion = [&](real_t t) -> gsMatrix<real_t>
        {
            const index_t udofs = aleSolver.getAssembler()->getUdofs();
            const gsDofMapper& mapper = aleSolver.getAssembler()->getMappers()[0];
            gsMatrix<> disp(2 * udofs, 1); 
            disp.setZero();
            
            // Small vertical motion for middle nodes
            real_t amplitude = 0.01 * std::sin(2 * M_PI * t / 0.1);
            
            for (index_t p = 0; p < 2; ++p)
            {
                const gsMatrix<>& coefs = fluidDomain.patch(p).coefs();
                for (index_t i = 0; i < coefs.rows(); ++i)
                {
                    if (mapper.is_free(i, p))
                    {
                        real_t x = coefs(i, 0);
                        // Only move interior nodes
                        if (x > 0.1 && x < 0.9)
                        {
                            index_t idx = mapper.index(i, p);
                            disp(idx + udofs) = amplitude; // y-displacement
                        }
                    }
                }
            }
            
            return disp;
        };
        aleSolver.setMeshUpdateFunction(smallMotion);
        
        for (int i = 0; i < 3; ++i)
        {
            aleSolver.nextIteration();
            gsMatrix<> vel = aleSolver.solutionCoefs(0);
            gsMatrix<> pres = aleSolver.solutionCoefs(1);
            gsField<> meshVelField = aleSolver.getMeshVelocityField();
            real_t maxMeshVel = meshVelField.coefficientVector().template lpNorm<gsEigen::Infinity>();
            
            gsInfo << "Step " << i << ": Max velocity = " << vel.lpNorm<gsEigen::Infinity>()
                   << ", Max pressure = " << pres.lpNorm<gsEigen::Infinity>()
                   << ", Max mesh vel = " << maxMeshVel << "\n";
        }
    }
    
    return 0;
}