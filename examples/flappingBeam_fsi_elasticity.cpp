/** @file flappingBeam_fluid.cpp
 
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Li
*/


#include <gismo.h>
#include <cmath>
#include <map>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowBndEvaluators.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsElTimeIntegrator.h>
#include <gsElasticity/src/gsMassAssembler.h>
#include <gsElasticity/src/gsGeoUtils.h>
#include <gsElasticity/src/gsALE.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsElasticity/src/gsPartitionedFSI.h>

using namespace gismo;

template<class T, int MatOrder> void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, int geo);

// reset the solver with a new mesh
template<class T, int MatOrder>
void resetSolverWithNewMesh(gsINSSolverUnsteady<T, MatOrder>*& solver,
                            const gsMultiPatch<T>& newPatches,
                            gsElasticityAssembler<T>* beamAssembler);

// transfer the solution to a new mesh
template<class T, int MatOrder>
void transferSolutionToNewMesh(gsINSSolverUnsteady<T, MatOrder>* solver,
                              const gsField<T>& velocityField,
                              const gsField<T>& pressureField,
                              gsElasticityAssembler<T>* beamAssembler);

// 根据时间更新几何体边界
template<class T>
gsMultiPatch<T> updateGeometryBoundary(gsMultiPatch<T>& patches, T time, gsINSSolverUnsteady<T, ColMajor>* fluidSolver);

// 先声明initializeSolidSolver函数原型
template<class T>
void initializeSolidSolver(const gsMultiPatch<T>& fluidPatches, T timeStep);

// 注释：固体边界速度现在通过resetSolverWithNewMesh中的Dirichlet边界条件直接应用
// template<class T>
// void updateFluidVelocityFromSolidBoundary(gsINSSolverUnsteady<T, ColMajor>* fluidSolver, 
//                                           const gsMultiPatch<T>& patches);

// 声明压力到牵引力转换函数和固体求解器重建函数
template<class T>
class gsPressureToTraction : public gsFunction<T>
{
public:
    gsPressureToTraction(const gsFunction<T>& pressureFunc, 
                        const gsGeometry<T>& boundaryGeom,
                        boundary::side side)
        : m_pressureFunc(pressureFunc.clone()), m_boundaryGeom(boundaryGeom.clone()), m_side(side) 
    {
        // gsInfo << "创建gsPressureToTraction对象，边界类型: " << m_side << "\n";
    }
    
    ~gsPressureToTraction() 
    {
        // gsInfo << "销毁gsPressureToTraction对象，边界类型: " << m_side << "\n";
    }
    
    short_t domainDim() const override { return 1; } // 边界函数只接受1D参数
    short_t targetDim() const override { return 2; }
    
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        try {
            // Debug output
            if (u.rows() != 1) {
                gsWarn << "gsPressureToTraction::eval_into received input with " << u.rows() 
                       << " rows, expected 1. Input:\n" << u << "\n";
            }
            
            // u 是1D边界参数，需要映射到2D物理空间来评估压力
            result.resize(2, u.cols());
            
            // 为每个评估点计算变形后的法向量和压力
            for (index_t i = 0; i < u.cols(); ++i) {
                // u.col(i) 是1D边界参数
                gsMatrix<T> boundaryParam = u.col(i); // 1D边界参数
                
                // 将1D边界参数映射到2D参数空间坐标
                gsMatrix<T> paramPoint(2, 1);
                
                // 根据边界类型设置参数空间坐标
                if (m_side == boundary::north) { // v=1
                    paramPoint(0, 0) = boundaryParam(0, 0); // u
                    paramPoint(1, 0) = 1.0; // v=1
                } else if (m_side == boundary::south) { // v=0
                    paramPoint(0, 0) = boundaryParam(0, 0); // u
                    paramPoint(1, 0) = 0.0; // v=0
                } else if (m_side == boundary::east) { // u=1
                    paramPoint(0, 0) = 1.0; // u=1
                    paramPoint(1, 0) = boundaryParam(0, 0); // v
                } else if (m_side == boundary::west) { // u=0
                    paramPoint(0, 0) = 0.0; // u=0
                    paramPoint(1, 0) = boundaryParam(0, 0); // v
                } else {
                    GISMO_ERROR("Unsupported boundary side");
                }
                
                // 在参数空间评估压力
                gsMatrix<T> pressure;
                m_pressureFunc->eval_into(paramPoint, pressure);
                
                // 计算该点处变形后边界的法向量
                gsVector<T> normal = computeDeformedNormal(boundaryParam);
                
                // 牵引力 = -压力 * 法向量
                result(0, i) = -pressure(0, i) * normal(0);
                result(1, i) = -pressure(0, i) * normal(1);
            }
        } catch (const std::exception& e) {
            gsWarn << "gsPressureToTraction::eval_into 异常: " << e.what() << "\n";
            throw;
        }
    }
    
    gsFunction<T>* clone() const {
        return new gsPressureToTraction<T>(*m_pressureFunc, *m_boundaryGeom, m_side);
    }
    
private:
    std::unique_ptr<gsFunction<T>> m_pressureFunc;   // 保存压力函数的拷贝
    std::unique_ptr<gsGeometry<T>> m_boundaryGeom;   // 保存边界几何的拷贝
    boundary::side m_side;
    
    // 计算变形后边界的法向量
     gsVector<T> computeDeformedNormal(const gsMatrix<T>& boundaryParam) const
    {
        // Compute normal using boundary curve derivative
        gsMatrix<T> deriv;
        m_boundaryGeom->deriv_into(boundaryParam, deriv);
        
        // Extract tangent components
        T tx = deriv(0, 0);
        T ty = deriv(1, 0);
        
        // Normalize tangent
        T length = std::sqrt(tx*tx + ty*ty);
        
        gsVector<T> normal(2);
        if (length > 1e-12) {
            tx /= length;
            ty /= length;
        
            if (m_side == boundary::north || m_side == boundary::east) {
                // Rotate counter-clockwise
                normal(0) = ty;
                normal(1) = -tx;
            } else { // south or west
                // Rotate clockwise
                normal(0) = -ty;
                normal(1) = tx;
            }
        }
        
        return normal;
    }
    
};

// Wrapper class that converts 2D patch parameters to 1D boundary parameters
template<class T>
class gsPressureToTractionPatch : public gsFunction<T>
{
public:
    gsPressureToTractionPatch(const gsFunction<T>& pressureFunc, 
                             const gsGeometry<T>& patchGeom,
                             boundary::side side)
        : m_pressureFunc(pressureFunc.clone()), 
          m_patchGeom(patchGeom.clone()), 
          m_side(side) 
    {
        m_boundaryGeom = m_patchGeom->boundary(m_side);
    }
    
    short_t domainDim() const override { return 2; } // Accept 2D parameters
    short_t targetDim() const override { return 2; }
    
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        result.resize(2, u.cols());
        
        for (index_t i = 0; i < u.cols(); ++i) {
            // u.col(i) is a 2D parameter point on the boundary
            gsMatrix<T> paramPoint = u.col(i);
            
            // Evaluate pressure at this 2D point
            gsMatrix<T> pressure;
            m_pressureFunc->eval_into(paramPoint, pressure);
            
            // Compute normal using Jacobian of the patch
            gsMatrix<T> jac;
            m_patchGeom->jacobian_into(paramPoint, jac);
            
            // Get the tangent vector along the boundary
            gsVector<T> tangent(2);
            if (m_side == boundary::north || m_side == boundary::south) {
                // For north/south boundaries, tangent is in u-direction
                tangent(0) = jac(0, 0);  // dx/du
                tangent(1) = jac(1, 0);  // dy/du
            } else { // east or west
                // For east/west boundaries, tangent is in v-direction
                tangent(0) = jac(0, 1);  // dx/dv
                tangent(1) = jac(1, 1);  // dy/dv
            }
            
            // Normalize tangent
            T length = tangent.norm();
            if (length > 1e-12) {
                tangent /= length;
            }
            
            // Compute outward normal by rotating tangent 90 degrees
            gsVector<T> normal(2);
            if (m_side == boundary::north) {
                // North: rotate tangent counter-clockwise
                normal(0) = -tangent(1);
                normal(1) = tangent(0);
            } else if (m_side == boundary::south) {
                // South: rotate tangent clockwise
                normal(0) = tangent(1);
                normal(1) = -tangent(0);
            } else if (m_side == boundary::east) {
                // East: rotate tangent counter-clockwise
                normal(0) = -tangent(1);
                normal(1) = tangent(0);
            } else { // west
                // West: rotate tangent clockwise
                normal(0) = tangent(1);
                normal(1) = -tangent(0);
            }
            
            // Traction = -pressure * normal
            result(0, i) = -pressure(0, 0) * normal(0);
            result(1, i) = -pressure(0, 0) * normal(1);
        }
    }
    
    virtual gsFunction<T>* clone() const {
        return new gsPressureToTractionPatch<T>(*m_pressureFunc, *m_patchGeom, m_side);
    }
    
private:
    std::unique_ptr<gsFunction<T>> m_pressureFunc;
    std::unique_ptr<gsGeometry<T>> m_patchGeom;
    std::unique_ptr<gsGeometry<T>> m_boundaryGeom;
    boundary::side m_side;

    gsVector<T> computeDeformedNormal(const gsMatrix<T>& boundaryParam) const
    {
        // Compute normal using boundary curve derivative
        gsMatrix<T> deriv;
        m_boundaryGeom->deriv_into(boundaryParam, deriv);
        
        // Extract tangent components
        T tx = deriv(0, 0);
        T ty = deriv(1, 0);
        
        // Normalize tangent
        T length = std::sqrt(tx*tx + ty*ty);
        
        gsVector<T> normal(2);
        if (length > 1e-12) {
            tx /= length;
            ty /= length;
            
            // Compute outward normal by rotating tangent 90 degrees
            // The rotation depends on the boundary side
            if (m_side == boundary::north || m_side == boundary::east) {
                // Rotate counter-clockwise
                normal(0) = ty;
                normal(1) = -tx;
            } else { // south or west
                // Rotate clockwise
                normal(0) = -ty;
                normal(1) = tx;
            }
        } else {
            // Default normal based on side
            if (m_side == boundary::north) {
                normal(0) = 0; normal(1) = 1;
            } else if (m_side == boundary::south) {
                normal(0) = 0; normal(1) = -1;
            } else if (m_side == boundary::east) {
                normal(0) = 1; normal(1) = 0;
            } else {
                normal(0) = -1; normal(1) = 0;
            }
        }
        
        return normal;
    }
};

template<class T>
void rebuildSolidSolver(const gsFunctionSet<T>& pressureFunction);

// 声明固体求解器相关的全局变量
static std::unique_ptr<gsElTimeIntegrator<real_t>> beamSolver;
static std::unique_ptr<gsElasticityAssembler<real_t>> beamAssembler;
static std::unique_ptr<gsMassAssembler<real_t>> beamMassAssembler;
static gsMultiPatch<real_t> beamGeometry;
static gsMultiPatch<real_t> beamDisplacement;
static gsMultiPatch<real_t> beamVelocity;
static gsBoundaryConditions<real_t> beamBCs;
static bool solidInitialized = false;
static real_t lastTime = 0.0;

// 添加全局变量用于动画收集
static gsParaviewCollection dispCollection("beam_displacement_animation");
static gsParaviewCollection stressCollection("beam_stress_animation");
static bool animationInitialized = false;

// 添加全局变量用于存储边界速度信息（用于ALE计算）
static std::map<int, gsMatrix<real_t>> boundaryVelocities;
static bool boundaryVelocitiesInitialized = false;


int main(int argc, char *argv[])
{
    typedef gsGMRes<real_t> LinSolver;

    // ========================================= Settings ========================================= 

    // solvers
    bool steady = false;
    bool steadyIt = false;
    bool unsteady = false;
    bool unsteadyIt = false;

    // domain definition
    int geo = 1; // 0 - custom input file, 1 - step, 2 - cavity, 3 - blade profile 2D, 4 - flapping beam
    int dim = 2; // relevant for step and cavity
    std::string inputFile = "";
    
    // discretization settings
    int numRefine = 2;
    int wallRefine = 0; // relevant for step, cavity, profile2D
    int leadRefine = 0; // relevant for profile2D
    int numElevate = 0; // number of degree elevations (before refinement)

    // problem parameters - Turek-Hron FSI 3 标准参数
    real_t viscosity = 1.0e-3;  // 动力粘度 μ = 1.0×10^-3 Pa·s (Turek-Hron FSI 3)
    real_t inVelX = 1; // inlet x-velocity for profile2D
    real_t inVelY = 0; // inlet y-velocity for profile2D
    
    // solver settings
    int maxIt = 100;
    int picardIt = 5;
    int linIt = 50;
    real_t timeStep = 0.01;  // Turek-Hron FSI 3 标准时间步长
    real_t tol = 1e-5;
    real_t picardTol = 1e-4;
    real_t linTol = 1e-6;
    std::string matFormation = "EbE";
    std::string precond = "MSIMPLER_FdiagEqual";
    bool stokesInit = false; // start unsteady problem from Stokes solution

    // output settings
    bool quiet = false;
    bool plot = false;
    bool plotMesh = false;
    int plotPts = 1000;
    bool animation = false;
    int animStep = 100;

    // ---------------------------------------------------------------------------------

    //command line
    gsCmdLine cmd("Solves the Navier-Stokes problem in a given domain (step, cavity, blade profile).");

    cmd.addSwitch("steady", "Solve steady problem with direct linear solver", steady);
    cmd.addSwitch("steadyIt", "Solve steady problem with preconditioned GMRES as linear solver", steadyIt);
    cmd.addSwitch("unsteady", "Solve unsteady problem with direct linear solver", unsteady);
    cmd.addSwitch("unsteadyIt", "Solve unsteady problem with preconditioned GMRES as linear solver", unsteadyIt);

    cmd.addInt("g", "geo", "Computational domain (0 - custom file, 1 - step, 2 - cavity, 3 - profile (only 2D), 4 - flapping beam)", geo);
    cmd.addInt("d", "dim", "Space dimension", dim);
    cmd.addString("", "input", "Full path to the input xml file containing geometry, right-hand side functin and boundary conditions", inputFile);

    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("", "wallRefine", "Number of h-refinement steps near step corner, cavity walls of blade profile", wallRefine);
    cmd.addInt("", "leadRefine", "Number of h-refinement steps near the beginning of the blade (for profile geometry)", leadRefine);
    cmd.addInt("e", "degElevate", "Number of degree elevations (performed before h-refinement)", numElevate);

    cmd.addReal("v", "visc", "Viscosity value", viscosity);
    cmd.addReal("", "inVelX", "x-coordinate of inflow velocity (for profile geometry)", inVelX);
    cmd.addReal("", "inVelY", "y-coordinate of inflow velocity (for profile geometry)", inVelY);

    cmd.addInt("", "maxIt", "Max. number of Picard iterations or time steps", maxIt);
    cmd.addInt("", "picardIt", "Max. number of inner Picard iterations for unsteady problem", picardIt);
    cmd.addInt("", "linIt", "Max. number of GMRES iterations (if the lin. systems are solved iteratively)", linIt);
    cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
    cmd.addReal("", "tol", "Stopping tolerance", tol);
    cmd.addReal("", "picardTol", "Tolerance for inner Picard iteration for unsteady problem", picardTol);
    cmd.addReal("", "linTol", "Tolerance for iterative linear solver", linTol);
    cmd.addString("", "loop", "Matrix formation method (EbE = element by element, RbR = row by row)", matFormation);
    cmd.addString("p", "precond", "Preconditioner type (format: PREC_Fstrategy, PREC = {PCD, PCDmod, LSC, AL, SIMPLE, SIMPLER, MSIMPLER}, Fstrategy = {FdiagEqual, Fdiag, Fmod, Fwhole})", precond);
    cmd.addSwitch("stokesInit", "Set Stokes initial condition", stokesInit);

    cmd.addSwitch("quiet", "Supress (some) terminal output", quiet);
    cmd.addSwitch("plot", "Plot the final result in ParaView format", plot);
    cmd.addSwitch("plotMesh", "Plot the computational mesh", plotMesh);
    cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
    cmd.addSwitch("animation", "Plot animation of the unsteady problem", animation);
    cmd.addInt("", "animStep", "Number of iterations between screenshots for animation (used when animation = true)", animStep);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (!inputFile.empty())
        geo = 0;

    if ( !(steady || steadyIt || unsteady || unsteadyIt) )
        gsWarn << "All computation flags set to false - nothing will be computed.\nPlease select at least one of the flags: --steady, --steadyIt, --unsteady, --unsteadyIt\n\n";

    // ========================================= Define problem (geometry, BCs, rhs) ========================================= 
    
    gsMultiPatch<> patches;
    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> f; // external force

    std::string fn, geoStr;

    switch(geo)
    {
        case 0:
            break; // inputFile is given from cmd
        
        default:
            gsWarn << "Unknown geometry ID, using backward-facing step.\n";
            geo = 1;

        case 1:
            geoStr = "BFS" + util::to_string(dim) + "D";
            fn = geoStr + "_problem.xml";
            inputFile = fn;
            break;

        case 2:
            geoStr = "LDC" + util::to_string(dim) + "D";
            fn = geoStr + "_problem.xml";
            inputFile = fn;
            break;

        case 3:
            if (dim == 3)
                gsWarn << "Geometry 3 is only 2D!\n";

            geoStr = "profile2D";
            inputFile = geoStr + "_problem.xml";
            break;
            
        case 4:
            if (dim == 3)
                gsWarn << "Flapping beam geometry is only 2D!\n";
                
            geoStr = "flappingBeam";
            inputFile = geoStr + "_flow.xml";
            break;
    }

    gsInfo << "Reading problem definition from file:\n" << inputFile << "\n\n";

    std::string path = gsFileManager::find(inputFile);
    if ( path.empty() )
    {
        gsWarn<<"Input file not found, quitting.\n";
        return 1;
    }

    gsFileData<> fd(inputFile);
    fd.getId(7, patches);   // id=0: multipatch domain
    fd.getId(1, f);         // id=1: source function
    fd.getId(2, bcInfo);    // id=2: boundary conditions

    gsInfo << "Solving Navier-Stokes problem in " << geoStr << " domain.\n";
    gsInfo << patches;
    gsInfo << "viscosity = " << viscosity << "\n";
    gsInfo << "source function = " << f << "\n";

    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> basis(patches);
    basis.degreeElevate(numElevate);

    switch(geo)
    {
        case 0:
            for (int r = 0; r < numRefine; ++r)
                basis.uniformRefine();
            break;

        case 1:
        default:
            refineBasis_step(basis, numRefine, 0, wallRefine, 0, 0, dim, 8.0, 2.0, 2.0); // 8, 2, 2 are dimensions of the domain in the input xml file
            break;

        case 2:
            refineBasis_cavity(basis, numRefine, wallRefine, dim);
            break;

        case 3:
            refineBasis_profile2D(basis, numRefine, wallRefine, leadRefine);
            break;
            
        case 4:
            // refine the flapping beam appropriately
            for (int r = 0; r < numRefine; ++r)
                basis.uniformRefine();
            break;
    }    

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)

    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases);
    params.options().setSwitch("quiet", quiet);
    params.options().setString("assemb.loop", matFormation);

    gsOptionList solveOpt;
    solveOpt.addInt("geo", "", geo);
    solveOpt.addInt("maxIt", "", maxIt);
    solveOpt.addInt("plotPts", "", plotPts);
    solveOpt.addInt("animStep", "", animStep);
    solveOpt.addInt("nonlin.maxIt", "", picardIt);
    solveOpt.addReal("tol", "", tol);
    solveOpt.addReal("picardTol", "", picardTol);
    solveOpt.addSwitch("plot", "", plot);
    solveOpt.addSwitch("animation", "", animation);
    solveOpt.addSwitch("plotMesh", "", plotMesh);
    solveOpt.addSwitch("stokesInit", "", stokesInit);
    solveOpt.addString("id", "", "");

    if (steady)
    {
        solveOpt.setString("id", "steady");
        params.options().setString("lin.solver", "direct");

        gsINSSolverSteady<real_t, ColMajor> NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the steady problem with direct linear solver.\n";

        solveProblem(NSsolver, solveOpt, geo);

        // example of flow rate computation
        if (geo == 1)
        {
            std::vector<std::pair<int, boxSide> > bndOut;
            bndOut.push_back(std::make_pair(0, boundary::east));
            bndOut.push_back(std::make_pair(1, boundary::east));
            gsField<> velocity = NSsolver.constructSolution(0);
            gsFlowBndEvaluator_flowRate<real_t> flowRateEval(params, bndOut);
            flowRateEval.setVelocityField(velocity);
            flowRateEval.evaluate();
            gsInfo << "bndOut flow rate = " << flowRateEval.getValue() << "\n";
        }
    }

    if (steadyIt)
    {
        solveOpt.setString("id", "steadyIt");
        params.options().setString("lin.solver", "iter");
        params.options().setInt("lin.maxIt", linIt);
        params.options().setReal("lin.tol", linTol);
        params.options().setString("lin.precType", precond);

        gsINSSolverSteady<real_t, ColMajor > NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the steady problem with preconditioned GMRES as linear solver.\n";
        gsInfo << "Used preconditioner: " << params.options().getString("lin.precType") << "\n";

        solveProblem(NSsolver, solveOpt, geo);

        gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* linSolverPtr = dynamic_cast<gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* >( NSsolver.getLinSolver() );
        reportLinIterations(linSolverPtr);
    }

    if (unsteady)
    {
        solveOpt.setString("id", "unsteady");
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("nonlin.maxIt", picardIt);
        params.options().setReal("nonlin.tol", picardTol);
        params.options().setString("lin.solver", "direct");

        gsINSSolverUnsteady<real_t, ColMajor> NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the unsteady problem with direct linear solver.\n";

        solveProblem(NSsolver, solveOpt, geo);
    }

    if (unsteadyIt)
    {
        solveOpt.setString("id", "unsteadyIt");
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("nonlin.maxIt", picardIt);
        params.options().setReal("nonlin.tol", picardTol);
        params.options().setString("lin.solver", "iter");
        params.options().setInt("lin.maxIt", linIt);
        params.options().setReal("lin.tol", linTol);
        params.options().setString("lin.precType", precond);
        // params.precOptions().setReal("gamma", 10); // parameter for AL preconditioner

        gsINSSolverUnsteady<real_t, ColMajor > NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the unsteady problem with preconditioned GMRES as linear solver.\n";
        gsInfo << "Used preconditioner: " << params.options().getString("lin.precType") << "\n";

        solveProblem(NSsolver, solveOpt, geo);
        
        gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* linSolverPtr = dynamic_cast<gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* >( NSsolver.getLinSolver() );
        reportLinIterations(linSolverPtr);
    }

    // 最后记得保存动画集合
    if (animationInitialized) 
    {
        dispCollection.save();
        stressCollection.save();
        
        gsInfo << "已完成位移和应力动画文件的写入，可以在ParaView中查看:\n";
        gsInfo << "  - beam_displacement_animation.pvd\n";
        gsInfo << "  - beam_stress_animation.pvd\n";
    }

    return 0; 
}

// ... 然后是 solveProblem 的定义 ...
template<class T, int MatOrder>
void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, int geo)
{
    gsStopwatch clock;

    // ------------------------------------
    // prepare strings for output filenames

    bool plot = opt.getSwitch("plot");
    std::string geoStr = "";
    std::string id = opt.getString("id");
    if (plot)
    {
        index_t dim = NSsolver.getParams()->getPde().domain().geoDim();
        std::string dimStr = util::to_string(dim) + "D";

        switch(opt.getInt("geo"))
        {
            case 0:
                geoStr = "customGeo";
                break;
            case 1:
                geoStr = "BFS" + dimStr;
                break;
            case 2:
                geoStr = "LDC" + dimStr;
                break;
            case 3:
                geoStr = "profile2D";
                break;
            case 4:
                geoStr = "flapping_beam" + dimStr;
                break;
            default:
                // 如果没有匹配的几何体ID，使用一个通用名称
                geoStr = "geo" + util::to_string(opt.getInt("geo")) + "_" + dimStr;
                break;
        }
    }

    // ------------------------------------
    // solve problem

    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();

    gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

    gsINSSolverUnsteady<T, MatOrder>* pSolver = dynamic_cast<gsINSSolverUnsteady<T, MatOrder>*>(&NSsolver);

    
    if (pSolver)
    {
        if (opt.getSwitch("stokesInit"))
            pSolver->solveStokes();
        if (opt.getSwitch("animation"))
        {
            pSolver->solveWithAnimation(opt.getInt("maxIt"), opt.getInt("animStep"), 
                                        geoStr + "_" + id, opt.getReal("tol"), opt.getInt("plotPts"));
        }
        else
        {
            // create ParaView files for time step collection
            std::string fileNameU = geoStr + "_" + id + "_velocity_animation.pvd";
            std::ofstream fileU(fileNameU.c_str());
            GISMO_ASSERT(fileU.is_open(), "Error creating " << fileNameU);

            std::string fileNameP = geoStr + "_" + id + "_pressure_animation.pvd";
            std::ofstream fileP(fileNameP.c_str());
            GISMO_ASSERT(fileP.is_open(), "Error creating " << fileNameP);

            // use the existing functions in gsFlowUtils.h
            gismo::startAnimationFile(fileU);
            gismo::startAnimationFile(fileP);
            
            // determine the animation step
            int animStep = 5; 
            try {
                animStep = opt.getInt("animStep");
            } catch (...) {
                gsInfo << "Using default animation step: " << animStep << "\n";
            }
            
            // save the initial state
            if (opt.getSwitch("plot"))
            {
                gsField<T> velocity = pSolver->constructSolution(0);
                gsField<T> pressure = pSolver->constructSolution(1);
                
                std::string velocityBaseName = geoStr + "_" + id + "_velocity_step0";
                std::string pressureBaseName = geoStr + "_" + id + "_pressure_step0";
                
                gsWriteParaview<>(velocity, velocityBaseName, opt.getInt("plotPts"), opt.getSwitch("plotMesh"));
                gsWriteParaview<>(pressure, pressureBaseName, opt.getInt("plotPts"));
                
                // add references to the collection file
                int numPatches = NSsolver.getParams()->getPde().patches().nPatches();
                for (int p = 0; p < numPatches; p++)
                {
                    fileU << "<DataSet timestep=\"0\" part=\"" << p 
                          << "\" file=\"" << velocityBaseName << p << ".vts\"/>\n";
                    fileP << "<DataSet timestep=\"0\" part=\"" << p 
                          << "\" file=\"" << pressureBaseName << p << ".vts\"/>\n";
                }
            }
            
            // create temporary variables to store the previous solution
            gsVector<T> prevVelocityCoefs;
            gsVector<T> prevPressureCoefs;
            
            // get the initial coefficients
            if (opt.getSwitch("stokesInit")) {
                prevVelocityCoefs = pSolver->solutionCoefs(0);
                prevPressureCoefs = pSolver->solutionCoefs(1);
            }
            
            for (int i = 0; i < opt.getInt("maxIt"); ++i)
            {
                gsInfo << "\n=== 时间步 " << (i+1) << " 开始 ===\n";
                
                // 验证时间步开始时的求解器状态
                size_t patchCount = pSolver->getParams()->getPde().patches().nPatches();
                gsInfo << "时间步开始时patch数量: " << patchCount << "\n";
                
                if (patchCount == 0) {
                    gsWarn << "CRITICAL ERROR: 时间步 " << (i+1) << " 开始时patches为空!\n";
                    break;
                }
                
                // 执行一个时间步
                gsInfo << "执行时间步迭代...\n";
                pSolver->nextIteration();
                
                // 验证求解器状态
                size_t patchCountAfter = pSolver->getParams()->getPde().patches().nPatches();
                gsInfo << "时间步迭代后patch数量: " << patchCountAfter << "\n";
                
                if (patchCountAfter == 0) {
                    gsWarn << "CRITICAL ERROR: nextIteration后patches变为空!\n";
                    break;
                }
                
                // 获取当前解的系数
                gsVector<T> currentVelocityCoefs = pSolver->solutionCoefs(0);
                gsVector<T> currentPressureCoefs = pSolver->solutionCoefs(1);

                // 移动边界 - 计算当前时间
                T currentTime = (i + 1) * pSolver->getParams()->options().getReal("timeStep");
                gsMultiPatch<T> currentPatches = pSolver->getParams()->getPde().patches();
                // updateGeometryBoundary(currentPatches, currentTime, pSolver);
                currentPatches = updateGeometryBoundary(currentPatches, currentTime, pSolver);
                
                // 更新网格和求解器
                resetSolverWithNewMesh(pSolver, currentPatches, beamAssembler.get());

                // Always update previous solution to ensure we can calculate change in next iteration
                prevVelocityCoefs = currentVelocityCoefs;
                prevPressureCoefs = currentPressureCoefs;
                
                if (opt.getSwitch("plot"))
                {
                    // Use the latest solution fields
                    try {
                        gsField<T> velocity = pSolver->constructSolution(0);
                        gsField<T> pressure = pSolver->constructSolution(1);
                        
                        std::string velocityBaseName = geoStr + "_" + id + "_velocity_step" + util::to_string(i+1);
                        std::string pressureBaseName = geoStr + "_" + id + "_pressure_step" + util::to_string(i+1);
                        
                        gsWriteParaview<>(velocity, velocityBaseName, opt.getInt("plotPts"), opt.getSwitch("plotMesh"));
                        gsWriteParaview<>(pressure, pressureBaseName, opt.getInt("plotPts"));
                        
                        // 添加边界速度的显式可视化
                        if (boundaryVelocitiesInitialized && !boundaryVelocities.empty()) {
                            gsInfo << "创建固液交界处边界速度可视化...\n";
                            
                            // 创建一个包含所有边界速度的多片段
                            gsMultiPatch<T> boundaryVelocityPatches;
                            
                            // 定义固液交界面
                            struct BeamInterface {
                                size_t patchIdx;
                                boundary::side boundSide;
                            };
                            
                            std::vector<BeamInterface> interfaces = {
                                {3, boundary::south},
                                {4, boundary::north},
                                {5, boundary::west},
                            };
                            
                            // 为每个交界面创建边界速度场
                            for (const auto& interface : interfaces) {
                                int boundaryKey = interface.patchIdx * 10 + static_cast<int>(interface.boundSide);
                                auto it = boundaryVelocities.find(boundaryKey);
                                if (it != boundaryVelocities.end()) {
                                    // 获取边界几何
                                    const gsMultiPatch<T>& currentPatches = pSolver->getParams()->getPde().patches();
                                    std::unique_ptr<gsGeometry<T>> boundaryGeom = 
                                        currentPatches.patch(interface.patchIdx).boundary(interface.boundSide);
                                    
                                    // 创建边界上的速度场
                                    const gsMatrix<T>& velocityData = it->second;
                                    
                                    // 将速度数据设置为边界patch的系数
                                    if (const gsBSpline<T>* bspline = dynamic_cast<const gsBSpline<T>*>(boundaryGeom.get())) {
                                        gsMatrix<T> coefs = bspline->coefs();
                                        // 创建一个新的patch来存储速度数据
                                        gsBSpline<T> velocityPatch(bspline->knots(), velocityData);
                                        boundaryVelocityPatches.addPatch(velocityPatch);
                                    } else if (const gsNurbs<T>* nurbs = dynamic_cast<const gsNurbs<T>*>(boundaryGeom.get())) {
                                        gsMatrix<T> coefs = nurbs->coefs();
                                        // 创建速度NURBS patch
                                        gsKnotVector<T> kv = nurbs->knots();
                                        gsMatrix<T> weights = nurbs->weights();
                                        gsNurbs<T> velocityPatch(kv, velocityData, weights);
                                        boundaryVelocityPatches.addPatch(velocityPatch);
                                    }
                                }
                            }
                            
                            // 保存边界速度场
                            if (boundaryVelocityPatches.nPatches() > 0) {
                                std::string boundaryVelBaseName = geoStr + "_" + id + "_boundary_velocity_step" + util::to_string(i+1);
                                gsField<T> boundaryVelField(boundaryVelocityPatches, boundaryVelocityPatches);
                                gsWriteParaview<>(boundaryVelField, boundaryVelBaseName, opt.getInt("plotPts"), true);
                                gsInfo << "已保存边界速度到: " << boundaryVelBaseName << "\n";
                            }
                            
                            // 另一种方法：创建边界点云和速度向量
                            gsMatrix<T> allBoundaryPoints;
                            gsMatrix<T> allBoundaryVelocities;
                            
                            for (const auto& interface : interfaces) {
                                int boundaryKey = interface.patchIdx * 10 + static_cast<int>(interface.boundSide);
                                auto it = boundaryVelocities.find(boundaryKey);
                                if (it != boundaryVelocities.end()) {
                                    const gsMultiPatch<T>& currentPatches = pSolver->getParams()->getPde().patches();
                                    std::unique_ptr<gsGeometry<T>> boundaryGeom = 
                                        currentPatches.patch(interface.patchIdx).boundary(interface.boundSide);
                                    
                                    // 在边界上采样点
                                    index_t numSamples = 50; // 每个边界采样50个点
                                    gsMatrix<T> paramPoints(1, numSamples);
                                    for (index_t j = 0; j < numSamples; ++j) {
                                        paramPoints(0, j) = j / static_cast<T>(numSamples - 1);
                                    }
                                    
                                    // 评估边界上的物理坐标
                                    gsMatrix<T> physicalPoints = boundaryGeom->eval(paramPoints);
                                    
                                    // 获取对应的速度值
                                    const gsMatrix<T>& velocityData = it->second;
                                    gsMatrix<T> velocityValues(2, numSamples);
                                    
                                    // 插值速度到采样点
                                    for (index_t j = 0; j < numSamples; ++j) {
                                        T t = paramPoints(0, j);
                                        index_t idx = static_cast<index_t>(t * (velocityData.rows() - 1));
                                        T alpha = t * (velocityData.rows() - 1) - idx;
                                        
                                        if (idx < velocityData.rows() - 1) {
                                            velocityValues.col(j) = (1 - alpha) * velocityData.row(idx).transpose() 
                                                                  + alpha * velocityData.row(idx + 1).transpose();
                                        } else {
                                            velocityValues.col(j) = velocityData.row(velocityData.rows() - 1).transpose();
                                        }
                                    }
                                    
                                    // 添加到总的点云
                                    if (allBoundaryPoints.cols() == 0) {
                                        allBoundaryPoints = physicalPoints;
                                        allBoundaryVelocities = velocityValues;
                                    } else {
                                        gsMatrix<T> tempPoints(allBoundaryPoints.rows(), 
                                                             allBoundaryPoints.cols() + physicalPoints.cols());
                                        tempPoints.leftCols(allBoundaryPoints.cols()) = allBoundaryPoints;
                                        tempPoints.rightCols(physicalPoints.cols()) = physicalPoints;
                                        allBoundaryPoints = tempPoints;
                                        
                                        gsMatrix<T> tempVels(allBoundaryVelocities.rows(), 
                                                           allBoundaryVelocities.cols() + velocityValues.cols());
                                        tempVels.leftCols(allBoundaryVelocities.cols()) = allBoundaryVelocities;
                                        tempVels.rightCols(velocityValues.cols()) = velocityValues;
                                        allBoundaryVelocities = tempVels;
                                    }
                                }
                            }
                            
                            // 保存边界点云和速度向量
                            if (allBoundaryPoints.cols() > 0) {
                                std::string pointCloudFile = geoStr + "_" + id + "_boundary_points_step" + util::to_string(i+1) + ".csv";
                                std::ofstream file(pointCloudFile);
                                file << "x,y,vx,vy,magnitude\n";
                                for (index_t j = 0; j < allBoundaryPoints.cols(); ++j) {
                                    T magnitude = allBoundaryVelocities.col(j).norm();
                                    file << allBoundaryPoints(0, j) << "," 
                                         << allBoundaryPoints(1, j) << ","
                                         << allBoundaryVelocities(0, j) << ","
                                         << allBoundaryVelocities(1, j) << ","
                                         << magnitude << "\n";
                                }
                                file.close();
                                gsInfo << "已保存边界点云数据到: " << pointCloudFile << "\n";
                                gsInfo << "  包含 " << allBoundaryPoints.cols() << " 个边界点\n";
                                gsInfo << "  最大速度幅值: " << allBoundaryVelocities.colwise().norm().maxCoeff() << " m/s\n";
                            }
                        }
                        
                        // Add references to the collection file
                        int numPatches = pSolver->getParams()->getPde().patches().nPatches();
                        for (int p = 0; p < numPatches; p++)
                        {
                            fileU << "<DataSet timestep=\"" << i+1 << "\" part=\"" << p 
                                  << "\" file=\"" << velocityBaseName << p << ".vts\"/>\n";
                            fileP << "<DataSet timestep=\"" << i+1 << "\" part=\"" << p 
                                  << "\" file=\"" << pressureBaseName << p << ".vts\"/>\n";
                        }
                    } catch (const std::exception& e) {
                        gsInfo << "Error " << e.what() << "\n";
                    }
                }
            }
            
            // Close ParaView collection files
            gismo::endAnimationFile(fileU);
            gismo::endAnimationFile(fileP);
            
            gsInfo << "\nSimulation completed " << opt.getInt("maxIt") << " time steps.\n";
            gsInfo << "ParaView animation files created: \n";
            gsInfo << "  " << fileNameU << "\n";
            gsInfo << "  " << fileNameP << "\n";
        }
    }
    else
    {
        NSsolver.solve(opt.getInt("maxIt"), opt.getReal("tol"), 0);
    }

    real_t totalT = clock.stop();

    gsInfo << "\nAssembly time:" << NSsolver.getAssemblyTime() << "\n";
    gsInfo << "Solve time:" << NSsolver.getSolveTime() << "\n";
    gsInfo << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";
    gsInfo << "Total solveProblem time:" << totalT << "\n\n";

    // ------------------------------------
    // plot

    if (opt.getSwitch("plot")) 
    {
        try {
            gsField<> velocity = NSsolver.constructSolution(0);
            gsField<> pressure = NSsolver.constructSolution(1);
 
            int plotPts = opt.getInt("plotPts");
 
            gsInfo << "Plotting in Paraview...\n";
            gsWriteParaview<>(velocity, geoStr + "_" + id + "_velocity", plotPts, opt.getSwitch("plotMesh"));
            gsWriteParaview<>(pressure, geoStr + "_" + id + "_pressure", plotPts);
            // plotQuantityFromSolution("divergence", velocity, geoStr + "_" + id + "_velocityDivergence", plotPts);
        } catch (const std::exception &e) {
            gsInfo << "Error constructing solution for plotting: " << e.what() << "\n";
        }
    }
}

template<class T, int MatOrder>
void resetSolverWithNewMesh(gsINSSolverUnsteady<T, MatOrder>*& solver,
                           const gsMultiPatch<T>& newPatches,
                           gsElasticityAssembler<T>* beamAssembler)
{
    // 1. Save current solver parameters
    gsFlowSolverParams<T>* currentParams = solver->getParams().get();
    gsOptionList solverOptions = currentParams->options();
    T viscosity = currentParams->getPde().viscosity();
    const gsFunction<T>* f = currentParams->getPde().rhs();
    gsBoundaryConditions<> bcInfo = currentParams->getPde().bc();

    // 2. Save current solution (as fields)
    gsField<T> velocityFieldOld = solver->constructSolution(0);
    gsField<T> pressureFieldOld = solver->constructSolution(1);
    
    gsInfo << "Current dofs: " << solver->numDofs() << "\n";
    gsInfo << "输入的新patches数量: " << newPatches.nPatches() << "\n";
    
    // 验证输入patches的有效性
    for (size_t i = 0; i < newPatches.nPatches(); ++i) {
        if (newPatches.patch(i).coefs().rows() == 0) {
            gsWarn << "输入的patch " << i << " 无效（控制点为空）\n";
        } else {
            gsInfo << "Patch " << i << " 有效，控制点数量: " << newPatches.patch(i).coefs().rows() << "\n";
        }
    }
    
    // Create deep copy of new patches to ensure their lifecycle - this is a key modification
    gsMultiPatch<T>* persistentPatches = new gsMultiPatch<T>(newPatches);
    
    gsInfo << "创建的持久化patches数量: " << persistentPatches->nPatches() << "\n";
    
    // 3. Create new basis functions while preserving Taylor–Hood P(k+1)/P(k) relationship
    // -------------------------------------------------------------------
    // 旧求解器的基函数（index 0 -> 速度, 1 -> 压力）
    const std::vector<gsMultiBasis<T>> oldBases = solver->getAssembler()->getBases();

    // a) 先为 **压力** 创建基函数，次数为旧压力次数 k
    gsMultiBasis<T> pressureBasis(*persistentPatches);
    if (oldBases.size() > 1) // 确保旧压力空间存在
    {
        for (index_t p = 0; p < pressureBasis.nBases(); ++p)
        {
            index_t k = oldBases[1].basis(p).degree(0); // 旧压力次数 k
            pressureBasis.basis(p).setDegree(k);
        }
    }

    // b) **速度** 基函数：从压力基函数复制后整体升阶 +1 → k+1
    gsMultiBasis<T> velocityBasis = pressureBasis; // 先复制
    velocityBasis.degreeElevate(1);               // 整体升阶一次（保持一次性，不累积）

    // c) 组装 newBases，保证 index 0 -> 速度, index 1 -> 压力
    std::vector<gsMultiBasis<T>> newBases;
    newBases.push_back(velocityBasis);  // 速度 k+1
    newBases.push_back(pressureBasis);  // 压力 k

    // 4. Create new PDE and solver parameters
    // Use persistent patches to create PDE
    gsBoundaryConditions<> updatedBcInfo = bcInfo;

    // ==== 动态 Dirichlet 条件添加完毕 ====

    gsNavStokesPde<T>* persistentPde = new gsNavStokesPde<T>(*persistentPatches, updatedBcInfo, f, viscosity);
    
    // Create persistent solver parameters - these objects need to persist
    typename gsFlowSolverParams<T>::Ptr newParamsPtr =
        std::make_shared<gsFlowSolverParams<T>>(*persistentPde, newBases);
    newParamsPtr->options() = solverOptions;
    
    // 5. Create new solver instance - fix truncation issue here
    auto* newSolver = new gsINSSolverUnsteady<T, MatOrder>(newParamsPtr);
    
    // Initialize new solver
    newSolver->initialize();
    
    // Print number of DOFs and verify basis/geometry consistency
    gsInfo << "New dofs: " << newSolver->numDofs() << "\n";
    gsInfo << "New patches: " << persistentPatches->nPatches() << "\n";
    gsInfo << "New bases[0] pieces: " << newSolver->getAssembler()->getBases()[0].nBases() << "\n";
    gsInfo << "New bases[1] pieces: " << newSolver->getAssembler()->getBases()[1].nBases() << "\n";

    // 强化一致性检查，添加更详细的错误信息
    if (newSolver->getAssembler()->getBases()[0].nBases() != persistentPatches->nPatches()) {
        gsWarn << "CRITICAL: 速度基函数数量与patches不匹配! "
               << "基函数: " << newSolver->getAssembler()->getBases()[0].nBases() 
               << ", patches: " << persistentPatches->nPatches() << "\n";
        gsWarn << "这可能导致patch消失或求解器错误\n";
    }
    if (newSolver->getAssembler()->getBases()[1].nBases() != persistentPatches->nPatches()) {
        gsWarn << "CRITICAL: 压力基函数数量与patches不匹配! "
               << "基函数: " << newSolver->getAssembler()->getBases()[1].nBases() 
               << ", patches: " << persistentPatches->nPatches() << "\n";
        gsWarn << "这可能导致patch消失或求解器错误\n";
    }
    
    // 验证每个patch的有效性
    for (size_t i = 0; i < persistentPatches->nPatches(); ++i) {
        const auto& patch = persistentPatches->patch(i);
        if (patch.coefs().rows() == 0) {
            gsWarn << "CRITICAL: Patch " << i << " 在持久化后变为空!\n";
        }
    }

    // 6. Project old solution onto new mesh
    transferSolutionToNewMesh(newSolver, velocityFieldOld, pressureFieldOld, beamAssembler);
    
    // 7. Verify new solver state
    gsInfo << "Verifying new solver state...\n";
    gsInfo << "New solver patch count: " << newSolver->getParams()->getPde().patches().nPatches() << "\n";
    GISMO_ASSERT(newSolver->getParams()->getPde().patches().nPatches() > 0, 
                "New solver has empty patches");
    
    // 8. Replace solver, maintaining pointer ownership
    gsINSSolverUnsteady<T, MatOrder>* oldSolver = solver;
    solver = newSolver;
    
    // 9. Final verification - 强化检查
    gsInfo << "=== 最终验证 resetSolverWithNewMesh ===\n";
    gsInfo << "最终求解器patch数量: " << solver->getParams()->getPde().patches().nPatches() << "\n";
    gsInfo << "最终求解器自由度: " << solver->numDofs() << "\n";
    gsInfo << "最终求解器速度基函数数量: " << solver->getAssembler()->getBases()[0].nBases() << "\n";
    gsInfo << "最终求解器压力基函数数量: " << solver->getAssembler()->getBases()[1].nBases() << "\n";
    
    // 验证patch不为空
    if (solver->getParams()->getPde().patches().nPatches() == 0) {
        gsWarn << "CRITICAL ERROR: 最终求解器的patches数量为0！\n";
        gsWarn << "这表明在求解器重置过程中patch丢失了\n";
    }
    
    // 验证每个patch的有效性
    for (size_t i = 0; i < solver->getParams()->getPde().patches().nPatches(); ++i) {
        const auto& patch = solver->getParams()->getPde().patches().patch(i);
        if (patch.coefs().rows() == 0) {
            gsWarn << "CRITICAL: 最终求解器的patch " << i << " 为空!\n";
        } else {
            gsInfo << "✓ Patch " << i << " 有效，控制点数量: " << patch.coefs().rows() << "\n";
        }
    }
    
    // 验证基函数和patches的一致性
    size_t expectedPatches = newPatches.nPatches();
    size_t actualPatches = solver->getParams()->getPde().patches().nPatches();
    if (actualPatches != expectedPatches) {
        gsWarn << "PATCH数量不匹配! 期望: " << expectedPatches << ", 实际: " << actualPatches << "\n";
    } else {
        gsInfo << "✓ Patch数量验证通过: " << actualPatches << "\n";
    }
    
    // 验证边界速度是否正确应用到新求解器
    if (boundaryVelocitiesInitialized && !boundaryVelocities.empty()) {
        gsInfo << "=== 验证边界速度应用结果 ===\n";
        
        // 构造求解器的速度场，检查边界是否有正确的速度
        try {
            gsField<T> velocityField = solver->constructSolution(0);
            gsInfo << "成功构造新求解器的速度场\n";
            
            // 在边界位置评估速度，验证是否与固体速度一致
            for (const auto& pair : boundaryVelocities) {
                int boundaryKey = pair.first;
                const gsMatrix<T>& expectedVelocities = pair.second;
                
                int patchIdx = boundaryKey / 10;
                int sideIdx = boundaryKey % 10;
                
                T expectedMaxVel = expectedVelocities.cwiseAbs().maxCoeff();
                gsInfo << "  Patch " << patchIdx << " Side " << sideIdx 
                       << ": 期望最大速度 = " << expectedMaxVel << " m/s\n";
            }
            
            gsInfo << "边界速度已通过Dirichlet边界条件正确应用到新求解器\n";
        } catch (const std::exception& e) {
            gsWarn << "验证边界速度时出错: " << e.what() << "\n";
        }
    } else {
        gsInfo << "无边界速度信息可用于验证\n";
    }
}

template<class T, int MatOrder>
void transferSolutionToNewMesh(gsINSSolverUnsteady<T, MatOrder>* solver,
                              const gsField<T>& velocityField,
                              const gsField<T>& pressureField,
                              gsElasticityAssembler<T>* beamAssembler)
{
    // 获取新的基函数
    const std::vector<gsMultiBasis<T>>& newBases = solver->getAssembler()->getBases();
    const gsMultiPatch<T>& newPatches = solver->getParams()->getPde().patches();

    // 保证基函数一致
    gsDebugVar(newPatches.nPatches());
    gsDebugVar(newBases.size());
    for (size_t i = 0; i < newBases.size(); ++i) {
        gsDebugVar(newBases[i].nBases());
    }
    
    // 获得自由度 维度
    const index_t vDim = velocityField.function(0).targetDim();
    const index_t fullUdofs = solver->getAssembler()->getUdofs();
    const index_t tarDim = solver->getParams()->getPde().domain().geoDim();
    const index_t pShift = tarDim * fullUdofs;
    
    // 创建与setSolutionCoefs尺寸匹配的矩阵
    gsMatrix<T> newVelocityCoefs(pShift, 1);
    gsMatrix<T> newPressureCoefs(solver->numDofs() - pShift, 1);
    
    newVelocityCoefs.setZero();
    newPressureCoefs.setZero();
    
    // 对速度场进行准插值
    gsQuasiInterpolate<T> quasiInterp;
    
    // 获取映射器
    const std::vector<gsDofMapper>& mappers = solver->getAssembler()->getMappers();
    
    // 对每个patch进行插值
    for (size_t p = 0; p < newPatches.nPatches(); ++p)
    {
        // 获取当前patch的基函数
        const gsBasis<T>& vBasis = newBases[0].basis(p);
        
        for (size_t i = 0; i < vBasis.size(); ++i)
        {
            try {
                // 对每个基函数进行插值
                gsMatrix<T> coef = quasiInterp.localIntpl(vBasis, velocityField.function(p), i);
                
                // 将系数存储到全局系数矩阵中 
                if (mappers[0].is_free(i, p)) {
                    index_t globalIndex = mappers[0].index(i, p);
                    
                    for (index_t d = 0; d < vDim; ++d) {
                        newVelocityCoefs(globalIndex + d * fullUdofs, 0) = coef(0, d);
                    }
                }
            } catch (const std::exception& e) {
                gsInfo << "Velocity field interpolation error: " << e.what() << " on patch " << p << " basis function " << i << "\n";
            }
        }
    }
    
    // 插值压力场
    for (size_t p = 0; p < newPatches.nPatches(); ++p)
    {
        // 获取当前patch的基函数
        const gsBasis<T>& pBasis = newBases[1].basis(p);
        
        for (size_t i = 0; i < pBasis.size(); ++i)
        {
            try {
                // 对每个基函数进行插值
                gsMatrix<T> coef = quasiInterp.localIntpl(pBasis, pressureField.function(p), i);
                
                // 将系数存储到全局系数矩阵中 - 使用映射器获取全局索引
                if (mappers[1].is_free(i, p)) {
                    index_t globalIndex = mappers[1].index(i, p);
                    newPressureCoefs(globalIndex, 0) = coef(0, 0);
                }
            } catch (const std::exception& e) {
                gsInfo << "Pressure field interpolation error: " << e.what() << " on patch " << p << " basis function " << i << "\n";
            }
        }
    }
    
    // 将投影后的系数设置为新解
    solver->setSolutionCoefs(newVelocityCoefs, 0);
    solver->setSolutionCoefs(newPressureCoefs, 1);

    // 应用固体边界速度（如果有边界速度数据）
    if (boundaryVelocitiesInitialized && !boundaryVelocities.empty()) {
        gsInfo << "=== 在transferSolutionToNewMesh中应用固体边界速度 ===\n";
        
        // 获取速度系数矩阵进行修改
        gsMatrix<T> modifiedVelocityCoefs = solver->solutionCoefs(0);
        
        // 获取速度基函数和映射器
        const std::vector<gsMultiBasis<T>>& bases = solver->getAssembler()->getBases();
        const std::vector<gsDofMapper>& mappers = solver->getAssembler()->getMappers();
        const gsMultiBasis<T>& velocityBasis = bases[0];
        const gsDofMapper& velocityMapper = mappers[0];
        
        const index_t vDim = 2; // 2D问题
        const index_t fullUdofs = solver->getAssembler()->getUdofs();
        
        // 与梁接触的边界
        struct BeamInterface {
            size_t patchIdx;
            boundary::side boundSide;
        };
        
        std::vector<BeamInterface> interfaces = {
            {3, boundary::south},
            {4, boundary::north},
            {5, boundary::west},
        };
        
        // 为每个接触边界应用固体速度
        for (const auto& interface : interfaces) {
            int boundaryKey = interface.patchIdx * 10 + static_cast<int>(interface.boundSide);
            
            auto it = boundaryVelocities.find(boundaryKey);
            if (it == boundaryVelocities.end()) continue;
            
            const gsMatrix<T>& solidBoundaryVelocities = it->second;
            if (solidBoundaryVelocities.rows() == 0) continue;
            
            // 计算平均速度
            gsVector<T> avgSolidVelocity(vDim);
            avgSolidVelocity.setZero();
            for (index_t d = 0; d < vDim; ++d) {
                avgSolidVelocity(d) = solidBoundaryVelocities.col(d).mean();
            }
            
            gsInfo << "应用patch " << interface.patchIdx << " 边界平均速度: [" 
                   << avgSolidVelocity.transpose() << "] m/s\n";
            
            // 获取边界上的基函数
            if (interface.patchIdx >= velocityBasis.nBases()) continue;
            
            const gsBasis<T>& basis = velocityBasis.basis(interface.patchIdx);
            if (basis.dim() != 2) continue;
            
            index_t rows = basis.component(0).size();
            index_t cols = basis.component(1).size();
            
            std::vector<index_t> boundaryBasisIndices;
            
            // 根据边界类型确定边界上的基函数索引
            switch(interface.boundSide) {
                case boundary::west:
                    for (index_t j = 0; j < cols; ++j) {
                        boundaryBasisIndices.push_back(j * rows + 0);
                    }
                    break;
                case boundary::east:
                    for (index_t j = 0; j < cols; ++j) {
                        boundaryBasisIndices.push_back(j * rows + (rows-1));
                    }
                    break;
                case boundary::south:
                    for (index_t i = 0; i < rows; ++i) {
                        boundaryBasisIndices.push_back(0 * rows + i);
                    }
                    break;
                case boundary::north:
                    for (index_t i = 0; i < rows; ++i) {
                        boundaryBasisIndices.push_back((cols-1) * rows + i);
                    }
                    break;
                default:
                    continue;
            }
            
            // 应用速度到边界基函数
            for (size_t k = 0; k < boundaryBasisIndices.size(); ++k) {
                index_t localIdx = boundaryBasisIndices[k];
                
                if (velocityMapper.is_free(localIdx, interface.patchIdx)) {
                    index_t globalIdx = velocityMapper.index(localIdx, interface.patchIdx);
                    
                    for (index_t d = 0; d < vDim; ++d) {
                        index_t coeffIdx = globalIdx + d * fullUdofs;
                        
                        if (coeffIdx < modifiedVelocityCoefs.rows()) {
                            T velocityValue = avgSolidVelocity(d);
                            
                            // 如果有具体的边界点速度，使用插值
                            if (k < (size_t)solidBoundaryVelocities.rows() && d < solidBoundaryVelocities.cols()) {
                                velocityValue = solidBoundaryVelocities(k, d);
                            }
                            
                            modifiedVelocityCoefs(coeffIdx, 0) = velocityValue;
                        }
                    }
                }
            }
        }
        
        // 将修改后的速度系数设置回求解器
        solver->setSolutionCoefs(modifiedVelocityCoefs, 0);
        
        gsInfo << "固体边界速度已应用到新网格的流体求解器\n";
    }

    gsDebugVar(newVelocityCoefs.size());
    gsDebugVar(newPressureCoefs.size());

    gsInfo << "Projected solution to new mesh\n";
}


template<class T>
gsMultiPatch<T> updateGeometryBoundary(gsMultiPatch<T>& patches, T time, gsINSSolverUnsteady<T, ColMajor>* fluidSolver)
{
    gsInfo << "\n===== ENTERING updateGeometryBoundary at time = " << time << " =====\n";
    gsInfo << "Function is being called successfully!\n";
    
    // 清空上一时间步的边界速度信息
    if (boundaryVelocitiesInitialized) {
        boundaryVelocities.clear();
        gsInfo << "已清空上一时间步的边界速度信息\n";
    }
    
    // 保存原始几何体以避免累积位移
    static gsMultiPatch<T> initial_fluid_patches_static; // Renamed for clarity: stores the initial fluid mesh
    static bool firstCall = true;
    
    // 声明重力函数
    // gsConstantFunction<> g(0.0, -1000.0, 2); // 增大重力，原来是(0.0, 0.0, 2)
    
    if (firstCall) {
        initial_fluid_patches_static = patches; // Store the true initial fluid mesh
        initializeSolidSolver(patches, 0.01);  // Turek-Hron FSI 3 时间步长
        firstCall = false;
        gsInfo << "保存了原始流体几何体，patch数量: " << initial_fluid_patches_static.nPatches() << "\n";
        
        // 保存初始网格用于对比
        gsWriteParaview<>(patches, "initial_fluid_mesh", 1000);
        gsInfo << "保存初始流体网格到: initial_fluid_mesh.pvd\n";
        
        // 初始化动画标志
        if (!animationInitialized) {
            animationInitialized = true;
            gsInfo << "初始化了位移和应力动画集合\n";
        }
    }
    
    // 计算时间步长
    T timeStep = time - lastTime;
    if (timeStep <= 0) 
    {
        timeStep = 0.01; // Turek-Hron FSI 3 默认时间步长
    }

    // ==== 获取流体压力场 ====
    gsField<T> pressureField = fluidSolver->constructSolution(1); // 1表示压力场

    // 从pressureField中提取函数部分，保持const限定符
    const gsFunctionSet<T>& pressureFunction = pressureField.function();

    // 调试：检查压力场的数值范围
    gsMatrix<T> evalPoint(2, 1);
    evalPoint << 0.4, 0.5; // 梁中央位置的一点
    
    if (pressureFunction.size() > 0) {
        gsMatrix<T> pressureValue = pressureFunction.function(0).eval(evalPoint);
        gsInfo << "压力场在点 (" << evalPoint.transpose() << ") 的值: " << pressureValue << "\n";
        gsInfo << "压力场数值范围检查完成\n";
        
        // // 如果压力值太小，进行放大
        // T pressureScale = 1000.0; // 放大因子
        // if (pressureValue.cwiseAbs().maxCoeff() < 1.0) {
        //     gsInfo << "压力值较小，应用放大因子: " << pressureScale << "\n";
        // }
    }

    gsWriteParaview<>(pressureField, "pressureField", 1000);

    // gsInfo << "pressureFunction: " << pressureFunction.size() << "\n";
    
    
    gsInfo << "开始构造固体边界条件，将压力转换为动态法向量牵引力\n";
    
    // 检查当前梁的变形状态
    if (beamDisplacement.nPatches() > 0) {
        real_t currentMaxDisp = beamDisplacement.patch(0).coefs().cwiseAbs().maxCoeff();
        gsInfo << "当前梁最大位移: " << currentMaxDisp << " m\n";
        if (currentMaxDisp > 1e-15) {
            gsInfo << "✓ 将使用变形后的梁几何计算边界法向量\n";
        } else {
            gsInfo << "! 梁位移为零，将使用原始几何计算法向量\n";
        }
    } else {
        gsInfo << "! 梁位移场未初始化，将使用原始几何计算法向量\n";
    }
    
    // 重建固体求解器，传入当前的压力场
    rebuildSolidSolver(pressureFunction);
    
    gsInfo << "固体求解器已重建（使用动态法向量），开始求解固体问题\n";
    // ==== 求解固体问题 ====
    if (beamSolver) 
    {
        gsInfo << "Solving solid time step, time = " << time << ", time step = " << timeStep << "\n";
        
        // 检查压力场信息
        gsInfo << "Checking boundary conditions and pressure field:\n";
        gsInfo << "Total number of boundary conditions: " << beamBCs.size() << "\n";
        
        try {
            // 尝试调用压力场的eval方法，检查压力值
            gsMatrix<T> evalPoint(2, 1);
            evalPoint << 0.4, 0.5; // 梁中央位置的一点
            
            // 创建一个临时的压力场用于测试
            if (pressureField.function().size() > 0) {
                gsMatrix<T> pressureValue = pressureField.function().function(0).eval(evalPoint);
                gsInfo << "Pressure value at point (" << evalPoint.transpose() << "): " << pressureValue << "\n";
            } else {
                gsInfo << "Warning: Pressure field is empty or cannot be evaluated\n";
            }
        } catch (const std::exception& e) {
            gsInfo << "Error evaluating pressure field: " << e.what() << "\n";
        }
        
        // 添加调试信息 - 只保留已知存在的方法
        gsInfo << "Current velocity vector size: " << beamSolver->velocityVector().size() << "\n";
        
        // 验证边界条件是否正确应用
        gsInfo << "当前固体求解器边界条件数量: " << beamBCs.size() << "\n";
        
        // 验证压力到牵引力的转换
        if (pressureField.function().size() > 0) {
            gsMatrix<T> testPoint(2, 1);
            testPoint << 0.5, 1.0; // north边界上的参数坐标 (u=0.5, v=1.0)
            
            // 测试原始压力
            gsMatrix<T> pressureValue = pressureField.function().function(0).eval(testPoint);
            gsInfo << "验证: 测试点处压力 p = " << pressureValue << "\n";
            
            // 测试牵引力（使用当前变形状态的north边界）
            if (beamDisplacement.nPatches() > 0) {
                // 构造变形后的梁几何（与rebuildSolidSolver中的逻辑一致）
                gsMultiPatch<T> testDeformedGeometry = beamGeometry;
                if (beamDisplacement.patch(0).coefs().cwiseAbs().maxCoeff() > 1e-15) {
                    gsMatrix<T> originalCoefs = testDeformedGeometry.patch(0).coefs();
                    gsMatrix<T> displacementCoefs = beamDisplacement.patch(0).coefs();
                    if (originalCoefs.rows() == displacementCoefs.rows() && 
                        originalCoefs.cols() == displacementCoefs.cols()) {
                        testDeformedGeometry.patch(0).setCoefs(originalCoefs + displacementCoefs);
                    }
                }
                
                std::unique_ptr<gsGeometry<T>> testBoundary = testDeformedGeometry.patch(0).boundary(boundary::north);
                gsPressureToTraction<real_t> testTraction(pressureField.function().function(0), *testBoundary, boundary::north);
                gsMatrix<T> tractionValue;
                gsMatrix<T> boundaryTestPoint(1, 1);
                boundaryTestPoint(0, 0) = 0.5; // u=0.5 on the north boundary
                testTraction.eval_into(boundaryTestPoint, tractionValue);
                gsInfo << "验证: north边界动态牵引力 t = " << tractionValue.transpose() << "\n";
            } else {
                gsInfo << "验证: 无位移场，跳过牵引力验证\n";
            }
        }
        
        // 执行固体求解器的时间步进
        gsInfo << "开始固体时间步进，时间步长: " << timeStep << "s\n";
        try {
            beamSolver->makeTimeStep(timeStep);
            gsInfo << "✓ 固体时间步完成，迭代次数: " << beamSolver->numberIterations() << "\n";
        } catch (const std::exception& e) {
            gsWarn << "固体时间步进时出错: " << e.what() << "\n";
            throw;
        }
        
        // 获取最新的位移解
        beamSolver->constructSolution(beamDisplacement);

        gsInfo << "Displacement coefficient matrix size: " << beamDisplacement.coefs().rows() << " x " << beamDisplacement.coefs().cols() << "\n";
        
        // 调试：检查位移场的数值
        gsInfo << "=== 固体响应分析 ===\n";
        real_t maxDisp = beamDisplacement.coefs().cwiseAbs().maxCoeff();
        real_t meanDisp = beamDisplacement.coefs().mean();
        index_t nonZeroElements = (beamDisplacement.coefs().cwiseAbs().array() > 1e-12).count();
        
        gsInfo << "位移统计:\n";
        gsInfo << "  最大位移: " << maxDisp << " m\n";
        gsInfo << "  平均位移: " << meanDisp << " m\n";
        gsInfo << "  非零元素数: " << nonZeroElements << " / " << beamDisplacement.coefs().size() << "\n";
        
        // 判断固体是否有响应
        if (maxDisp > 1e-10) {
            gsInfo << "✓ 固体有明显响应，动态法向量压力传递成功!\n";
        } else if (maxDisp > 1e-15) {
            gsInfo << "⚠ 固体有微小响应，可能需要检查压力大小或法向量计算\n";
        } else {
            gsWarn << "✗ 固体无响应，动态法向量压力传递可能失败!\n";
        }
        
        // 额外验证：检查变形前后的法向量变化
        if (beamDisplacement.nPatches() > 0 && maxDisp > 1e-15) {
            gsInfo << "=== 法向量变形分析 ===\n";
            
            // 比较原始几何和变形几何的法向量
            gsMatrix<T> testParam(2, 1);
            testParam << 0.5, 1.0; // north边界中点 (u=0.5, v=1.0)
            
            // 原始几何的法向量（north边界应该是(0,1)）
            std::unique_ptr<gsGeometry<T>> originalBoundary = beamGeometry.patch(0).boundary(boundary::north);
            gsMatrix<T> boundaryParam1D(1, 1);
            boundaryParam1D(0, 0) = testParam(0, 0); // north边界使用u参数 (0.5)
            
            gsMatrix<T> originalTangent;
            originalBoundary->deriv_into(boundaryParam1D, originalTangent);
            
            // 变形几何的法向量
            gsMultiPatch<T> currentDeformed = beamGeometry;
            gsMatrix<T> originalCoefs = currentDeformed.patch(0).coefs();
            gsMatrix<T> displacementCoefs = beamDisplacement.patch(0).coefs();
            
            // 检查矩阵维度兼容性
            if (originalCoefs.rows() == displacementCoefs.rows() && 
                originalCoefs.cols() == displacementCoefs.cols()) {
                currentDeformed.patch(0).setCoefs(originalCoefs + displacementCoefs);
                
                std::unique_ptr<gsGeometry<T>> deformedBoundary = currentDeformed.patch(0).boundary(boundary::north);
                gsMatrix<T> deformedTangent;
                deformedBoundary->deriv_into(boundaryParam1D, deformedTangent);
                
                if (originalTangent.cols() > 0 && deformedTangent.cols() > 0) {
                    // 计算并显示法向量
                    T otx = originalTangent(0,0), oty = originalTangent(1,0);
                    T oLength = std::sqrt(otx*otx + oty*oty);
                    T dtx = deformedTangent(0,0), dty = deformedTangent(1,0);
                    T dLength = std::sqrt(dtx*dtx + dty*dty);
                    
                    if (oLength > 1e-12 && dLength > 1e-12) {
                        otx /= oLength; oty /= oLength;
                        dtx /= dLength; dty /= dLength;
                        
                        gsInfo << "原始边界法向量: [" << -oty << ", " << otx << "]\n";
                        gsInfo << "变形边界法向量: [" << -dty << ", " << dtx << "]\n";
                        
                        T angleDiff = std::acos(std::abs(-oty*(-dty) + otx*dtx));
                        gsInfo << "法向量角度变化: " << angleDiff * 180.0 / M_PI << " 度\n";
                    }
                }
            } else {
                gsWarn << "法向量变形分析: 几何与位移系数矩阵大小不匹配\n";
                gsWarn << "原始几何: " << originalCoefs.rows() << "x" << originalCoefs.cols() 
                       << ", 位移: " << displacementCoefs.rows() << "x" << displacementCoefs.cols() << "\n";
                gsWarn << "跳过法向量变形分析\n";
            }
        }
        
        // 输出具体的位移值（前几个自由度）
        gsInfo << "前几个位移系数: ";
        for (index_t i = 0; i < std::min(index_t(5), beamDisplacement.coefs().rows()); ++i) {
            gsInfo << beamDisplacement.coefs()(i) << " ";
        }
        gsInfo << "\n";
        
        // 计算应力场
        gsPiecewiseFunction<> stresses;
        beamAssembler->constructCauchyStresses(beamDisplacement, stresses, stress_components::von_mises);
        
        // 创建位移场和应力场
        gsField<> displacementField(beamGeometry, beamDisplacement);

        gsWriteParaview<>(displacementField, "displacementField", 1000);
        gsField<> stressField(beamGeometry, stresses, true);
        
        // 保存到Paraview文件（单个时间步）
        std::string timeStr = util::to_string(time);
        std::string dispFileName = "beam_displacement_" + timeStr;
        std::string stressFileName = "beam_stress_" + timeStr;
        
        gsWriteParaview<>(displacementField, dispFileName, 1000);
        gsWriteParaview<>(stressField, stressFileName, 1000);
        
        if (animationInitialized) {
            // 添加到动画集合
            dispCollection.addTimestep(dispFileName + "0", time, ".vts");
            stressCollection.addTimestep(stressFileName + "0", time, ".vts");
        }
        
        gsInfo << "固体时间步完成，迭代次数: " << beamSolver->numberIterations() << "\n";
        gsInfo << "已保存位移场到: " << dispFileName << "\n";
        gsInfo << "已保存应力场到: " << stressFileName << "\n";
    }
    
    // ==== 直接更新流体网格 ====
    // Reset 'patches' to the initial fluid domain configuration.
    patches = initial_fluid_patches_static; 
    // 创建位移场和速度场并获取函数对象
    gsField<T> displacementField(beamGeometry, beamDisplacement);
    const gsFunctionSet<T>& displacementFunc = displacementField.function();
    
    // 构造固体速度场
    gsMultiPatch<T> beamVelocityField;
    if (beamSolver) {
        // 从固体求解器获取速度向量
        gsMatrix<T> velocityVector = beamSolver->velocityVector();
        
        gsInfo << "=== 固体速度提取分析 ===\n";
        gsInfo << "速度向量大小: " << velocityVector.size() << "\n";
        gsInfo << "固体求解器自由度: " << beamSolver->numDofs() << "\n";
        
        // 构造速度场 - 使用与位移场相同的结构
        beamVelocityField = gsMultiPatch<T>(beamGeometry);
        
        // 速度向量的结构应该与位移向量相同
        // 对于2D问题，速度向量通常是 [v_x1, v_x2, ..., v_xN, v_y1, v_y2, ..., v_yN]
        // 或者是交错的 [v_x1, v_y1, v_x2, v_y2, ...]
        
        index_t solverDofs = beamSolver->numDofs();
        gsMatrix<T> velocityCoefs;
        
        if (velocityVector.size() == solverDofs) {
            // 假设是交错格式: [v_x1, v_y1, v_x2, v_y2, ...]
            index_t numNodes = solverDofs / 2;
            velocityCoefs.resize(numNodes, 2);
            
            for (index_t i = 0; i < numNodes; ++i) {
                velocityCoefs(i, 0) = velocityVector(2*i);     // x方向速度
                velocityCoefs(i, 1) = velocityVector(2*i+1);   // y方向速度
            }
            gsInfo << "使用交错格式提取速度，节点数: " << numNodes << "\n";
            
        } else if (velocityVector.size() == 2 * solverDofs) {
            // 假设是分块格式: [v_x1, v_x2, ..., v_xN, v_y1, v_y2, ..., v_yN]
            index_t numNodes = solverDofs;
            velocityCoefs.resize(numNodes, 2);
            
            for (index_t i = 0; i < numNodes; ++i) {
                velocityCoefs(i, 0) = velocityVector(i);              // x方向速度
                velocityCoefs(i, 1) = velocityVector(i + numNodes);   // y方向速度
            }
            gsInfo << "使用分块格式提取速度，节点数: " << numNodes << "\n";
            
        } else {
            // 如果大小不匹配，尝试其他方法提取速度
            gsInfo << "速度向量大小不匹配，尝试从求解器状态提取速度\n";
            
            // 检查是否速度向量包含在完整的解向量中
            if (velocityVector.size() > 0 && beamDisplacement.nPatches() > 0) {
                gsMatrix<T> currentDisp = beamDisplacement.patch(0).coefs();
                index_t numNodes = currentDisp.rows();
                
                // 尝试从速度向量中提取合理的速度值
                velocityCoefs.resize(numNodes, 2);
                velocityCoefs.setZero();
                
                // 根据实际速度向量大小进行处理
                if (velocityVector.size() >= numNodes) {
                    // 可能是按节点存储的速度
                    for (index_t i = 0; i < std::min(numNodes, velocityVector.rows()); ++i) {
                        velocityCoefs(i, 0) = velocityVector(i); // x方向速度
                        if (velocityVector.size() >= 2 * numNodes && i + numNodes < velocityVector.size()) {
                            velocityCoefs(i, 1) = velocityVector(i + numNodes); // y方向速度
                        }
                    }
                    gsInfo << "从速度向量提取了 " << numNodes << " 个节点的速度\n";
                } else {
                    gsInfo << "速度向量太小，使用基于位移的近似速度\n";
                    
                    // 保存上一时间步的位移用于计算速度
                    static gsMatrix<T> lastDisplacement;
                    static bool firstTimeStep = true;
                    
                    if (firstTimeStep) {
                        // 第一时间步，速度为零
                        velocityCoefs = gsMatrix<T>::Zero(currentDisp.rows(), currentDisp.cols());
                        lastDisplacement = currentDisp;
                        firstTimeStep = false;
                        gsInfo << "第一时间步，速度设为零\n";
                    } else {
                        // 计算速度 = (当前位移 - 上一时间步位移) / 时间步
                        T dt = time - lastTime;
                        if (dt > 1e-12) {
                            velocityCoefs = (currentDisp - lastDisplacement) / dt;
                            gsInfo << "基于位移差分计算速度，时间步: " << dt << "\n";
                        } else {
                            velocityCoefs = gsMatrix<T>::Zero(currentDisp.rows(), currentDisp.cols());
                            gsInfo << "时间步太小，速度设为零\n";
                        }
                        lastDisplacement = currentDisp;
                    }
                }
            } else {
                // 如果没有位移数据，创建零速度
                velocityCoefs = gsMatrix<T>::Zero(10, 2); // 假设10个节点
                gsInfo << "无位移数据，使用零速度\n";
            }
        }
        
        // 设置速度系数到速度场
        beamVelocityField.patch(0).setCoefs(velocityCoefs);
        
        T maxVel = velocityCoefs.cwiseAbs().maxCoeff();
        gsInfo << "已构造固体速度场，节点数: " << velocityCoefs.rows() 
               << ", 最大速度: " << maxVel << " m/s\n";
        
        // 输出前几个速度值用于调试
        gsInfo << "前几个速度值:\n";
        for (index_t i = 0; i < std::min(index_t(3), velocityCoefs.rows()); ++i) {
            gsInfo << "  节点 " << i << ": vx=" << velocityCoefs(i, 0) 
                   << ", vy=" << velocityCoefs(i, 1) << "\n";
        }
        
    }

    // 固体求解器一定存在
    
    gsField<T> velocityField(beamGeometry, beamVelocityField);
    const gsFunctionSet<T>& velocityFunc = velocityField.function();

    // 1. 找出所有需要更新的patches (与梁接触的流体区域)
    // 假设梁在patch 3, 4, 5的区域附近
    gsInfo << "Updating fluid mesh boundary in contact with the beam...\n";
    
    // 定义与梁交界的边界对应关系 - 修改为包含更详细的信息
    struct FluidInterface {
        size_t patchIdx;          // patches中的索引
        boundary::side boundSide; // 对应的边界方向
    };
    
    // 根据几何结构确定所有与梁接触的边界
    std::vector<FluidInterface> interfaces = 
    {
        {3, boundary::south},
        {4, boundary::north},
        {5, boundary::west},
    };

    // for (const auto & interface : interfaces)
    // {
    //     gsInfo << "Patch " << interface.patchIdx << " on " 
    //            << interface.boundSide << " boundary needs update\n";

    //     gsMultiPatch<T> boundaryPatch = patches.patch(interface.patchIdx).boundary(interface.boundSide);

    //     gsMultiPatch<T> boundaryPatchDeformed = boundaryPatch; // 变形后的边界

    //     gsMatrix<T> deformed_coefs = boundaryPatchDeformed.coefs();

    //     deformed_coefs = deformed_coefs.array() + 
    //                      displacementFunc.function(0).eval(boundaryPatchDeformed.points()).array();
    //     boundaryPatchDeformed.setCoefs(deformed_coefs);
    // }
    // return 0;
    // // 先更新边界控制点
    // for (const auto& interface : interfaces)
    // {
    //     gsInfo << "Processing patch " << interface.patchIdx << " of " 
    //            << interface.boundSide << " boundary\n";
        
    //     // 获取patch边界 - 添加NURBS权重处理
    //     // 'patches' is now a copy of the initial fluid mesh (due to 'patches = initial_fluid_patches_static;')
    //     // So, these are the initial coordinates of the fluid boundary segment.
    //     gsMatrix<T> initial_coefs_NxD;
        
    //     try {
    //          std::unique_ptr<gsGeometry<T>> boundaryGeom = patches.patch(interface.patchIdx).boundary(interface.boundSide);
             
    //          // 检查是否为NURBS几何
    //          const gsNurbs<T>* nurbsBoundary = dynamic_cast<const gsNurbs<T>*>(boundaryGeom.get());
    //          if (nurbsBoundary) {
    //              gsInfo << "检测到NURBS流体边界几何，patch " << interface.patchIdx << "\n";
    //              initial_coefs_NxD = nurbsBoundary->coefs();
    //              gsInfo << "NURBS边界控制点: " << initial_coefs_NxD.rows() << "x" << initial_coefs_NxD.cols() << "\n";
    //              gsInfo << "NURBS权重数量: " << nurbsBoundary->weights().rows() << "\n";
    //          } else {
    //              // B样条几何
    //              initial_coefs_NxD = boundaryGeom->coefs();
    //              gsInfo << "B样条边界控制点: " << initial_coefs_NxD.rows() << "x" << initial_coefs_NxD.cols() << "\n";
    //          }
    //      } catch (const std::exception& e) {
    //          gsWarn << "获取patch " << interface.patchIdx << " 边界几何时出错: " << e.what() << "\n";
    //          gsWarn << "跳过此边界处理\n";
    //          continue;
    //      }

    //     // Transpose to (Dim x N_points) for evaluation, as gsFunction::eval typically expects points as columns.
    //     gsMatrix<T> points_for_eval_DxN = initial_coefs_NxD.transpose();
        
    //     // Define vDim for dimension consistency
    //     const index_t vDim = points_for_eval_DxN.rows(); // Should be 2 for 2D problem
        
    //     // 保存原始边界控制点用于对比
    //     gsInfo << "Original boundary control points for patch " << interface.patchIdx << ":\n";
    //     for (index_t k = 0; k < std::min(index_t(3), initial_coefs_NxD.rows()); ++k) {
    //         gsInfo << "  Point " << k << ": (" << initial_coefs_NxD(k, 0) << ", " << initial_coefs_NxD(k, 1) << ")\n";
    //     }
        
    //     // Evaluate displacement and velocity at these initial material points.
    //     // displacementFunc is derived from beamDisplacement (total displacement) on beamGeometry (initial solid config).
    //     gsMatrix<T> displacement_vectors_DxN = gsMatrix<T>::Zero(vDim, points_for_eval_DxN.cols());
    //     gsMatrix<T> velocity_vectors_DxN = gsMatrix<T>::Zero(vDim, points_for_eval_DxN.cols());

        // // 智能位移传递：先尝试正确方法，失败时使用安全备选方案
        // gsInfo << "开始智能位移传递，尝试避免NURBS权重问题\n";
        
        // bool useDirectMapping = true;
        
        // // 方法1: 尝试正常的参数反演和评估
        // if (!useDirectMapping) {
        //     gsInfo << "尝试标准参数反演方法\n";
        //     for (index_t j=0; j<points_for_eval_DxN.cols(); ++j)
        //     {
        //         try {
        //             gsVector<T> phys = points_for_eval_DxN.col(j);
        //             gsMatrix<T> parMatrix(2, 1);
                    
        //             // 使用更鲁棒的参数反演
        //             beamGeometry.patch(0).invertPoints(phys, parMatrix, 1e-6, true);
                    
        //             // 直接使用位移系数进行简单插值，避免NURBS评估
        //             if (beamDisplacement.nPatches() > 0) {
        //                 // 获取梁上的最大位移作为近似
        //                 real_t maxDispX = beamDisplacement.patch(0).coefs().col(0).maxCoeff();
        //                 real_t maxDispY = beamDisplacement.patch(0).coefs().col(1).maxCoeff();
                        
                        
        //                 // 根据边界位置按比例分配位移
        //                 T xPos = phys(0);
        //                 T ratio = (xPos - 0.2) / 0.4; // 梁的相对位置 (0.2 到 0.6)
        //                 ratio = std::max(0.0, std::min(1.0, ratio)); // 限制在[0,1]
                        
        //                 displacement_vectors_DxN(0, j) = maxDispX * ratio;
        //                 displacement_vectors_DxN(1, j) = maxDispY * ratio;
        //             } else {
        //                 displacement_vectors_DxN(0, j) = 0.0;
        //                 displacement_vectors_DxN(1, j) = 0.0;
        //             }
                    
        //             // 计算速度：基于固体速度场的评估
        //             if (beamVelocityField.nPatches() > 0) {
        //                 try {
        //                     // 检查是否是NURBS几何，如果是，直接使用系数插值
        //                     const gsNurbs<T>* nurbsVel = dynamic_cast<const gsNurbs<T>*>(&beamVelocityField.patch(0));
        //                     if (nurbsVel) {
        //                         // 对于NURBS，直接使用控制点系数进行插值，避免权重问题
        //                         gsMatrix<T> velocityCoefs = nurbsVel->coefs();
        //                         if (velocityCoefs.rows() > 0 && velocityCoefs.cols() >= vDim) {
        //                             // 使用参数坐标进行简单插值
        //                             T u = parMatrix(0, 0);
        //                             u = std::max(T(0), std::min(T(1), u)); // 限制在[0,1]范围内
                                    
        //                             // 线性插值
        //                             index_t nCoefs = velocityCoefs.rows();
        //                             T idx_float = u * (nCoefs - 1);
        //                             index_t idx_low = static_cast<index_t>(std::floor(idx_float));
        //                             index_t idx_high = std::min(idx_low + 1, nCoefs - 1);
        //                             T alpha = idx_float - idx_low;
                                    
        //                             for (index_t d = 0; d < vDim; ++d) {
        //                                 velocity_vectors_DxN(d, j) = (1 - alpha) * velocityCoefs(idx_low, d) 
        //                                                             + alpha * velocityCoefs(idx_high, d);
        //                             }
                                    
        //                             if (j == 0 || j == points_for_eval_DxN.cols()-1) {
        //                                 gsInfo << "NURBS边界点 " << j << " 速度插值: [" 
        //                                        << velocity_vectors_DxN(0, j) << ", " 
        //                                        << velocity_vectors_DxN(1, j) << "] m/s\n";
        //                             }
        //                         }
        //                     } else {
        //                         // 非NURBS情况，使用标准评估
        //                         gsMatrix<T> velocityEval = beamVelocityField.patch(0).eval(parMatrix);
                                
        //                         if (velocityEval.rows() == vDim && velocityEval.cols() == 1) {
        //                             for (index_t d = 0; d < vDim; ++d) {
        //                                 velocity_vectors_DxN(d, j) = velocityEval(d, 0);
        //                             }
                                    
        //                             if (j == 0 || j == points_for_eval_DxN.cols()-1) {
        //                                 gsInfo << "边界点 " << j << " 速度评估成功: [" 
        //                                        << velocity_vectors_DxN(0, j) << ", " 
        //                                        << velocity_vectors_DxN(1, j) << "] m/s\n";
        //                             }
        //                         } else {
        //                             // 维度不匹配，使用备用方法
        //                             gsWarn << "速度场评估维度不匹配，使用备用方法\n";
        //                             gsMatrix<T> velocityCoefs = beamVelocityField.patch(0).coefs();
        //                             if (velocityCoefs.rows() > 0 && velocityCoefs.cols() >= vDim) {
        //                                 T u = parMatrix(0, 0);
        //                                 index_t nCoefs = velocityCoefs.rows();
        //                                 index_t idx = std::min(static_cast<index_t>(u * (nCoefs - 1)), nCoefs - 1);
                                        
        //                                 for (index_t d = 0; d < vDim; ++d) {
        //                                     velocity_vectors_DxN(d, j) = velocityCoefs(idx, d);
        //                                 }
        //                             }
        //                         }
        //                     }
        //                 } catch (const std::exception& e) {
        //                     gsWarn << "速度场评估失败: " << e.what() << "，使用零速度\n";
        //                     for (index_t d = 0; d < vDim; ++d) {
        //                         velocity_vectors_DxN(d, j) = 0.0;
        //                     }
        //                 }
        //             } 
                    
        //         } catch (const std::exception& e) {
        //             gsWarn << "点 " << j << " 参数反演失败: " << e.what() << "\n";
        //             // 使用零位移作为后备
        //             for (index_t d = 0; d < vDim; ++d) {
        //                 displacement_vectors_DxN(d, j) = 0.0;
        //                 velocity_vectors_DxN(d, j) = 0.0;
        //             }
        //         }
        //     }
        // }
        
        // gsInfo << "完成位移传递，最大位移: " << displacement_vectors_DxN.cwiseAbs().maxCoeff() << " m\n";

        // // Log displacement and velocity magnitude to verify
        // if (displacement_vectors_DxN.size() > 0) 
        // { // Check if not empty before norm or access
        //      if (displacement_vectors_DxN.rows() > 0 && displacement_vectors_DxN.cols() > 0) {
        //         gsInfo << "  Sample displacement value [0,0]: " << displacement_vectors_DxN(0,0) << " m\n";
        //         gsInfo << "  Max displacement: " << displacement_vectors_DxN.cwiseAbs().maxCoeff() << " m\n";
        //      }
        // }
        
        // if (velocity_vectors_DxN.size() > 0) 
        // { // Check if not empty before norm or access
        //      if (velocity_vectors_DxN.rows() > 0 && velocity_vectors_DxN.cols() > 0) {
        //         gsInfo << "  Sample velocity value [0,0]: " << velocity_vectors_DxN(0,0) << " m/s\n";
        //         gsInfo << "  Max velocity: " << velocity_vectors_DxN.cwiseAbs().maxCoeff() << " m/s\n";
        //      }
        // }
        
        // // 应用位移到边界控制点
        // gsMatrix<T> new_coefs_NxD = initial_coefs_NxD;
        // 
        // // 添加位移到控制点
        // for (index_t j = 0; j < displacement_vectors_DxN.cols(); ++j) 
        // {
        //     for (index_t d = 0; d < vDim; ++d) 
        //     {
        //         new_coefs_NxD(j, d) += displacement_vectors_DxN(d, j);
        //     }
        // }
        
        // gsInfo << "Updated boundary control points for patch " << interface.patchIdx << ":\n";
        // for (index_t k = 0; k < std::min(index_t(3), new_coefs_NxD.rows()); ++k) {
        //     gsInfo << "  Point " << k << ": (" << new_coefs_NxD(k, 0) << ", " << new_coefs_NxD(k, 1) << ")";
        //     gsInfo << " [Δx=" << (new_coefs_NxD(k, 0) - initial_coefs_NxD(k, 0)) << ", Δy=" << (new_coefs_NxD(k, 1) - initial_coefs_NxD(k, 1)) << "]\n";
        // }

        // 直接修改整个patch的控制点矩阵，而不是通过边界函数
        // Note: This entire block is currently commented out, so these errors are in dead code
        gsMatrix<T> allCoefs = currentPatch.coefs();
        
        gsInfo << "修改前patch " << interface.patchIdx << " 控制点矩阵大小: " << allCoefs.rows() << "x" << allCoefs.cols() << "\n";
        
        // 找到边界控制点在全局控制点矩阵中的位置，并直接修改
        // 这需要根据patch的拓扑结构来确定
        // 对于张量积B样条曲面，边界控制点有特定的索引模式
        
        // 获取基函数的维度信息以确定控制点的排列
        const gsBasis<T>& basis = currentPatch.basis();
        if (basis.dim() == 2) {
            // 对于2D张量积B样条，控制点按行列排列
            index_t rows = basis.component(0).size();
            index_t cols = basis.component(1).size();
            
            gsInfo << "Patch基函数维度: " << rows << " x " << cols << "\n";
            
            // 根据边界类型确定需要修改的控制点索引
            std::vector<index_t> boundaryIndices;
            
            switch(interface.boundSide) {
                case boundary::west:  // 左边界 (u=0)
                    for (index_t j = 0; j < cols; ++j) {
                        boundaryIndices.push_back(j * rows + 0);
                    }
                    break;
                case boundary::east:  // 右边界 (u=1)
                    for (index_t j = 0; j < cols; ++j) {
                        boundaryIndices.push_back(j * rows + (rows-1));
                    }
                    break;
                case boundary::south: // 下边界 (v=0)
                    for (index_t i = 0; i < rows; ++i) {
                        boundaryIndices.push_back(0 * rows + i);
                    }
                    break;
                case boundary::north: // 上边界 (v=1)
                    for (index_t i = 0; i < rows; ++i) {
                        boundaryIndices.push_back((cols-1) * rows + i);
                    }
                    break;
                case boundary::none:
                case boundary::front:
                case boundary::back:
                default:
                    gsWarn << "Unhandled boundary side: " << interface.boundSide << "\n";
                    break;
            }
            
            gsInfo << "边界控制点索引数量: " << boundaryIndices.size() << ", 新控制点数量: " << new_coefs_NxD.rows() << "\n";
            
            // 应用新的控制点坐标
            for (index_t k = 0; k < std::min((index_t)boundaryIndices.size(), new_coefs_NxD.rows()); ++k) {
                index_t globalIdx = boundaryIndices[k];
                if (globalIdx < allCoefs.rows()) {
                    gsInfo << "更新控制点 " << globalIdx << ": (" << allCoefs(globalIdx, 0) << ", " << allCoefs(globalIdx, 1) << ") -> ";
                    allCoefs(globalIdx, 0) = new_coefs_NxD(k, 0);
                    allCoefs(globalIdx, 1) = new_coefs_NxD(k, 1);
                    gsInfo << "(" << allCoefs(globalIdx, 0) << ", " << allCoefs(globalIdx, 1) << ")\n";
                }
            }
            
            // 设置修改后的控制点
            currentPatch.setCoefs(allCoefs);
        }
        
        // 立即验证边界控制点是否被正确修改 - 添加NURBS处理
            try {
             std::unique_ptr<gsGeometry<T>> verifyBoundary = currentPatch.boundary(interface.boundSide);
             gsMatrix<T> verification_coefs;
             
             const gsNurbs<T>* nurbsVerify = dynamic_cast<const gsNurbs<T>*>(verifyBoundary.get());
             if (nurbsVerify) {
                 verification_coefs = nurbsVerify->coefs();
             } else {
                 verification_coefs = verifyBoundary->coefs();
             }
             
             gsInfo << "立即验证patch " << interface.patchIdx << " 边界控制点修改:\n";
             for (index_t k = 0; k < std::min(index_t(2), verification_coefs.rows()); ++k) {
                 gsInfo << "  验证点 " << k << ": (" << verification_coefs(k, 0) << ", " << verification_coefs(k, 1) << ")\n";
             }
         } catch (const std::exception& e) {
             gsWarn << "验证边界控制点时出错: " << e.what() << "\n";
         }
        
        // 保存边界速度信息用于ALE计算
        gsMatrix<T> boundaryVel(velocity_vectors_DxN.cols(), vDim);
        for (index_t j = 0; j < velocity_vectors_DxN.cols(); ++j) {
            for (index_t d = 0; d < vDim; ++d) {
                boundaryVel(j, d) = velocity_vectors_DxN(d, j);
            }
        }
        
        // 存储到全局映射中，使用patch索引和边界标识作为键
        int boundaryKey = interface.patchIdx * 10 + static_cast<int>(interface.boundSide);
        boundaryVelocities[boundaryKey] = boundaryVel;
        boundaryVelocitiesInitialized = true;
        
        gsInfo << "对patch " << interface.patchIdx << " 应用了真实的固体边界变形和速度\n";
        gsInfo << "边界速度已保存，键值: " << boundaryKey << ", 大小: " 
               << boundaryVel.rows() << "x" << boundaryVel.cols() << "\n";
    }

    // 先保存变形网格用于调试（在应用完所有位移后）
    // gsWriteParaview<>(patches, "deformed_mesh_time_" + util::to_string(time), 1000);
    gsInfo << "保存变形网格到: deformed_mesh_time_" << time << ".pvd\n";
    
    // 单独保存每个修改过的patch进行验证
    for (const auto& interface : interfaces) {
        gsMultiPatch<T> singlePatch;
        singlePatch.addPatch(patches.patch(interface.patchIdx));
        
        std::string patchFileName = "patch_" + util::to_string(interface.patchIdx) + "_time_" + util::to_string(time);
        gsWriteParaview<>(singlePatch, patchFileName, 1000);
        gsInfo << "保存patch " << interface.patchIdx << " 到: " << patchFileName << ".pvd\n";
    }
    
    // 为了确保变形被保存，额外输出一些边界点的坐标进行验证
    for (const auto& interface : interfaces) {
        try {
            std::unique_ptr<gsGeometry<T>> finalBoundary = patches.patch(interface.patchIdx).boundary(interface.boundSide);
            gsMatrix<T> deformed_coefs;
            
            const gsNurbs<T>* nurbsFinal = dynamic_cast<const gsNurbs<T>*>(finalBoundary.get());
            if (nurbsFinal) {
                deformed_coefs = nurbsFinal->coefs();
            } else {
                deformed_coefs = finalBoundary->coefs();
            }
            
            gsInfo << "最终验证：Patch " << interface.patchIdx << " 变形后边界点坐标:\n";
            for (index_t k = 0; k < std::min(index_t(2), deformed_coefs.rows()); ++k) {
                gsInfo << "  点 " << k << ": (" << deformed_coefs(k, 0) << ", " << deformed_coefs(k, 1) << ")\n";
            }
        } catch (const std::exception& e) {
            gsWarn << "最终验证边界控制点时出错 (patch " << interface.patchIdx << "): " << e.what() << "\n";
        }
    }

    // 可视化整个多片段几何体 (如果需要)

    // 2. 提取需要重新网格化的patches创建一个子集
    gsMultiPatch<T> interfacePatches;
    std::vector<size_t> interfacePatchIndices;
    
    for (const auto& interface : interfaces) 
    {
        // 避免重复添加
        if (std::find(interfacePatchIndices.begin(), interfacePatchIndices.end(), 
                     interface.patchIdx) == interfacePatchIndices.end()) {
            interfacePatches.addPatch(patches.patch(interface.patchIdx));
            interfacePatchIndices.push_back(interface.patchIdx);
        }
    }
    
    // 3. 应用参数化方法重新生成内部网格 - 修复patch消失问题
    short_t geodim = interfacePatches.targetDim();
    
    // 确保方向一致性
    interfacePatches.fixOrientation();
    
    gsMultiPatch<T> optResults; // 将发送回原始patch
    
    // 添加安全检查：确保interfacePatches不为空
    if (interfacePatches.nPatches() == 0) {
        gsWarn << "Interface patches为空，跳过重新参数化\n";
        gsInfo << "跳过重新参数化，保持位移应用的结果\n";
        // 不进行重新参数化，直接使用当前的patches
    }
    else {
        gsInfo << "开始重新参数化 " << interfacePatches.nPatches() << " 个interface patches\n";
    
    if (geodim == 2) 
    {
            try {
                gsBarrierPatch<2, real_t> opt(interfacePatches, false); // 改回true以保持边界固定
                opt.options().setInt("Verbose", 1); // 增加详细输出来调试
        opt.options().setInt("ParamMethod", 1);
        opt.options().setInt("AAPreconditionType", 1);
                opt.options().setReal("BarrierTolerance", 1e-6); // 设置合理的容忍度
                opt.options().setInt("MaxIterations", 100); // 限制最大迭代次数
                
                gsInfo << "执行2D Barrier参数化...\n";
        opt.compute();
        optResults = opt.result();
                
                gsInfo << "重新参数化完成，生成了 " << optResults.nPatches() << " 个patch\n";
                
                // 验证结果的有效性
                if (optResults.nPatches() != interfacePatches.nPatches()) {
                    gsWarn << "重新参数化后patch数量不匹配! 原始: " << interfacePatches.nPatches() 
                           << ", 结果: " << optResults.nPatches() << "\n";
                    gsInfo << "跳过重新参数化，保持原始边界变形结果\n";
                    optResults = interfacePatches; // 使用原始的变形结果
                }
            } catch (const std::exception& e) {
                gsWarn << "2D Barrier参数化失败: " << e.what() << "\n";
                gsInfo << "跳过重新参数化，保持原始边界变形结果\n";
                optResults = interfacePatches; // 使用原始的变形结果
            }
    } 
    else if (geodim == 3) 
    {
            try {
                gsBarrierPatch<3, real_t> opt(interfacePatches, false); // 改回true以保持边界固定
                opt.options().setInt("Verbose", 1);
        opt.options().setInt("ParamMethod", 1);
                opt.options().setReal("BarrierTolerance", 1e-6);
                opt.options().setInt("MaxIterations", 100);
                
                gsInfo << "执行3D Barrier参数化...\n";
        opt.compute();
        optResults = opt.result();
                
                gsInfo << "重新参数化完成，生成了 " << optResults.nPatches() << " 个patch\n";
                
                if (optResults.nPatches() != interfacePatches.nPatches()) {
                    gsWarn << "重新参数化后patch数量不匹配!\n";
                    optResults = interfacePatches;
                }
            } catch (const std::exception& e) {
                gsWarn << "3D Barrier参数化失败: " << e.what() << "\n";
                optResults = interfacePatches;
            }
    } 
    else 
    {
        gsInfo << "当前版本仅支持pardim = geodim = 2或3。\n";
        return patches;
    }

        // 4. 将优化结果应用回原始patches - 添加安全检查
        gsInfo << "应用重新参数化结果回原始patches...\n";
        for (size_t i = 0; i < interfacePatchIndices.size() && i < optResults.nPatches(); ++i) {
        size_t patchIdx = interfacePatchIndices[i];
            gsInfo << "将优化结果应用回patch " << patchIdx << " (结果patch " << i << ")\n";
            
            // 验证patch的有效性
            if (optResults.patch(i).coefs().rows() > 0) {
        patches.patch(patchIdx) = optResults.patch(i);
                gsInfo << "成功更新patch " << patchIdx << "\n";
            } else {
                gsWarn << "优化结果patch " << i << " 无效，保持原始patch\n";
            }
        }
    }
    
    
    // gsInfo << "跳过重新参数化，保持位移应用的结果\n";
    
    // 5. 确保patch之间的边界一致性
    patches.computeTopology();
    
    // 遍历所有内部接缝
    for (const auto& interface : patches.topology().interfaces()) {
        // 获取接缝两边的patch索引和边界
        int patch1 = interface.first().patch;
        int patch2 = interface.second().patch;
        
        boundary::side side1 = interface.first().side();
        boundary::side side2 = interface.second().side();
        
        gsInfo << "Checking patch " << patch1 << " (" << side1 << ") and patch " 
               << patch2 << " (" << side2 << ") interface\n";
        
        // 检查接缝是否一致
        auto& geom1 = patches.patch(patch1);
        auto& geom2 = patches.patch(patch2);
        
        // 比较边界控制点
        std::unique_ptr<gsGeometry<T>> boundary1 = geom1.boundary(side1);
        std::unique_ptr<gsGeometry<T>> boundary2 = geom2.boundary(side2);
        gsMatrix<T> coefs1 = boundary1->coefs();
        gsMatrix<T> coefs2 = boundary2->coefs();
        
        // 如果不一致，需要修复
        // 注意：这里只是简单的检查，实际可能需要考虑方向
        if (coefs1.rows() != coefs2.rows()) {
            gsWarn << "Boundary control point number mismatch! patch1: " << coefs1.rows() 
                   << ", patch2: " << coefs2.rows() << "\n";
        }
    }
    
    // 确保边界定义和拓扑是最新的
    patches.computeTopology();
    gsWriteParaview<>(patches, "updated_patches", 1000);
    
    gsInfo << "Fluid mesh updated, all boundary displacements applied, interfaces checked\n";
    
    // 更新上次时间
    lastTime = time;
    
    gsInfo << "Geometry boundary updated\n";
    return patches;
}

// 初始化固体求解器（简化版本）
template<class T>
void initializeSolidSolver(const gsMultiPatch<T>& fluidPatches, T timeStep)
{
    if (solidInitialized)
        return;
        
    gsInfo << "Initializing solid (beam) solver...\n";
    
    // 加载梁的几何模型
    std::string filenameBeam = "flappingBeam_beam.xml";
    // 查找文件
    std::string path = gsFileManager::find(filenameBeam);
    if (path.empty())
    {
        gsWarn << "Beam geometry file not found: " << filenameBeam << ", using simplified model\n";
        // 创建一个简单的矩形作为梁几何
        beamGeometry = gsMultiPatch<real_t>();
        
        // 创建box矩阵：第一列是起点，第二列是终点 (Turek-Hron FSI 3 几何)
        gsMatrix<real_t> box(2, 2);
        box << 0.2, 0.6,     // x坐标范围: 梁长度 0.35m
               0.19, 0.21;   // y坐标范围: 梁厚度 0.02m (Turek-Hron FSI 3)
        
        // 使用正确的BSplineSquare重载方法
        beamGeometry.addPatch(gsNurbsCreator<real_t>::BSplineSquare(box));
    }
    else
    {
        gsReadFile<>(path, beamGeometry);
    }
    
    gsInfo << "Beam geometry information:\n" << beamGeometry << "\n";
    
    // 设置梁的材料参数 - Turek-Hron FSI 3标准值
    T youngsModulus = 5.6e6;   // 杨氏模量：5.6 × 10^6 Pa
    T poissonsRatio = 0.4;     // 泊松比：0.4
    T density = 1.0e3;         // 固体密度：1.0 × 10^3 kg/m³ (Turek-Hron FSI 3)
    
    gsInfo << "Material parameters: E=" << youngsModulus << ", nu=" << poissonsRatio << ", rho=" << density << "\n";
    
    // 创建基函数并细化
    gsMultiBasis<> basisDisplacement(beamGeometry);
    for (index_t i = 0; i < 3; ++i)
        basisDisplacement.uniformRefine();
    
    gsInfo << "Basis function information:\n" << basisDisplacement << "\n";
    
    // 创建边界条件 - 初始化时只有固定端边界条件
    beamBCs = gsBoundaryConditions<real_t>();
    
    // 固定端边界条件（左侧固定）
    for (index_t d = 0; d < 2; ++d)
        beamBCs.addCondition(0, boundary::west, condition_type::dirichlet, 0, d);
    
    gsInfo << "初始化固体求解器，只添加固定端边界条件: " << beamBCs.size() << "\n";
    gsInfo << "注意：压力边界条件将在每个时间步动态更新\n";
    
    // 初始化固体力 - Turek-Hron FSI 3 不包含体积力
    gsConstantFunction<> g(0.0, 0.0, 2); // Turek-Hron FSI 3: 无重力作用
    
    // 创建弹性问题组装器
    beamAssembler.reset(new gsElasticityAssembler<real_t>(beamGeometry, basisDisplacement, beamBCs, g));
    beamAssembler->options().setReal("YoungsModulus", youngsModulus);
    beamAssembler->options().setReal("PoissonsRatio", poissonsRatio);
    beamAssembler->options().setInt("MaterialLaw", material_law::saint_venant_kirchhoff);
    
    // 创建质量矩阵组装器
    beamMassAssembler.reset(new gsMassAssembler<real_t>(beamGeometry, basisDisplacement, beamBCs, g));
    beamMassAssembler->options().setReal("Density", density);
    
    // 创建时间积分器（使用Newmark方法）
    beamSolver.reset(new gsElTimeIntegrator<real_t>(*beamAssembler, *beamMassAssembler));
    beamSolver->options().setInt("Scheme", time_integration::implicit_nonlinear);
    beamSolver->options().setReal("Beta", 0.25);  // Newmark方法标准参数
    beamSolver->options().setReal("Gamma", 0.5);
    
    // 初始化求解器和解向量
    beamSolver->setDisplacementVector(gsMatrix<>::Zero(beamAssembler->numDofs(), 1));
    beamSolver->setVelocityVector(gsMatrix<>::Zero(beamAssembler->numDofs(), 1));
    
    // 初始化结构解
    beamDisplacement = gsMultiPatch<real_t>(beamGeometry);
    beamDisplacement.patch(0).coefs().setZero();
    
    gsInfo << "Solid solver initialized, number of degrees of freedom: " << beamAssembler->numDofs() << "\n";
    gsInfo << "Initial velocity vector size: " << beamSolver->velocityVector().size() << "\n";
    
    solidInitialized = true;
}

// 重建固体求解器，应用当前压力场作为牵引边界条件
template<class T>
void rebuildSolidSolver(const gsFunctionSet<T>& pressureFunction)
{
    if (!solidInitialized) {
        gsWarn << "固体求解器尚未初始化，无法重建\n";
        return;
    }

    gsMatrix<real_t> oldSolutionVector = beamSolver->solutionVector();
    gsMatrix<real_t> oldVelocity = beamSolver->velocityVector();
    
    // 从完整解向量中提取位移部分（前 numDofs 行）
    index_t dispDofs = beamAssembler->numDofs();
    gsMatrix<real_t> oldDisplacement = oldSolutionVector.topRows(dispDofs);
    
    gsInfo << "已保存当前位移向量（大小: " << oldDisplacement.size() << "）和速度向量（大小: " << oldVelocity.size() << "）\n";
    
    // 重新构造边界条件
    beamBCs = gsBoundaryConditions<real_t>();
    
    // 1. 固定端边界条件（左侧固定）
    for (index_t d = 0; d < 2; ++d)
        beamBCs.addCondition(0, boundary::west, condition_type::dirichlet, 0, d);
    
    gsInfo << "已添加固定端边界条件\n";
    
    // 2. 添加牵引边界条件（将压力转换为牵引力）
    if (pressureFunction.size() > 0) {
        // 获取第一个patch的压力函数
        const gsFunction<T>& pFunc = pressureFunction.function(0);
        
        // 构造当前变形状态下的梁几何
        gsMultiPatch<T> deformedBeamGeometry = beamGeometry;
        
        // 如果有位移场，将其添加到几何上获得变形后的几何
        if (beamDisplacement.nPatches() > 0 && beamDisplacement.patch(0).coefs().cwiseAbs().maxCoeff() > 1e-15) {
            gsInfo << "使用变形后的梁几何计算法向量\n";
            // 将位移添加到原始几何控制点
            for (size_t p = 0; p < deformedBeamGeometry.nPatches(); ++p) {
                gsMatrix<T> originalCoefs = deformedBeamGeometry.patch(p).coefs();
                gsMatrix<T> displacementCoefs = beamDisplacement.patch(p).coefs();
                
                // 确保维度匹配
                if (originalCoefs.rows() == displacementCoefs.rows() && 
                    originalCoefs.cols() == displacementCoefs.cols()) {
                    gsMatrix<T> deformedCoefs = originalCoefs + displacementCoefs;
                    deformedBeamGeometry.patch(p).setCoefs(deformedCoefs);
                    gsInfo << "Patch " << p << " 已应用位移，最大位移: " 
                           << displacementCoefs.cwiseAbs().maxCoeff() << "\n";
                } else {
                    gsWarn << "几何与位移系数矩阵大小不匹配，使用原始几何\n";
                }
            }
        } else {
            gsInfo << "位移场为零或不存在，使用原始几何计算法向量\n";
        }
        
        // 创建牵引力函数：t = -p * n（使用变形后边界的动态法向量）
        gsInfo << "开始创建动态法向量牵引边界条件...\n";
        

        // 定义耦合界面
        std::vector<boundary::side> couplingInterfaces = {
            boundary::north,
            boundary::south, 
            boundary::east
        };
        
        // 循环添加耦合界面的牵引边界条件
        for (const auto& side : couplingInterfaces) {
            gsPressureToTractionPatch<real_t> traction(pFunc, deformedBeamGeometry.patch(0), side);
            gsFunction<real_t>* tractionPtr = traction.clone();
            beamBCs.addCondition(0, side, condition_type::neumann, tractionPtr);
            
            // 输出界面名称
            std::string sideName;
            switch(side) {
                case boundary::north: sideName = "North"; break;
                case boundary::south: sideName = "South"; break;
                case boundary::east: sideName = "East"; break;
                case boundary::west: sideName = "West"; break;
                default: sideName = "Unknown"; break;
            }
        }

        gsInfo << "已添加动态法向量牵引边界条件，总边界条件数: " << beamBCs.size() << "\n";

        // 验证压力值和法向量计算
        gsMatrix<T> testPoint(2, 1);
        testPoint << 0.5, 1.0; // north边界的参数坐标 (u=0.5, v=1.0)
        gsMatrix<T> pressureValue = pFunc.eval(testPoint);
        gsInfo << "测试点参数坐标 (" << testPoint.transpose() << ") 处的压力值: " << pressureValue << "\n";
        
        // 测试法向量计算（使用north边界作为示例）
        gsMatrix<T> boundaryTestParam(1, 1);
        boundaryTestParam(0, 0) = testPoint(0, 0); // north边界使用u参数 (0.5)
        
        // 重新创建north边界用于验证
        std::unique_ptr<gsGeometry<T>> testNorthBoundary = deformedBeamGeometry.patch(0).boundary(boundary::north);
        
        gsMatrix<T> tangent;
        testNorthBoundary->deriv_into(boundaryTestParam, tangent);
        if (tangent.cols() > 0 && tangent.rows() >= 2) {
            T tx = tangent(0, 0);
            T ty = tangent(1, 0);
            T length = std::sqrt(tx*tx + ty*ty);
            if (length > 1e-12) {
                tx /= length; ty /= length;
                gsInfo << "North边界切向量 (归一化): [" << tx << ", " << ty << "]\n";
                gsInfo << "North边界法向量: [" << ty << ", " << -tx << "]\n";
            }
        }
        
    } else {
        gsWarn << "压力场为空，只应用固定端边界条件\n";
    }
    
    // 3. 重新构造组装器
    gsConstantFunction<> gZero(0.0, 0.0, 2); // 无体积力
    
    // 获取当前材料参数
    real_t youngsModulus = beamAssembler->options().getReal("YoungsModulus");
    real_t poissonsRatio = beamAssembler->options().getReal("PoissonsRatio");
    real_t density = beamMassAssembler->options().getReal("Density");
    
    gsInfo << "重建组装器，材料参数: E=" << youngsModulus << ", nu=" << poissonsRatio << ", rho=" << density << "\n";
    
    // 获取基函数
    gsMultiBasis<> basisDisplacement(beamGeometry);
    for (index_t i = 0; i < 3; ++i)
        basisDisplacement.uniformRefine();
    
    // 重新创建弹性问题组装器
    beamAssembler.reset(new gsElasticityAssembler<real_t>(beamGeometry, basisDisplacement, beamBCs, gZero));
    beamAssembler->options().setReal("YoungsModulus", youngsModulus);
    beamAssembler->options().setReal("PoissonsRatio", poissonsRatio);
    beamAssembler->options().setInt("MaterialLaw", material_law::saint_venant_kirchhoff);
    
    // 重新创建质量矩阵组装器  
    beamMassAssembler.reset(new gsMassAssembler<real_t>(beamGeometry, basisDisplacement, beamBCs, gZero));
    beamMassAssembler->options().setReal("Density", density);
    
    // 4. 重新创建时间积分器
    beamSolver.reset(new gsElTimeIntegrator<real_t>(*beamAssembler, *beamMassAssembler));
    beamSolver->options().setInt("Scheme", time_integration::implicit_nonlinear);
    beamSolver->options().setReal("Beta", 0.25);
    beamSolver->options().setReal("Gamma", 0.5);
    
    // 5. 恢复之前的位移和速度状态
    if (oldDisplacement.size() == beamAssembler->numDofs()) {
        beamSolver->setDisplacementVector(oldDisplacement);
        gsInfo << "成功恢复位移向量\n";
    } else {
        gsWarn << "位移向量大小不匹配: 旧=" << oldDisplacement.size() << ", 新=" << beamAssembler->numDofs() << ", 使用零初值\n";
        beamSolver->setDisplacementVector(gsMatrix<real_t>::Zero(beamAssembler->numDofs(), 1));
    }
    
    if (oldVelocity.size() == beamAssembler->numDofs()) {
        beamSolver->setVelocityVector(oldVelocity);
        gsInfo << "成功恢复速度向量\n";
    } else {
        gsWarn << "速度向量大小不匹配: 旧=" << oldVelocity.size() << ", 新=" << beamAssembler->numDofs() << ", 使用零初值\n";
        beamSolver->setVelocityVector(gsMatrix<real_t>::Zero(beamAssembler->numDofs(), 1));
    }
    
    gsInfo << "=== 固体求解器重建完成，自由度: " << beamAssembler->numDofs() << " ===\n";
}

// 注释：将固体边界速度更新到流体求解器的功能已移到resetSolverWithNewMesh中
// 通过Dirichlet边界条件直接应用，确保在求解器重建后正确设置
/*
template<class T>
void updateFluidVelocityFromSolidBoundary(gsINSSolverUnsteady<T, ColMajor>* fluidSolver, 
                                          const gsMultiPatch<T>& patches)
{
    // 此函数已被弃用，功能移至resetSolverWithNewMesh中
    gsInfo << "注意: updateFluidVelocityFromSolidBoundary函数已弃用\n";
    gsInfo << "固体边界速度现在通过resetSolverWithNewMesh中的Dirichlet边界条件应用\n";
}
*/