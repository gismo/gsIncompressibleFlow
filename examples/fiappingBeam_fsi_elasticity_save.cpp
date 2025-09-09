/** @file flappingBeam_fluid.cpp
 
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Li
*/


#include <gismo.h>
#include <cmath>

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
#include <gsElasticity/src/gsPartitionedFSI.h>

using namespace gismo;

template<class T, int MatOrder> void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, int geo);

// reset the solver with a new mesh
template<class T, int MatOrder>
void resetSolverWithNewMesh(gsINSSolverUnsteady<T, MatOrder>*& solver,
                            const gsMultiPatch<T>& newPatches,
                            gsElTimeIntegrator<T>* solidSolver = nullptr);

// transfer the solution to a new mesh
template<class T, int MatOrder>
void transferSolutionToNewMesh(gsINSSolverUnsteady<T, MatOrder>* solver,
                              const gsField<T>& velocityField,
                              const gsField<T>& pressureField,
                              const gsField<T>& solidVelocityField = gsField<T>());

// 根据时间更新几何体边界
template<class T>
gsMultiPatch<T> updateGeometryBoundary(gsMultiPatch<T>& patches, T time, gsINSSolverUnsteady<T, ColMajor>* fluidSolver);

// 先声明initializeSolidSolver函数原型
template<class T>
void initializeSolidSolver(const gsMultiPatch<T>& fluidPatches, T timeStep);

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
    
    short_t domainDim() const override { return m_pressureFunc->domainDim(); }
    short_t targetDim() const override { return 2; }
    
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        try {
            gsMatrix<T> pressure;
            m_pressureFunc->eval_into(u, pressure);
            result.resize(2, u.cols());
            
            // 为每个评估点计算变形后的法向量
            for (index_t i = 0; i < u.cols(); ++i) {
                gsMatrix<T> point = u.col(i);
                
                // 计算该点处变形后边界的法向量
                gsVector<T> normal = computeDeformedNormal(point);
                
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
    gsVector<T> computeDeformedNormal(const gsMatrix<T>& u) const
    {
        gsVector<T> normal(2);
        
        // 边界几何是1D的，需要将2D参数点映射到1D边界参数
        gsMatrix<T> boundaryParam(1, 1);
        
        // 根据边界类型选择合适的参数分量
        switch(m_side) {
            case boundary::north:
            case boundary::south:
                // 对于north/south边界，使用u参数（第一个分量）
                boundaryParam(0, 0) = u(0, 0);
                break;
            case boundary::east:
            case boundary::west:
                // 对于east/west边界，使用v参数（第二个分量）
                boundaryParam(0, 0) = u(1, 0);
                break;
            default:
                // 默认使用第一个参数
                boundaryParam(0, 0) = u(0, 0);
                break;
        }
        
        // 获取边界曲线在参数处的切向量
        gsMatrix<T> tangent;
        m_boundaryGeom->deriv_into(boundaryParam, tangent);
        
        if (tangent.cols() > 0 && tangent.rows() >= 2) {
            T tx = tangent(0, 0);
            T ty = tangent(1, 0);
            
            // 计算切向量的长度
            T length = std::sqrt(tx*tx + ty*ty);
            if (length > 1e-12) {
                // 归一化切向量
                tx /= length;
                ty /= length;
                
                // 法向量 = 切向量逆时针旋转90度
                // 对于外法向量，需要根据边界方向调整
                switch(m_side) {
                    case boundary::north:
                        normal(0) = ty;   // 外法向量指向上方
                        normal(1) = -tx;
                        break;
                    case boundary::south:
                        normal(0) = -ty;  // 外法向量指向下方
                        normal(1) = tx;
                        break;
                    case boundary::east:
                        normal(0) = ty;   // 外法向量指向右方
                        normal(1) = -tx;
                        break;
                    case boundary::west:
                        normal(0) = -ty;  // 外法向量指向左方
                        normal(1) = tx;
                        break;
                    default:
                        // 默认情况，使用逆时针旋转90度
                        normal(0) = -ty;
                        normal(1) = tx;
                        break;
                }
            } else {
                // 如果切向量长度为0，使用默认法向量
                gsWarn << "零长度切向量，使用默认法向量\n";
                getDefaultNormal(normal);
            }
        } else {
            gsWarn << "无法计算切向量，使用默认法向量\n";
            getDefaultNormal(normal);
        }
        
        return normal;
    }
    
    // 获取默认法向量（当计算失败时的后备方案）
    void getDefaultNormal(gsVector<T>& normal) const
    {
        switch(m_side) {
            case boundary::north:
                normal(0) = 0.0; normal(1) = 1.0;
                break;
            case boundary::south:
                normal(0) = 0.0; normal(1) = -1.0;
                break;
            case boundary::east:
                normal(0) = 1.0; normal(1) = 0.0;
                break;
            case boundary::west:
                normal(0) = -1.0; normal(1) = 0.0;
                break;
            default:
                normal(0) = 0.0; normal(1) = 1.0;
                break;
        }
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

// 添加全局变量用于动画文件
static std::ofstream dispPvdFile;
static std::ofstream stressPvdFile;
static bool pvdFilesInitialized = false;

// 移除ALE相关变量
// static std::unique_ptr<gsALE<real_t>> aleModule;
// static gsMultiPatch<real_t> aleDisplacement;
// static gsMultiPatch<real_t> aleVelocity;
// static gsMultiPatch<real_t> aleGeometry;

// 移除FSI模块
// static std::unique_ptr<gsPartitionedFSI<real_t>> fsiModule;

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
    int animStep = 5;

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

    // 最后记得关闭PVD文件
    if (pvdFilesInitialized) 
    {
        // 添加XML文件结束标签
        dispPvdFile << "</Collection>\n</VTKFile>\n";
        stressPvdFile << "</Collection>\n</VTKFile>\n";
        
        // 关闭文件
        dispPvdFile.close();
        stressPvdFile.close();
        
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
                updateGeometryBoundary(currentPatches, currentTime, pSolver);
                currentPatches = updateGeometryBoundary(currentPatches, currentTime, pSolver);
                
                // 更新网格和求解器
                resetSolverWithNewMesh(pSolver, currentPatches, beamSolver);
        

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
                           gsElTimeIntegrator<T>* solidSolver)
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

    // 3. Obtain current solid velocity
    gsField<T> solidVelocityField = solidSolver->constructSolution(0);

    gsDebug << "Current solid velocity field:\n" << solidVelocityField << "\n";

    
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

    // 验证每个patch的有效性
    for (size_t i = 0; i < persistentPatches->nPatches(); ++i) {
        const auto& patch = persistentPatches->patch(i);
        if (patch.coefs().rows() == 0) {
            gsWarn << "CRITICAL: Patch " << i << " 在持久化后变为空!\n";
        }
    }

    // 6. Project old solution onto new mesh with solid velocity boundary conditions
    transferSolutionToNewMesh(newSolver, velocityFieldOld, pressureFieldOld, solidVelocityField);
    
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
    
    gsInfo << "=== resetSolverWithNewMesh 完成 ===\n";
    
    // Note: We don't delete oldSolver in case it's externally managed
    // We also don't delete persistentPatches and persistentPde as they're managed by the solver
}

template<class T, int MatOrder>
void transferSolutionToNewMesh(gsINSSolverUnsteady<T, MatOrder>* solver,
                              const gsField<T>& velocityField,
                              const gsField<T>& pressureField,
                              const gsField<T>& solidVelocityField = gsField<T>())
{
    // 获取新的基函数和补丁
    const std::vector<gsMultiBasis<T>>& newBases = solver->getAssembler()->getBases();
    const gsMultiPatch<T>& newPatches = solver->getParams()->getPde().patches();

    // 添加调试信息，检查一致性
    gsDebugVar(newPatches.nPatches());
    gsDebugVar(newBases.size());
    for (size_t i = 0; i < newBases.size(); ++i) {
        gsDebugVar(newBases[i].nBases());
    }
    
    // 获取目标维度和自由度数量
    const index_t vDim = velocityField.function(0).targetDim();
    const index_t fullUdofs = solver->getAssembler()->getUdofs();
    const index_t tarDim = solver->getParams()->getPde().domain().geoDim();
    const index_t pShift = tarDim * fullUdofs;
    
    // 创建与setSolutionCoefs期望的尺寸匹配的矩阵
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
        if (p >= newBases[0].nBases()) {
            gsInfo << "Invalid patch index " << p << " for velocity basis. Skipping." << "\n";
            continue;
        }
        
        if (p >= velocityField.nPatches()) {
            gsInfo << "Skipping velocity interpolation for patch " << p << " as it does not exist in the old field.\n";
            continue;
        }

        // 获取当前patch的基函数
        const gsBasis<T>& vBasis = newBases[0].basis(p);
        
        for (size_t i = 0; i < vBasis.size(); ++i)
        {
            try {
                // 对每个基函数进行插值
                gsMatrix<T> coef = quasiInterp.localIntpl(vBasis, velocityField.function(p), i);
                
                // 将系数存储到全局系数矩阵中 - 使用映射器获取全局索引
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
        if (p >= newBases[1].nBases()) 
        {
            gsInfo << "Invalid patch index " << p << " for pressure basis. Skipping." << "\n";
            continue;
        }
        
        if (p >= pressureField.nPatches()) {
            gsInfo << "Skipping pressure interpolation for patch " << p << " as it does not exist in the old field.\n";
            continue;
        }
        
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

    // 如果提供了固体速度场，更新流体边界速度
    if (solidVelocityField.nPatches() > 0)
    {
        gsInfo << "Updating fluid boundary velocity from solid velocity field\n";
        
        // 获取流体-固体界面的边界条件信息
        const gsBoundaryConditions<>& bc = solver->getParams()->getPde().bc();
        
        // 遍历所有边界条件，找到流体-固体界面
        for (index_t p = 0; p < newPatches.nPatches(); ++p)
        {
            // 获取当前patch的边界条件
            for (boxSide side : newPatches.patch(p).boundary())
            {
                patchSide ps(p, side);
                
                // 检查是否是流体-固体界面（通常是Dirichlet边界条件）
                if (bc.getConditionFromSide(ps, 0) == condition_type::dirichlet)
                {
                    // 获取该边界上的基函数
                    const gsBasis<T>& vBasis = newBases[0].basis(p);
                    gsMatrix<T> bnd = vBasis.boundary(side);
                    
                    // 获取边界上的自由度索引
                    gsMatrix<index_t> bndDofs = vBasis.boundaryOffset(side, 0);
                    
                    // 在边界上评估固体速度
                    gsMatrix<T> evalPoints;
                    vBasis.eval_into(bnd, evalPoints);
                    
                    // 将物理坐标转换为固体域的参数坐标
                    gsMatrix<T> solidParams;
                    if (solidVelocityField.nPatches() > 0)
                    {
                        // 假设固体只有一个patch
                        solidVelocityField.patch(0).invertPoints(evalPoints, solidParams);
                        
                        // 在这些点上评估固体速度
                        gsMatrix<T> solidVel = solidVelocityField.function(0).eval(solidParams);
                        
                        // 更新流体速度系数
                        for (index_t i = 0; i < bndDofs.size(); ++i)
                        {
                            if (mappers[0].is_free(bndDofs(i), p))
                            {
                                index_t globalIndex = mappers[0].index(bndDofs(i), p);
                                
                                for (index_t d = 0; d < vDim; ++d)
                                {
                                    newVelocityCoefs(globalIndex + d * fullUdofs, 0) = solidVel(d, i);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // 更新速度系数
        solver->setSolutionCoefs(newVelocityCoefs, 0);
    }

    gsDebugVar(newVelocityCoefs.size());
    gsDebugVar(newPressureCoefs.size());
    // 删除或注释掉以下行，避免在此时构建解
    // gsField<T> velocityResult = solver->constructSolution(0);
    // gsField<T> pressureResult = solver->constructSolution(1);
    // gsWriteParaview<>(velocityResult, "velocityResult", 1000, false);
    // gsWriteParaview<>(pressureResult, "pressureResult", 1000, false);

    gsInfo << "Projected solution to new mesh\n";
}


template<class T>
gsMultiPatch<T> updateGeometryBoundary(gsMultiPatch<T>& patches, T time, gsINSSolverUnsteady<T, ColMajor>* fluidSolver)
{
    gsInfo << "\n===== ENTERING updateGeometryBoundary at time = " << time << " =====\n";
    gsInfo << "Function is being called successfully!\n";
    
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
        
        // 初始化PVD文件用于动画
        if (!pvdFilesInitialized) {
            // 创建位移动画文件
            dispPvdFile.open("beam_displacement_animation.pvd");
            GISMO_ASSERT(dispPvdFile.is_open(), "无法创建位移动画文件");
            dispPvdFile << "<?xml version=\"1.0\"?>\n"
                      << "<VTKFile type=\"Collection\" version=\"0.1\">\n"
                      << "<Collection>\n";
            
            // 创建应力动画文件
            stressPvdFile.open("beam_stress_animation.pvd");
            GISMO_ASSERT(stressPvdFile.is_open(), "无法创建应力动画文件");
            stressPvdFile << "<?xml version=\"1.0\"?>\n"
                        << "<VTKFile type=\"Collection\" version=\"0.1\">\n"
                        << "<Collection>\n";
            
            pvdFilesInitialized = true;
            gsInfo << "初始化了位移和应力动画文件\n";
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
                
                std::unique_ptr<gsGeometry<T>> testBoundary(testDeformedGeometry.patch(0).boundary(boundary::north));
                gsPressureToTraction<real_t> testTraction(pressureField.function().function(0), *testBoundary, boundary::north);
                gsMatrix<T> tractionValue;
                testTraction.eval_into(testPoint, tractionValue);
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
            std::unique_ptr<gsGeometry<T>> originalBoundary(beamGeometry.patch(0).boundary(boundary::north));
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
                
                std::unique_ptr<gsGeometry<T>> deformedBoundary(currentDeformed.patch(0).boundary(boundary::north));
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
        
        if (pvdFilesInitialized) {
            // 添加到PVD文件（动画集合）
            dispPvdFile << "<DataSet timestep=\"" << time << "\" part=\"0\" file=\"" 
                       << dispFileName << "0.vts\"/>\n";
            stressPvdFile << "<DataSet timestep=\"" << time << "\" part=\"0\" file=\"" 
                        << stressFileName << "0.vts\"/>\n";
            
            // 确保写入到文件（刷新缓冲区）
            dispPvdFile.flush();
            stressPvdFile.flush();
        }
        
        gsInfo << "固体时间步完成，迭代次数: " << beamSolver->numberIterations() << "\n";
        gsInfo << "已保存位移场到: " << dispFileName << "\n";
        gsInfo << "已保存应力场到: " << stressFileName << "\n";
    }
    
    // ==== 直接更新流体网格 ====
    // Reset 'patches' to the initial fluid domain configuration.
    patches = initial_fluid_patches_static; 
    // 创建位移场并获取函数对象
    gsField<T> displacementField(beamGeometry, beamDisplacement);
    const gsFunctionSet<T>& displacementFunc = displacementField.function();

    // 1. 找出所有需要更新的patches (与梁接触的流体区域)
    // 假设梁在patch 3, 4, 5的区域附近
    gsInfo << "Updating fluid mesh boundary in contact with the beam...\n";
    
    // 定义与梁交界的边界对应关系 - 修改为包含更详细的信息
    struct BeamInterface {
        size_t patchIdx;          // patches中的索引
        boundary::side boundSide; // 对应的边界方向
        // bool isCorner;            // 是否是转角处(可能需要特殊处理)
    };
    
    // 根据几何结构确定所有与梁接触的边界
    std::vector<BeamInterface> interfaces = 
    {
        {3, boundary::south},
        {4, boundary::north},
        {5, boundary::west},
    };
    
    // 先更新边界控制点
    for (const auto& interface : interfaces)
    {
        gsInfo << "Processing patch " << interface.patchIdx << " of " 
               << interface.boundSide << " boundary\n";
        
        // 获取patch边界
        // 'patches' is now a copy of the initial fluid mesh (due to 'patches = initial_fluid_patches_static;')
        // So, these are the initial coordinates of the fluid boundary segment.
        gsMatrix<T> initial_coefs_NxD = patches.patch(interface.patchIdx).boundary(interface.boundSide)->coefs();

        // Transpose to (Dim x N_points) for evaluation, as gsFunction::eval typically expects points as columns.
        gsMatrix<T> points_for_eval_DxN = initial_coefs_NxD.transpose();
        
        // Define vDim for dimension consistency
        const index_t vDim = points_for_eval_DxN.rows(); // Should be 2 for 2D problem
        
        // 保存原始边界控制点用于对比
        gsInfo << "Original boundary control points for patch " << interface.patchIdx << ":\n";
        for (index_t k = 0; k < std::min(index_t(3), initial_coefs_NxD.rows()); ++k) {
            gsInfo << "  Point " << k << ": (" << initial_coefs_NxD(k, 0) << ", " << initial_coefs_NxD(k, 1) << ")\n";
        }
        
        // Evaluate displacement at these initial material points.
        // displacementFunc is derived from beamDisplacement (total displacement) on beamGeometry (initial solid config).
        gsMatrix<T> displacement_vectors_DxN = gsMatrix<T>::Zero(vDim, points_for_eval_DxN.cols());

        gsVector<T> phys(2), par(2);
        for (index_t j=0; j<points_for_eval_DxN.cols(); ++j)
        {
            phys = points_for_eval_DxN.col(j);
            // Use invertPoints instead of invMap
            gsMatrix<T> parMatrix(2, 1);
            beamGeometry.patch(0).invertPoints(phys, parMatrix, 1e-8);
            par = parMatrix.col(0);
            
            gsMatrix<T> vv = displacementFunc.eval(parMatrix);
            if (vv.rows() >= vDim && vv.cols() > 0) {
                for (index_t d = 0; d < vDim; ++d) {
                    displacement_vectors_DxN(d, j) = vv(d, 0);
                }
            }
        }

        // Log displacement magnitude to verify
        if (displacement_vectors_DxN.size() > 0) 
        { // Check if not empty before norm or access
            //  gsInfo << "Displacement norm for fluid patch " << interface.patchIdx << " side " << interface.boundSide 
            //         << ": " << displacement_vectors_DxN.norm() << " (原始位移已放大5倍)\n";
             if (displacement_vectors_DxN.rows() > 0 && displacement_vectors_DxN.cols() > 0) {
                gsInfo << "  Sample displacement value [0,0]: " << displacement_vectors_DxN(0,0) << "\n";
                gsInfo << "  Max displacement: " << displacement_vectors_DxN.cwiseAbs().maxCoeff() << "\n";
             }
        }
        
        // 应用位移到边界控制点
        gsMatrix<T> new_coefs_NxD = initial_coefs_NxD;
        
        // 添加位移到控制点
        for (index_t j = 0; j < displacement_vectors_DxN.cols(); ++j) 
        {
            for (index_t d = 0; d < vDim; ++d) 
            {
                new_coefs_NxD(j, d) += displacement_vectors_DxN(d, j);
            }
        }
        
        gsInfo << "Updated boundary control points for patch " << interface.patchIdx << ":\n";
        for (index_t k = 0; k < std::min(index_t(3), new_coefs_NxD.rows()); ++k) {
            gsInfo << "  Point " << k << ": (" << new_coefs_NxD(k, 0) << ", " << new_coefs_NxD(k, 1) << ")";
            gsInfo << " [Δx=" << (new_coefs_NxD(k, 0) - initial_coefs_NxD(k, 0)) << ", Δy=" << (new_coefs_NxD(k, 1) - initial_coefs_NxD(k, 1)) << "]\n";
        }

        // 直接修改整个patch的控制点矩阵，而不是通过边界函数
        gsGeometry<T>& currentPatch = patches.patch(interface.patchIdx);
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
        
        // 立即验证边界控制点是否被正确修改
        gsMatrix<T> verification_coefs = currentPatch.boundary(interface.boundSide)->coefs();
        gsInfo << "立即验证patch " << interface.patchIdx << " 边界控制点修改:\n";
        for (index_t k = 0; k < std::min(index_t(2), verification_coefs.rows()); ++k) {
            gsInfo << "  验证点 " << k << ": (" << verification_coefs(k, 0) << ", " << verification_coefs(k, 1) << ")\n";
        }
        
        gsInfo << "对patch " << interface.patchIdx << " 应用了真实的固体边界变形\n";
    }

    // 先保存变形网格用于调试（在应用完所有位移后）
    gsWriteParaview<>(patches, "deformed_mesh_time_" + util::to_string(time), 1000);
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
        gsMatrix<T> deformed_coefs = patches.patch(interface.patchIdx).boundary(interface.boundSide)->coefs();
        gsInfo << "最终验证：Patch " << interface.patchIdx << " 变形后边界点坐标:\n";
        for (index_t k = 0; k < std::min(index_t(2), deformed_coefs.rows()); ++k) {
            gsInfo << "  点 " << k << ": (" << deformed_coefs(k, 0) << ", " << deformed_coefs(k, 1) << ")\n";
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
        gsMatrix<T> coefs1 = geom1.boundary(side1)->coefs();
        gsMatrix<T> coefs2 = geom2.boundary(side2)->coefs();
        
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
    
    gsInfo << "=== 重建固体求解器以应用新的压力场 ===\n";
    
    // 保存当前的位移和速度状态
    // 注意：位移数据从solutionVector中提取，因为displacementVector()方法被注释掉了
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
        
        try {
            // north边界
            std::unique_ptr<gsGeometry<T>> northBoundary(deformedBeamGeometry.patch(0).boundary(boundary::north));
            gsPressureToTraction<real_t> tractionNorth(pFunc, *northBoundary, boundary::north);
            gsFunction<real_t>* tractionNorthPtr = tractionNorth.clone();
            beamBCs.addCondition(0, boundary::north, condition_type::neumann, tractionNorthPtr);
            gsInfo << "✓ North边界牵引条件已添加\n";
            
            // south边界
            std::unique_ptr<gsGeometry<T>> southBoundary(deformedBeamGeometry.patch(0).boundary(boundary::south));
            gsPressureToTraction<real_t> tractionSouth(pFunc, *southBoundary, boundary::south);
            gsFunction<real_t>* tractionSouthPtr = tractionSouth.clone();
            beamBCs.addCondition(0, boundary::south, condition_type::neumann, tractionSouthPtr);
            gsInfo << "✓ South边界牵引条件已添加\n";
            
            // east边界
            std::unique_ptr<gsGeometry<T>> eastBoundary(deformedBeamGeometry.patch(0).boundary(boundary::east));
            gsPressureToTraction<real_t> tractionEast(pFunc, *eastBoundary, boundary::east);
            gsFunction<real_t>* tractionEastPtr = tractionEast.clone();
            beamBCs.addCondition(0, boundary::east, condition_type::neumann, tractionEastPtr);
            gsInfo << "✓ East边界牵引条件已添加\n";
        } catch (const std::exception& e) {
            gsWarn << "创建牵引边界条件时出错: " << e.what() << "\n";
            throw;
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
        std::unique_ptr<gsGeometry<T>> testNorthBoundary(deformedBeamGeometry.patch(0).boundary(boundary::north));
        
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