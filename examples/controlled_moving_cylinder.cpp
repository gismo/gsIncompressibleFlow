/** @file controlled_moving_cylinder.cpp
    @brief ALE算例：经典受控振荡圆柱（圆柱以给定规律在y向振荡）

    几何采用基准通道+圆柱（flow_around_cylinder.xml），外边界固定，入口抛物线速度，出口自然边界，
    圆柱壁面通过ALE网格位移施加随动边界（u = w）。

    本文件属于G+Smo库。
    MPL-2.0 许可；详情见 http://mozilla.org/MPL/2.0/。
*/

#include <gismo.h>
#include <fstream>
#include <cstdio>
#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    // ---------------------- 问题参数 ----------------------
    real_t meanVelocity = 1.0;     // 入口平均速度（基准2D-1常用）
    real_t rho = 1.0;              // 流体密度
    real_t Re = 100.0;             // 雷诺数（由几何L=0.1、U~0.2可估）
    real_t Lc = 0.1;               // 特征长度（圆柱直径）clear
    real_t nu = meanVelocity * Lc / Re; // 动力粘度（或运动粘度，取决于无量纲化）

    // 时间参数
    real_t T = 4.0;                // 总时间（增加到4秒，看完整的4个周期）
    real_t dt = 0.01;              // 时间步长（增大到0.01，减少计算量）
    index_t nTimeSteps = static_cast<index_t>(std::ceil(T / dt));
    index_t outputInterval = 2;    // 每N步输出一次（减少输出频率）

    // 圆柱受控振荡参数（y向） y_c(t) = A * sin(2π f t)
    real_t cylA = 0.4;            // 震幅（增大到半径0.05的160%，更明显）
    real_t cylF = 10.0;             // 频率（Hz，增加到1Hz使运动更快）
    gsVector<real_t> cylCenter0(2); cylCenter0 << 0.2, 0.2; // 初始中心
    real_t cylR = 0.05;            // 半径

    // 命令行
    gsCmdLine cmd("Controlled moving cylinder with ALE.");
    cmd.addReal("u","umean","Mean inflow velocity",meanVelocity);
    cmd.addReal("v","nu","Viscosity (nu)",nu);
    cmd.addReal("t","time","Total time",T);
    cmd.addReal("s","step","Time step",dt);
    cmd.addInt ("o","output","Output interval",outputInterval);
    cmd.addReal("A","ampl","Cylinder oscillation amplitude",cylA);
    cmd.addReal("f","freq","Cylinder oscillation frequency",cylF);
    std::string geomPath;
    bool enableParam = false;    // 是否启用参数化（默认关闭）
    int  paramEveryN = 1;        // 参数化频率：每N步执行一次（1=每步，>1=降频，0=关闭）
    int  numRefine   = 0;        // 均匀h细化次数（对速度与压力几何同步细化）
    cmd.addString("g","geom","Geometry XML path (flow_around_cylinder.xml)",geomPath);
    cmd.addSwitch("P","param","Enable mesh re-parameterization", enableParam);
    cmd.addInt   ("F","paramfreq","Mesh re-parameterization frequency (every N steps, 0=off)", paramEveryN);
    cmd.addInt   ("r","refine","Number of uniform h-refinements (velocity & pressure)", numRefine);
    bool stokesInit = false;     // 是否用Stokes做初值（默认关闭，直接用NS时间步推进）
    cmd.addSwitch("S","stokesInit","Use Stokes solve to initialize at t=0", stokesInit);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // 依据命令行参数更新步数
    nTimeSteps = static_cast<index_t>(std::ceil(T / dt));
    
    gsInfo << "=== Controlled Moving Cylinder (ALE) ===\n";
    gsInfo << "nu: " << nu << ", mean U: " << meanVelocity << ", dt: " << dt
           << ", T: " << T << ", steps: " << nTimeSteps << "\n";

    // ---------------------- 读取几何（通道+圆柱） ----------------------
    // 使用现有基准几何：2.2 x 0.41 通道，圆柱 r=0.05, center=(0.2,0.2)
    gsMultiPatch<> fluidDomain;
    {
        auto fileExists = [](const std::string &p)->bool{
            std::ifstream ifs(p.c_str());
            return static_cast<bool>(ifs);
        };

        std::string geomFile;
        if (!geomPath.empty() && fileExists(geomPath))
        {
            geomFile = geomPath;
        }
        else
        {
            // 从常见位置尝试查找（相对构建目录/源码根）
            const char* candidates[] = {
                "../optional/gsElasticity/filedata/flow_around_cylinder.xml",
                "../../optional/gsElasticity/filedata/flow_around_cylinder.xml",
                "../../../optional/gsElasticity/filedata/flow_around_cylinder.xml",
                "optional/gsElasticity/filedata/flow_around_cylinder.xml",
                "flow_around_cylinder.xml"
            };
            for (const char* c : candidates)
            {
                if (fileExists(c)) { geomFile = c; break; }
            }
        }

        if (geomFile.empty())
        {
            gsWarn << "Cannot locate flow_around_cylinder.xml. Please pass with -g <path>.\n";
            return -1;
        }
        gsInfo << "Using geometry: " << geomFile << "\n";
        gsReadFile<>(geomFile, fluidDomain);
    }

    // ---------------------- 构建离散空间（Taylor–Hood） ----------------------
    // 为保证几何控制点与速度DOF索引一致性：
    // - 压力基于原始几何（P1）
    // - 将几何升阶后作为速度几何（P2）以匹配速度DOF
    gsMultiPatch<> pressureDomain = fluidDomain;    // 基于原始几何的压力域
    pressureDomain.degreeElevate(2);                // 压力几何升一阶：Pk -> Pk+1
    fluidDomain.degreeElevate(3);                   // 速度几何升两阶：Pk -> Pk+2（仍保持速度比压力高一阶）

    // 升阶后需重建拓扑，确保接口正确黏合
    fluidDomain.computeTopology();
    pressureDomain.computeTopology();

    // ---------------------- 均匀h细化（可选）----------------------
    if (numRefine > 0)
    {
        gsInfo << "Applying uniform h-refinement times = " << numRefine << "...\n";
        for (int rr = 0; rr < numRefine; ++rr)
        {
            fluidDomain.uniformRefine();
            pressureDomain.uniformRefine();
        }
        fluidDomain.computeTopology();
        pressureDomain.computeTopology();
    }

    gsMultiBasis<> basisVelocity(fluidDomain);      // 速度P2
    gsMultiBasis<> basisPressure(pressureDomain);   // 压力P1

    // 保存原始网格，用于在ALE位移函数中基于原始控制点判定边界
    gsMultiPatch<> originalFluidDomain = fluidDomain;

    // ---------------------- 自动识别圆柱边界所在的参数边 ----------------------
    std::vector<int> circleSide(fluidDomain.nPatches(), -1);
    auto detectCircleSide = [&](index_t p)->int
    {
        const gsGeometry<>& geo = fluidDomain.patch(p);
        const int sides[4] = {boundary::west, boundary::east, boundary::south, boundary::north};
        real_t bestScore = std::numeric_limits<real_t>::infinity();
        int bestSide = -1;
        for (int sIdx = 0; sIdx < 4; ++sIdx)
        {
            int side = sides[sIdx];
            // 在该边上均匀采样
            const int ns = 16;
            gsMatrix<> uv(2, ns);
            for (int k = 0; k < ns; ++k)
            {
                real_t t = static_cast<real_t>(k) / static_cast<real_t>(ns - 1);
                if (side == boundary::west)  { uv(0, k) = 0.0; uv(1, k) = t; }
                if (side == boundary::east)  { uv(0, k) = 1.0; uv(1, k) = t; }
                if (side == boundary::south) { uv(0, k) = t;  uv(1, k) = 0.0; }
                if (side == boundary::north) { uv(0, k) = t;  uv(1, k) = 1.0; }
            }
            gsMatrix<> XY;
            geo.eval_into(uv, XY); // 2 x ns
            real_t avgDiff = 0.0;
            for (int k = 0; k < ns; ++k)
            {
                real_t dx = XY(0, k) - cylCenter0[0];
                real_t dy = XY(1, k) - cylCenter0[1];
                real_t r  = std::sqrt(dx * dx + dy * dy);
                avgDiff += std::abs(r - cylR);
            }
            avgDiff /= static_cast<real_t>(ns);
            if (avgDiff < bestScore)
            {
                bestScore = avgDiff;
                bestSide = side;
            }
        }
        // 阈值：距离平均偏差足够小，认为该边是圆柱边界
        if (bestScore < 1e-3)
            return bestSide;
        return -1;
    };
    for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
        circleSide[p] = detectCircleSide(p);

    // 若需要边界层加密，可在此处基于 circleSide 对相应参数方向插入分级结点；
    // 当前示例先提供统一h细化（-r N）。

    // ---------------------- 边界条件 ----------------------
    gsBoundaryConditions<> bcInfo;
    
    // 入口抛物线速度 U(y) = U_mean * 6 y (H-y)/H^2，H=0.41
    std::string uInExpr = util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2";
    gsFunctionExpr<> inletVec(uInExpr, "0", 2);   // 向量函数 [u_in(y), 0]
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    // 圆柱边界(流固界面)的法向速度：u_wall = [0, 2*pi*f*A*cos(2*pi*f*t)]
    // 使用表达式函数，t 由时间循环中动态设置
    const real_t wallVyAmp = 2.0 * M_PI * cylF * cylA;
    std::string wallVyExpr = util::to_string(wallVyAmp) + "*cos(2*pi*" + util::to_string(cylF) + "*t)";
    gsFunctionExpr<> cylWallVelX("0", 2);        // x 分量恒为 0 (依赖 x,y )
    gsFunctionExpr<> cylWallVelY(wallVyExpr, 2);  // y 分量为随时间变化的速度 (依赖 x,y )
    gsConstantFunction<> zeroP(0.0, 2);

    // 基于包围盒自动识别外边界：左=入口，上下=墙，右=自然
    real_t xmin =  std::numeric_limits<real_t>::infinity();
    real_t xmax = -std::numeric_limits<real_t>::infinity();
    real_t ymin =  std::numeric_limits<real_t>::infinity();
    real_t ymax = -std::numeric_limits<real_t>::infinity();
    for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
    {
        const gsMatrix<>& C = fluidDomain.patch(p).coefs();
        for (index_t i = 0; i < C.rows(); ++i)
        {
            xmin = std::min(xmin, C(i,0)); xmax = std::max(xmax, C(i,0));
            ymin = std::min(ymin, C(i,1)); ymax = std::max(ymax, C(i,1));
        }
    }
    const int sides[4] = {boundary::west, boundary::east, boundary::south, boundary::north};
    for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
    {
        const gsGeometry<>& geo = fluidDomain.patch(p);
        for (int sIdx = 0; sIdx < 4; ++sIdx)
        {
            int side = sides[sIdx];
            const int ns = 12;
            gsMatrix<> uv(2, ns);
            for (int k = 0; k < ns; ++k)
            {
                real_t t = static_cast<real_t>(k) / static_cast<real_t>(ns - 1);
                if (side == boundary::west)  { uv(0, k) = 0.0; uv(1, k) = t; }
                if (side == boundary::east)  { uv(0, k) = 1.0; uv(1, k) = t; }
                if (side == boundary::south) { uv(0, k) = t;  uv(1, k) = 0.0; }
                if (side == boundary::north) { uv(0, k) = t;  uv(1, k) = 1.0; }
            }
            gsMatrix<> XY; geo.eval_into(uv, XY);

            auto avgAbs = [&](int dim, real_t val){
                real_t a = 0.0; for (int k = 0; k < ns; ++k) a += std::abs(XY(dim,k) - val); return a/static_cast<real_t>(ns);
            };
            bool isLeft   = (avgAbs(0, xmin) < 1e-8);
            bool isRight  = (avgAbs(0, xmax) < 1e-8);
            bool isBottom = (avgAbs(1, ymin) < 1e-8);
            bool isTop    = (avgAbs(1, ymax) < 1e-8);

            if (isLeft)
            {
                bcInfo.addCondition(p, side, condition_type::dirichlet, &inletVec, 0);
            }
            else if (isBottom || isTop)
            {
                for (index_t d = 0; d < 2; ++d)
                    bcInfo.addCondition(p, side, condition_type::dirichlet, &zeroVel, 0, d);
            }
            // 右侧自然边界：不加条件
        }
    }

    // 注意：圆柱边界不显式施加速度Dirichlet，而是通过ALE移动边界并在对流项中使用相对速度(u-w)
    // 来实现无滑移的移动壁面条件。若在此处对圆柱加速度Dirichlet，会把边界自由度固定，
    // 使得网格位移函数无法驱动边界控制点，从而“圆柱不移动”。

    // ---------------------- 定义PDE与求解器参数 ----------------------
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &zeroVel, nu);
    
    std::vector<gsMultiBasis<>> discreteBases;
    discreteBases.push_back(basisVelocity);  // 速度P2
    discreteBases.push_back(basisPressure);  // 压力P1
    
    gsFlowSolverParams<real_t> params(nsPde, discreteBases);
    params.options().setInt ("nonlin.maxIt", 100);
    params.options().setReal("nonlin.tol", 1e-6);
    params.options().setInt ("lin.maxIt", 200);
    params.options().setReal("lin.tol", 1e-8);
    params.options().setReal("timeStep", dt);
    params.options().setSwitch("quiet", false);

    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));
    solver.initialize();
    solver.setALEActive(true);
    
    // 网格优化配置（参数化）：只对圆柱周围的4个patch做重参数化，跳过最右侧patch
    solver.setMeshOptimization(enableParam);
    solver.setMeshOptimizationCrossPatch(true);
    solver.getMeshOptOptions().setInt("Verbose", 0);           // 减少调试输出
    solver.getMeshOptOptions().setInt("ParamMethod", 1);        // VariationalHarmonicPatch 方法
    solver.getMeshOptOptions().setInt("AAPreconditionType", 0); // 标准预处理

    // 仅选择 circleSide[p] != -1 的patch进行参数化（即圆柱四周的四块），避免右侧大域被参数化
    {
        std::vector<index_t> selected;
        for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
            if (circleSide[p] != -1)
                selected.push_back(p);
        if (!selected.empty())
            solver.setMeshOptimizationPatches(selected);
        gsInfo << "[Param] Selected patches for reparameterization: ";
        for (auto pid : selected) gsInfo << pid << " ";
        gsInfo << "\n";
    }

    // 若启用参数化，则在求解初始 Stokes 前先做一次参数化，
    // 并同步 ALE 历史，避免第0步与第1步之间的跳变
    if (enableParam)
    {
        gsInfo << "[Param] Applying initial mesh optimization at t=0...\n";
        solver.applyInitialMeshOptimization();
    }
    
    // 方案3：使用弹性方程方法（如rotation_square_proper_ale.cpp中的方法）
    // solver.setMeshOptimization(true);
    // solver.getMeshOptOptions().setInt("Verbose", 0);
    // solver.getMeshOptOptions().setInt("ParamMethod", 2);  // ElasticPatch方法
    // solver.getMeshOptOptions().setInt("AAPreconditionType", 0);

    // ---------------------- 网格运动：仅驱动圆柱边界；外边界默认不动 ----------------------
    // 注意：ALE求解器在第n步时传递的时间是(n-1)*dt，需要修正
    auto meshMotion = [&](real_t t) -> gsMatrix<real_t>
    {
        // 使用求解器传递的当前时间 t（与内部时间循环同步）
        real_t current_t = t;
        
        // 圆柱中心当前 y 位移（刚体平移）
        real_t dy_disp = cylA * std::sin(2.0 * M_PI * cylF * current_t);
        
        // 调试输出：显示当前时间和位移
        static int callCount = 0;
        if (callCount % 10 == 0) // 每10次调用输出一次，避免输出过多
        {
            gsInfo << "[ALE Motion] t=" << current_t 
                   << ", dy_disp=" << dy_disp << " (A=" << cylA << ", f=" << cylF << ")\n";
        }
        callCount++;

        const index_t udofs = solver.getAssembler()->getUdofs();
        gsMatrix<> disp(2 * udofs, 1); 
        disp.setZero();

        const gsDofMapper& uMapper = solver.getAssembler()->getMappers()[0];

        // 仅对圆柱边界施加位移（其余默认0，由网格优化平滑扩展）。
        index_t movedCount = 0;
        index_t totalChecked = 0;
        index_t freeCount = 0;
        
        // 准确选择圆柱参数边界上的控制点：使用已识别的 circleSide[p]
        for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
        {
            if (circleSide[p] == -1)
                continue; // 非圆柱边界的patch不直接施加边界位移

            // 获取该patch在圆柱对应参数边上的控制点索引
            gsMatrix<index_t> bIdx = basisVelocity.basis(p).boundary(circleSide[p]);
            const gsMatrix<>& origCoefs = originalFluidDomain.patch(p).coefs();
            std::vector<char> moved(origCoefs.rows(), 0);

            for (index_t k = 0; k < bIdx.rows(); ++k)
            {
                index_t i = bIdx(k,0);
                totalChecked++;

                if (!uMapper.is_free(i, p))
                    continue;

                freeCount++;
                const real_t x = origCoefs(i,0);
                const real_t y = origCoefs(i,1);

                index_t base = uMapper.index(i, p);
                if (base < udofs && base + udofs < disp.rows())
                {
                    disp(base) = 0.0;              // 刚体平移：x不动
                    disp(base + udofs) = dy_disp;  // y整体平移
                    ++movedCount;
                    moved[i] = 1;

                    if (callCount == 0 && movedCount <= 10)
                    {
                        gsInfo << "[Debug] Applied rigid translation on circle CP (" << x << ", " << y
                               << ") patch=" << p << ", id=" << i << ", base=" << base
                               << ", dy_disp=" << dy_disp << "\n";
                    }
                }
            }

            // 兜底：若某些圆柱边界角点未被 boundary(side) 捕获，则用几何半径判据补齐（仅该patch）
            const real_t onCircleTol = 5e-6; // 严格一些，避免误选
            for (index_t i = 0; i < origCoefs.rows(); ++i)
            {
                if (moved[i]) continue;
                if (!uMapper.is_free(i, p)) continue;

                const real_t x = origCoefs(i,0);
                const real_t y = origCoefs(i,1);
                const real_t dx = x - cylCenter0[0];
                const real_t dy = y - cylCenter0[1];
                const real_t r  = std::sqrt(dx*dx + dy*dy);
                if (std::abs(r - cylR) < onCircleTol)
                {
                    index_t base = uMapper.index(i, p);
                    if (base < udofs && base + udofs < disp.rows())
                    {
                        disp(base) = 0.0;
                        disp(base + udofs) = dy_disp;
                        ++movedCount;
                        moved[i] = 1;
                    }
                }
            }
        }
        // 可选：首步打印一次，便于确认确实检出圆柱边界控制点
        static bool firstCall = true;
        if (firstCall)
        {
            gsInfo << "[ALE] Total points checked: " << totalChecked << "\n";
            gsInfo << "[ALE] Free DOF points: " << freeCount << "\n";
            gsInfo << "[ALE] cylinder boundary DOFs moved = " << movedCount << "\n";
            gsInfo << "[ALE] Cylinder center: (" << cylCenter0[0] << ", " << cylCenter0[1] 
                   << "), radius=" << cylR << ", tolerance=" << 5e-3 << "\n";
            gsInfo << "[ALE] Total velocity DOFs: " << udofs << "\n";
            gsInfo << "[ALE] Displacement vector size: " << disp.rows() << "\n";
            firstCall = false;
        }
        return disp;
    };
    solver.setMeshUpdateFunction(meshMotion);
    
    // 验证网格更新函数是否设置成功
    gsInfo << "[ALE Status] Mesh update function set successfully\n";
    
    // ---------------------- 初始Stokes、时间步进与输出 ----------------------
    if (stokesInit)
    {
        gsInfo << "Solving initial Stokes...\n";
        // 若使用显式壁面速度Dirichlet，这里需要设置时间变量。
        // 当前采用ALE相对速度方式，无需设置。
        solver.solveStokes();
    }
    else
    {
        gsInfo << "Skip Stokes init: start directly with Navier–Stokes time stepping.\n";
    }
    
    unsigned samp = 800; // 采样密度
    
    // 创建ParaView动画集合
    gsParaviewCollection collectionVelocity("cmc_velocity_animation");
    gsParaviewCollection collectionPressure("cmc_pressure_animation");
    gsParaviewCollection collectionMeshDisp("cmc_mesh_displacement_animation");
    gsParaviewCollection collectionMeshVel("cmc_mesh_velocity_animation");
    gsParaviewCollection collectionCurrentMesh("cmc_current_mesh_animation");
    
    // 初始化动画帧计数器
    index_t frameIndex = 0;
    
    for (index_t step = 0; step <= nTimeSteps; ++step)
    {
        real_t time = step * dt;
        gsInfo << "\n[Step " << step << "/" << nTimeSteps << "] t=" << time << "\n";
        
        // 调试：验证时间是否正确传递给运动函数
        real_t expected_disp = cylA * std::sin(2.0 * M_PI * cylF * time);
        gsInfo << "[Time Loop] Expected cylinder displacement: " << expected_disp << "\n";

        // 若使用显式壁面速度Dirichlet，这里需随时间更新；当前无需。

        if (step > 0)
        {
            // 根据频率开关本步是否执行参数化
            bool doParamThisStep = enableParam && (paramEveryN != 0) && (paramEveryN == 1 || (step % paramEveryN == 0));
            solver.setMeshOptimization(doParamThisStep);
            if (doParamThisStep)
                gsInfo << "[Param] Optimization enabled at step " << step << " (every " << paramEveryN << ")\n";
            solver.nextIteration();
        }

        if (step % outputInterval == 0)
        {
            // 构造并输出解场
            gsField<> velocity = solver.constructSolution(0);
            gsField<> pressure = solver.constructSolution(1);
            
            // 使用固定宽度的文件名格式，确保正确排序
            char vname[100], pname[100];
            sprintf(vname, "cmc_velocity_%03d", (int)frameIndex);
            sprintf(pname, "cmc_pressure_%03d", (int)frameIndex);
            
            // 写入单独文件
            gsWriteParaview(velocity, std::string(vname), samp, false);
            gsWriteParaview(pressure, std::string(pname), samp, false);
            
            // 添加到动画集合（按多patch写出的实际文件名：vname + patchIndex + .vts）
            for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
            {
                collectionVelocity.addPart(std::string(vname) + std::to_string(p) + ".vts", time, std::to_string(p));
                collectionPressure.addPart(std::string(pname) + std::to_string(p) + ".vts", time, std::to_string(p));
            }
            
            if (solver.isALEActive())
            {
                gsField<> meshDisp = solver.getMeshDisplacementField();
                gsField<> meshVel  = solver.getMeshVelocityField();
                char mname[100], vname2[100], gname[100];
                sprintf(mname,  "cmc_mesh_disp_%03d", (int)frameIndex);
                sprintf(vname2, "cmc_mesh_vel_%03d",  (int)frameIndex);
                sprintf(gname,  "cmc_current_mesh_%03d", (int)frameIndex);
                
                gsWriteParaview(meshDisp, std::string(mname), samp, false);
                gsWriteParaview(meshVel,  std::string(vname2), samp, false);
                for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
                {
                    collectionMeshDisp.addPart(std::string(mname)  + std::to_string(p) + ".vts", time, std::to_string(p));
                    collectionMeshVel .addPart(std::string(vname2) + std::to_string(p) + ".vts", time, std::to_string(p));
                }
                
                // 输出当前网格（变形后的几何）
                gsField<> currentMesh = solver.constructSolution(0); // 重用速度场的几何
                gsWriteParaview(currentMesh, std::string(gname), samp, false);
                for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
                    collectionCurrentMesh.addPart(std::string(gname) + std::to_string(p) + ".vts", time, std::to_string(p));
                
                gsInfo << "[Output] Frame " << frameIndex << " (Step " << step << ", t=" << time 
                       << "): files saved and added to animation collections\n";
            }
            
            frameIndex++;
        }
    }

    // 保存动画集合文件
    gsInfo << "\n=== Saving animation collections ===\n";
    collectionVelocity.save();
    collectionPressure.save();
    collectionMeshDisp.save();
    collectionCurrentMesh.save();
    collectionMeshVel.save();
    
    gsInfo << "模拟完成!\n";
    return 0;
}
