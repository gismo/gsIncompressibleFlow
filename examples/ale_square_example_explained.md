# ALE Square Example 详细解释文档

这个示例展示了如何使用G+Smo中的ALE（Arbitrary Lagrangian-Eulerian）求解器来模拟带有移动网格的不可压缩流动问题。

## 文件头部说明 (1-16行)

```cpp
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
```

**说明**：
- 文件注释表明这是一个用于流固耦合（FSI）仿真的ALE示例
- 包含了必要的头文件：主G+Smo库、INS求解器和ALE求解器
- 使用gismo命名空间简化代码

## 物理参数设置 (17-34行)

```cpp
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
```

**参数解释**：
- **Re = 100**: 雷诺数，控制流动的惯性力与粘性力之比
- **rho = 1.0**: 流体密度（标准化值）
- **mu = rho/Re**: 动力粘度，由雷诺数定义计算得出
- **U = 1.0**: 特征速度（用于入口边界条件）
- **L = 1.0**: 特征长度（域的尺寸）
- **T = 10.0**: 总仿真时间
- **dt = 0.01**: 时间步长
- **numRefine = 2**: 网格细化次数
- **degree = 2**: B样条基函数的阶数

## 输出基本信息 (35-39行)

```cpp
    gsInfo << "=== ALE FSI Example ===\n";
    gsInfo << "Reynolds number: " << Re << "\n";
    gsInfo << "Time step: " << dt << "\n";
    gsInfo << "Number of time steps: " << nTimeSteps << "\n\n";
```

**功能**：打印仿真的基本参数信息

## 创建计算域 (40-64行)

### 固体域创建 (41-45行)
```cpp
    // Create artificial solid (unit square)
    gsMultiPatch<> solidDomain;
    solidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1,0,0));

    gsMultiBasis<> solidBasis(solidDomain);
```

**说明**：创建一个单位正方形作为"人工固体"域（虽然在本例中未使用）

### 流体域创建 (46-64行)
```cpp
    // Create fluid domain (unit square)
    gsMultiPatch<> fluidDomain;
    // Upper row
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1,-1,1));
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1,0, 1));
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1,1));
    // Middle row
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1,-1,0));
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1,1,0));
    
    // Lower row
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1,-1,-1));
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1,0,-1));
    fluidDomain.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1,-1));

    fluidDomain.computeTopology();
    
    gsWriteParaview(fluidDomain, "fluid_domain", 1000);
```

**域的布局**：
```
[P0] [P1] [P2]  <- 上排 (y=1)
[P3] [空] [P4]  <- 中排 (y=0)，中间有孔
[P5] [P6] [P7]  <- 下排 (y=-1)
```

**关键点**：
- 创建了8个patch，形成一个3×3网格，中间缺少一个patch（形成孔洞）
- `computeTopology()`: 计算patch之间的连接关系
- `gsWriteParaview()`: 输出初始网格到ParaView文件

## 网格细化和基函数设置 (66-73行)

```cpp
    // Refine mesh
    for (index_t i = 0; i < numRefine; ++i)
        fluidDomain.uniformRefine();
    
    // Create basis
    gsMultiBasis<> basis(fluidDomain);
    basis.setDegree(degree);
```

**功能**：
- 对所有patch进行均匀细化（2次）
- 创建多patch基函数，设置阶数为2

## 边界条件设置 (74-103行)

### 定义边界函数 (77-79行)
```cpp
    // Velocity boundary conditions
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    gsConstantFunction<> inletVel(U, 0.0, 2);
```

**说明**：
- `zeroVel`: 零速度(0,0)，用于固壁边界
- `inletVel`: 入口速度(U,0)=(1,0)，水平流入

### 设置边界条件 (81-95行)
```cpp
    // Fixed walls (top and bottom)
    bcInfo.addCondition(6, boundary::south, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(5, boundary::south, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(4, boundary::south, condition_type::dirichlet, &zeroVel, 0);

    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, &zeroVel, 0);
    
    // Inlet (left)
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &inletVel, 0);
    bcInfo.addCondition(7, boundary::west, condition_type::dirichlet, &inletVel, 0);
    bcInfo.addCondition(6, boundary::west, condition_type::dirichlet, &inletVel, 0);
```

**边界条件布局**：
- 上下壁面：无滑移条件（零速度）
- 左侧入口：均匀流入（速度=1）
- 右侧出口：自然边界条件（未显式设置）
- 内部孔洞边界：作为FSI界面，将动态更新

## 创建PDE和求解器 (105-134行)

### 创建Navier-Stokes PDE (107行)
```cpp
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &zeroVel, mu);
```

**参数**：
- `fluidDomain`: 计算域
- `bcInfo`: 边界条件
- `&zeroVel`: 体积力（这里为零）
- `mu`: 动力粘度

### 设置离散基函数 (110-116行)
```cpp
    std::vector<gsMultiBasis<>> discreteBases;
    discreteBases.push_back(basis);  // Velocity basis
    discreteBases.push_back(basis);  // Pressure basis
    
    // For Taylor-Hood elements, elevate velocity space
    discreteBases[0].degreeElevate(1);
```

**Taylor-Hood元素**：
- 速度使用P3（阶数提升1）
- 压力使用P2
- 满足inf-sup条件，保证数值稳定性

### 设置求解器参数 (117-126行)
```cpp
    gsFlowSolverParams<real_t> params(nsPde, discreteBases);
    
    // Solver options
    params.options().setInt("nonlin.maxIt", 100);     // Max Picard iterations
    params.options().setReal("nonlin.tol", 1e-6);     // Picard tolerance
    params.options().setInt("lin.maxIt", 200);        // Max linear solver iterations
    params.options().setReal("lin.tol", 1e-8);        // Linear solver tolerance
    params.options().setReal("timeStep", dt);         // Time step size
    params.options().setSwitch("quiet", false);       // Enable verbose output
```

**参数说明**：
- 非线性求解：最多100次Picard迭代，容差1e-6
- 线性求解：最多200次迭代，容差1e-8
- 启用详细输出

### 创建并初始化ALE求解器 (128-134行)
```cpp
    // Create ALE solver
    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));
    
    // Initialize solver
    solver.initialize();
    
    // Activate ALE after initialization
    solver.setALEActive(true);
```

## 定义网格运动 (136-167行)

```cpp
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
```

**网格运动说明**：
- 定义了一个lambda函数，实现圆周运动
- 半径R=0.05，角频率w=2π（频率1Hz）
- 只移动内部点，边界点保持固定
- 使用DOF映射器确保只更新自由度

## 求解初始Stokes问题 (169-171行)

```cpp
    // Solve Stokes problem for initial conditions
    gsInfo << "Solving initial Stokes problem...\n";
    solver.solveStokes();
```

**功能**：求解稳态Stokes问题作为初始条件

## 时间步进循环 (173-217行)

```cpp
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
```

**每个时间步的操作**：
1. 调用`nextIteration()`推进一个时间步
2. 提取速度和压力系数
3. 获取网格位移场
4. 输出统计信息（最大值范数）

### 周期性输出结果 (195-216行)
```cpp
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
```

**输出说明**：
- 每10个时间步输出一次
- 输出速度场、压力场和网格位移场
- 使用1000×1000的采样密度
- 文件名包含时间步编号

## 程序结束 (219-223行)

```cpp
    gsInfo << "\nSimulation completed!\n";
    
    // No need to delete - nsPde is a stack object
    return 0;
}
```

## 总结

这个示例展示了ALE方法的核心功能：
1. **多patch域设置**：创建带孔洞的复杂几何
2. **Taylor-Hood元素**：确保数值稳定性
3. **动态网格**：通过lambda函数定义网格运动
4. **时间积分**：使用隐式时间步进
5. **可视化输出**：生成ParaView文件便于后处理

这是一个完整的ALE流动仿真框架，可以作为更复杂FSI问题的基础。