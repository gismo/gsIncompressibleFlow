# ALE Implementation Headers Documentation

本文档详细介绍了gsIncompressibleFlow模块中ALE（Arbitrary Lagrangian-Eulerian）实现的核心头文件。

## 目录
1. [gsINSTermsALE.h](#gsINSTermsALE) - ALE弱形式项
2. [gsINSVisitorsALE.h](#gsINSVisitorsALE) - ALE访问器
3. [gsINSAssemblerALE.h](#gsINSAssemblerALE) - ALE装配器
4. [gsINSSolverALE.h](#gsINSSolverALE) - ALE求解器
5. [类层次结构](#class-hierarchy)
6. [使用流程](#usage-flow)

## 类结构总览

```
┌─────────────────────────────────────────────────────────────────┐
│                     gsINSSolverUnsteadyALE                      │
│  主求解器类，管理整个ALE求解流程                                 │
│  - 网格更新函数 (m_meshUpdateFunc)                              │
│  - 网格优化控制 (m_useMeshOptimization)                         │
│  - 原始网格存储 (m_originalMesh)                                │
└────────────────────────┬────────────────────────────────────────┘
                         │ 包含
┌────────────────────────▼────────────────────────────────────────┐
│                  gsINSAssemblerUnsteadyALE                      │
│  装配器类，处理矩阵组装                                          │
│  - 网格速度系数 (m_meshVelCoefs)                                │
│  - 网格位移系数 (m_meshDispCoefs)                               │
│  - ALE访问器 (m_visitorUUnonlinALE)                            │
└────────────────────────┬────────────────────────────────────────┘
                         │ 包含
┌────────────────────────▼────────────────────────────────────────┐
│                  gsINSVisitorUUnonlinALE                        │
│  访问器类，管理ALE项                                            │
│  - ALE对流项 (m_aleConvectionTerm)                             │
│  - 网格速度场 (m_meshVelField)                                 │
└────────────────────────┬────────────────────────────────────────┘
                         │ 包含
┌────────────────────────▼────────────────────────────────────────┐
│                  gsINSTerm_ALEConvection                        │
│  ALE对流项，实现数学公式                                        │
│  - 计算相对速度 (u - u_mesh)                                   │
│  - 组装弱形式积分                                               │
└─────────────────────────────────────────────────────────────────┘
```

---

## 1. gsINSTermsALE.h {#gsINSTermsALE}

### 概述
定义了ALE对流项的弱形式实现，这是ALE方法的核心数学表达。

### 主要类：gsINSTerm_ALEConvection

#### 类定义 (行21-22)
```cpp
template <class T>
class gsINSTerm_ALEConvection : public gsFlowTermNonlin<T>
```

#### 核心功能
- **ALE对流项**：实现 `((u - u_mesh) · ∇φ_trial) * φ_test` 的弱形式 (行18)
- **相对速度计算**：`u_rel = u - u_mesh` (行96)
- **物理空间梯度转换**：将参数空间梯度转换到物理空间 (行119)

#### 关键成员变量 (行30-36)
- `m_tarDim`: 目标维度（通常为2或3）
- `m_meshVelVals`: 积分点上的网格速度值
- `m_meshVelField`: 网格速度场指针

#### 核心方法

##### computeMeshVelocity() (行57-70)
```cpp
void computeMeshVelocity(const gsMapData<T>& mapData)
```
- 在积分点上计算网格速度
- 若未设置网格速度场，默认为零 (行68)

##### assemble() (行78-134)
```cpp
virtual void assemble(const gsMapData<T>& mapData, 
                     const gsVector<T>& quWeights, 
                     const std::vector<gsMatrix<T>>& testFunData, 
                     const std::vector<gsMatrix<T>>& trialFunData, 
                     gsMatrix<T>& localMat)
```
- 组装ALE对流项的局部矩阵
- 计算相对速度 `relativeVel = m_solUVals - m_meshVelVals` (行96)
- 在物理空间计算梯度 (行119)
- 积分 `∫(φ_test * ((u-w)·∇φ_trial)) dΩ` (行131-132)

#### 设计特点
- **灵活性**：支持动态设置网格速度场
- **效率**：只在需要时计算网格速度
- **通用性**：支持2D和3D问题

---

## 2. gsINSVisitorsALE.h {#gsINSVisitorsALE}

### 概述
提供ALE访问器类，管理ALE对流项并集成到装配过程中。

### 主要类：gsINSVisitorUUnonlinALE

#### 类定义 (行22-24)
```cpp
template <class T, int MatOrder = RowMajor>
class gsINSVisitorUUnonlinALE : public gsINSVisitorUU<T, MatOrder>
```

#### 核心功能
- **管理ALE项**：创建和管理`gsINSTerm_ALEConvection`实例
- **接口适配**：连接ALE项与标准INS装配框架
- **状态管理**：更新网格速度场和当前解

#### 构造函数
提供两种构造方式：
1. **推荐方式**（带paramsPtr）(行43-58)：
```cpp
gsINSVisitorUUnonlinALE(typename gsFlowSolverParams<T>::Ptr paramsPtr,
                       const std::vector<gsDofMapper>& dofMappers,
                       index_t targetDim = 2,
                       const gsField<T>* meshVelField = nullptr)
```

2. **兼容方式**（不带paramsPtr）(行61-75)：
```cpp
gsINSVisitorUUnonlinALE(const std::vector<gsDofMapper>& dofMappers,
                       index_t targetDim = 2,
                       const gsField<T>* meshVelField = nullptr)
```

#### 关键方法

##### setMeshVelocityField() (行83-89)
```cpp
void setMeshVelocityField(const gsField<T>* meshVelField)
```
- 更新网格速度场
- 同步更新内部ALE对流项 (行88)

##### setCurrentSolution() (行92-100)
```cpp
void setCurrentSolution(const gsField<T>& solution)
```
- 设置当前速度解
- 传递给所有非线性项

#### 设计特点
- **继承复用**：最大化利用基类功能
- **封装性**：隐藏ALE项的创建细节
- **一致性**：保持与标准INS访问器相同的接口

---

## 3. gsINSAssemblerALE.h {#gsINSAssemblerALE}

### 概述
扩展非稳态INS装配器，添加ALE功能，处理网格运动和速度计算。

### 主要类：gsINSAssemblerUnsteadyALE

#### 类定义 (行21-23)
```cpp
template <class T, int MatOrder = RowMajor>
class gsINSAssemblerUnsteadyALE : public gsINSAssemblerUnsteady<T, MatOrder>
```

#### 核心功能
- **网格管理**：存储和更新网格位移/速度
- **ALE装配**：集成ALE访问器进行矩阵组装
- **场构造**：创建网格速度和位移场

#### 关键成员变量 (行39-51)
```cpp
gsMatrix<T> m_meshVelCoefs;      // 网格速度系数 (行40)
gsMatrix<T> m_meshDispCoefs;     // 当前网格位移 (行41)
gsMatrix<T> m_meshDispOld;       // 上一时刻网格位移 (行42)
gsINSVisitorUUnonlinALE<T, MatOrder>* m_visitorUUnonlinALE;  // ALE访问器 (行45)
bool m_isALEActive;              // ALE激活标志 (行48)
gsField<T>* m_tempMeshVelField;  // 临时网格速度场 (行51)
```

#### 核心方法

##### initialize() (行70-86)
```cpp
virtual void initialize() override
```
- 初始化基类 (行73)
- 分配网格系数存储空间 (行76-79)
- 创建ALE访问器 (行82-83)

##### updateMesh() (行108-131)
```cpp
void updateMesh(const gsMatrix<T>& meshDispNew)
```
- 计算网格速度：`v_mesh = (disp_new - disp_old) / dt` (行123)
- 更新位移存储 (行130)
- 扩展向量到完整尺寸（包含压力部分为零）(行121-127)

##### assembleNonlinearPart() (行154-195)
```cpp
virtual void assembleNonlinearPart() override
```
- 如果ALE激活，使用ALE访问器 (行156)
- 创建临时网格速度场 (行159-160)
- 组装ALE对流项 (行181, 187)

##### getMeshVelocityField() / getMeshDisplacementField() (行134-143)
```cpp
gsField<T> getMeshVelocityField() const      // 行134-137
gsField<T> getMeshDisplacementField() const  // 行140-143
```
- 从系数构造完整的场表示
- 用于可视化和数据交换

#### 设计特点
- **增量更新**：只存储自由度的位移/速度
- **内存效率**：按需创建临时场
- **向后兼容**：ALE未激活时退化为标准装配器

---

## 4. gsINSSolverALE.h {#gsINSSolverALE}

### 概述
最高层的ALE求解器类，提供完整的ALE流动求解功能，包括网格优化和FSI耦合支持。

### 主要类：gsINSSolverUnsteadyALE

#### 类定义 (行24-26)
```cpp
template <class T = real_t, int MatOrder = RowMajor>
class gsINSSolverUnsteadyALE : public gsINSSolverUnsteady<T, MatOrder>
```

#### 核心功能
- **完整ALE求解**：时间步进、网格更新、流动求解
- **网格优化**：集成gsBarrierPatch
- **动态边界**：支持旋转域
- **FSI准备**：提供流固耦合接口

#### 关键成员变量 (行38-61)
```cpp
bool m_isALEActive;                          // ALE激活标志 (行38)
std::function<gsMatrix<T>(T)> m_meshUpdateFunc;  // 网格更新函数 (行41)
bool m_useMeshOptimization;                  // 网格优化标志 (行44)
bool m_useDynamicBoundaryMapping;            // 动态边界映射标志 (行47)
gsMultiPatch<T> m_originalMesh;              // 原始网格 (行53)
gsMatrix<T> m_previousDisp;                  // 前一时刻位移 (行56)
T m_rotationPeriod;                          // 旋转周期 (行59)
gsVector<T> m_rotationCenter;                // 旋转中心 (行60)
```

#### 核心方法

##### setALEActive() (行87-99)
```cpp
void setALEActive(bool active)
```
- 激活/停用ALE (行91)
- 首次激活时存储原始网格 (行96)
- 同步更新装配器状态 (行91)

##### setMeshUpdateFunction() (行105-109)
```cpp
void setMeshUpdateFunction(std::function<gsMatrix<T>(T)> func)
```
- 设置网格运动规律
- 函数输入：时间t
- 函数输出：网格位移

##### nextIteration() (行183-200)
```cpp
virtual void nextIteration() override
```
- ALE时间步进主函数
- 调用顺序：
  1. 应用网格位移 (行189)
  2. 网格优化（可选）(行194)
  3. 求解流动问题 (行199)

##### applyMeshDisplacement() (行203-244)
```cpp
void applyMeshDisplacement()
```
- 将位移应用到实际几何
- 计算增量位移 (行241)
- 更新装配器中的网格 (行242)

##### optimizeMesh() (行247-290)
```cpp
void optimizeMesh()
```
- 使用gsBarrierPatch优化网格质量 (行275)
- 支持标准和动态边界映射模式 (行254)
- 异常处理：优化失败时继续计算 (行286-289)

#### 高级功能

##### 动态边界映射
- 用于旋转域问题
- 自动计算当前旋转角度
- 集成gsBarrierPatchDynamic

##### 网格优化选项
```cpp
gsOptionList m_meshOptOptions;
```
- Verbose：输出详细程度
- ParamMethod：参数化方法
- AAPreconditionType：预条件类型

### 辅助类：gsFSIHelper (行296-387)

#### 概述
为流固耦合提供工具函数（框架预留，待实现）。

#### 主要功能
- **接口定义**：设置流固界面
- **位移传递**：结构到流体网格
- **力传递**：流体到结构
- **谐波扩展**：边界位移到域内

#### 静态工具方法 (行365-370)
```cpp
static gsMatrix<T> computeMeshVelocity(const gsMatrix<T>& dispNew,
                                      const gsMatrix<T>& dispOld,
                                      T dt)
```
- 从位移历史计算网格速度
- 使用后向差分格式 (行369)

---

## 5. 类层次结构 {#class-hierarchy}

### 继承关系图
```
gsFlowTermNonlin<T>
    └── gsINSTerm_ALEConvection<T>

gsINSVisitorUU<T, MatOrder>
    └── gsINSVisitorUUnonlinALE<T, MatOrder>

gsINSAssemblerUnsteady<T, MatOrder>
    └── gsINSAssemblerUnsteadyALE<T, MatOrder>

gsINSSolverUnsteady<T, MatOrder>
    └── gsINSSolverUnsteadyALE<T, MatOrder>

独立辅助类：
    gsFSIHelper<T>
```

### 详细类结构图

```
┌──────────────────────────────────────────────────────────────┐
│  gsINSTerm_ALEConvection<T>                                  │
├──────────────────────────────────────────────────────────────┤
│ - m_tarDim: index_t                                          │
│ - m_meshVelVals: gsMatrix<T>                                 │
│ - m_meshVelField: const gsField<T>*                         │
├──────────────────────────────────────────────────────────────┤
│ + setMeshVelocityField(field): void                          │
│ + computeMeshVelocity(mapData): void                         │
│ + assemble(mapData, quWeights, ...): void override          │
└──────────────────────────────────────────────────────────────┘
                            ▲
                            │ 被包含
┌──────────────────────────────────────────────────────────────┐
│  gsINSVisitorUUnonlinALE<T, MatOrder>                       │
├──────────────────────────────────────────────────────────────┤
│ - m_tarDim: index_t                                          │
│ - m_meshVelField: const gsField<T>*                         │
│ - m_aleConvectionTerm: gsINSTerm_ALEConvection<T>*          │
├──────────────────────────────────────────────────────────────┤
│ + setMeshVelocityField(field): void                          │
│ + setCurrentSolution(solution): void                         │
└──────────────────────────────────────────────────────────────┘
                            ▲
                            │ 被包含
┌──────────────────────────────────────────────────────────────┐
│  gsINSAssemblerUnsteadyALE<T, MatOrder>                     │
├──────────────────────────────────────────────────────────────┤
│ - m_meshVelCoefs: gsMatrix<T>                                │
│ - m_meshDispCoefs: gsMatrix<T>                               │
│ - m_meshDispOld: gsMatrix<T>                                 │
│ - m_visitorUUnonlinALE: gsINSVisitorUUnonlinALE<T>*         │
│ - m_isALEActive: bool                                        │
│ - m_tempMeshVelField: gsField<T>*                           │
├──────────────────────────────────────────────────────────────┤
│ + initialize(): void override                                 │
│ + updateMesh(meshDispNew): void                              │
│ + assembleNonlinearPart(): void override                     │
│ + getMeshVelocityField(): gsField<T>                        │
│ + getMeshDisplacementField(): gsField<T>                    │
└──────────────────────────────────────────────────────────────┘
                            ▲
                            │ 被包含
┌──────────────────────────────────────────────────────────────┐
│  gsINSSolverUnsteadyALE<T, MatOrder>                        │
├──────────────────────────────────────────────────────────────┤
│ - m_isALEActive: bool                                        │
│ - m_meshUpdateFunc: std::function<gsMatrix<T>(T)>           │
│ - m_useMeshOptimization: bool                                │
│ - m_useDynamicBoundaryMapping: bool                          │
│ - m_originalMesh: gsMultiPatch<T>                           │
│ - m_previousDisp: gsMatrix<T>                                │
│ - m_rotationPeriod: T                                        │
│ - m_rotationCenter: gsVector<T>                             │
├──────────────────────────────────────────────────────────────┤
│ + setALEActive(active): void                                 │
│ + setMeshUpdateFunction(func): void                          │
│ + nextIteration(): void override                             │
│ + applyMeshDisplacement(): void                              │
│ + optimizeMesh(): void                                       │
│ + setMeshOptimization(enable): void                          │
│ + setRotationParameters(period, center): void                │
└──────────────────────────────────────────────────────────────┘
```

### 设计模式
- **继承扩展**：通过继承添加ALE功能
- **组合模式**：求解器包含装配器，装配器包含访问器
- **策略模式**：可切换标准/ALE模式

---

## 6. 使用流程 {#usage-flow}

### 基本使用步骤

1. **创建求解器**
```cpp
gsINSSolverUnsteadyALE<> solver(paramsPtr);
solver.initialize();
```

2. **激活ALE**
```cpp
solver.setALEActive(true);
```

3. **设置网格运动**
```cpp
solver.setMeshUpdateFunction([](real_t t) {
    return computeDisplacement(t);
});
```

4. **可选：启用网格优化**
```cpp
solver.setMeshOptimization(true);
solver.getMeshOptOptions().setInt("Verbose", 1);
```

5. **时间步进**
```cpp
for (int step = 0; step < nSteps; ++step) {
    solver.nextIteration();
}
```

### 数据流
```
用户定义网格运动
    ↓
gsINSSolverUnsteadyALE::nextIteration()
    ↓
applyMeshDisplacement()
    ↓
gsINSAssemblerUnsteadyALE::updateMesh()
    ↓
gsINSAssemblerUnsteadyALE::assembleNonlinearPart()
    ↓
gsINSVisitorUUnonlinALE (with mesh velocity)
    ↓
gsINSTerm_ALEConvection::assemble()
    ↓
组装ALE对流项矩阵
```

### 关键设计决策

1. **最小侵入性**：通过继承而非修改基类
2. **灵活性**：支持动态切换ALE/标准模式
3. **效率**：只在需要时创建额外数据结构
4. **可扩展性**：预留FSI和其他耦合接口
5. **鲁棒性**：网格优化失败时的容错处理

---

## 总结

这套ALE实现提供了一个完整、高效、灵活的移动网格流动求解框架。通过分层设计，从底层的弱形式项到高层的求解器，每一层都有明确的职责和接口。这种设计使得代码易于理解、维护和扩展，特别适合处理流固耦合、旋转机械等涉及移动边界的复杂流动问题。