// 压力诊断工具
#include <gismo.h>

using namespace gismo;

template<class T>
void diagnosePressureField(const gsField<T>& pressureField,
                          const gsMultiPatch<T>& fluidGeometry,
                          const std::vector<std::pair<index_t, boundary::side>>& interfaces)
{
    gsInfo << "\n=== 压力场诊断 ===\n";
    
    // 1. 全局压力统计
    T globalMin = 1e10, globalMax = -1e10, globalSum = 0;
    index_t totalPoints = 0;
    
    for (size_t p = 0; p < pressureField.patches().nPatches(); ++p)
    {
        // 在patch内部采样
        gsMatrix<T> pts = gsPointGrid<T>(pressureField.patches().patch(p).support(), 10);
        gsMatrix<T> vals;
        pressureField.function(p).eval_into(pts, vals);
        
        T patchMin = vals.minCoeff();
        T patchMax = vals.maxCoeff();
        T patchAvg = vals.mean();
        
        globalMin = std::min(globalMin, patchMin);
        globalMax = std::max(globalMax, patchMax);
        globalSum += vals.sum();
        totalPoints += vals.cols();
        
        gsInfo << "Patch " << p << " 压力统计:\n";
        gsInfo << "  范围: [" << patchMin << ", " << patchMax << "]\n";
        gsInfo << "  平均: " << patchAvg << "\n";
    }
    
    gsInfo << "\n全局压力统计:\n";
    gsInfo << "  范围: [" << globalMin << ", " << globalMax << "]\n";
    gsInfo << "  平均: " << globalSum / totalPoints << "\n";
    
    // 2. 边界压力统计
    gsInfo << "\n边界压力统计:\n";
    for (size_t i = 0; i < interfaces.size(); ++i)
    {
        index_t patchId = interfaces[i].first;
        boundary::side side = interfaces[i].second;
        
        // 获取边界几何
        gsGeometry<T>* bndGeom = fluidGeometry.patch(patchId).boundary(side);
        gsMatrix<T> support = bndGeom->support();
        
        // 在边界上采样
        gsMatrix<T> bndParams(1, 20);
        for (index_t j = 0; j < 20; ++j)
        {
            bndParams(0, j) = support(0,0) + (support(0,1) - support(0,0)) * j / 19.0;
        }
        
        // 转换为patch参数
        gsMatrix<T> patchParams(2, 20);
        for (index_t j = 0; j < 20; ++j)
        {
            switch (side)
            {
                case boundary::south:
                    patchParams(0, j) = bndParams(0, j);
                    patchParams(1, j) = 0.0;
                    break;
                case boundary::north:
                    patchParams(0, j) = bndParams(0, j);
                    patchParams(1, j) = 1.0;
                    break;
                case boundary::west:
                    patchParams(0, j) = 0.0;
                    patchParams(1, j) = bndParams(0, j);
                    break;
                case boundary::east:
                    patchParams(0, j) = 1.0;
                    patchParams(1, j) = bndParams(0, j);
                    break;
            }
        }
        
        // 评估压力
        gsMatrix<T> pressureVals;
        pressureField.function(patchId).eval_into(patchParams, pressureVals);
        
        gsInfo << "界面 " << i << " (patch " << patchId << ", side " << side << "):\n";
        gsInfo << "  压力范围: [" << pressureVals.minCoeff() << ", " << pressureVals.maxCoeff() << "]\n";
        gsInfo << "  平均压力: " << pressureVals.mean() << "\n";
        
        // 计算压力梯度
        T maxGrad = 0;
        for (index_t j = 1; j < pressureVals.cols(); ++j)
        {
            T grad = math::abs(pressureVals(0, j) - pressureVals(0, j-1));
            maxGrad = std::max(maxGrad, grad);
        }
        gsInfo << "  最大压力梯度: " << maxGrad << "\n";
        
        delete bndGeom;
    }
    
    // 3. 检查压力单位和合理性
    gsInfo << "\n压力合理性检查:\n";
    T typicalFluidDensity = 1000.0;  // kg/m^3 (水)
    T typicalVelocity = 1.0;         // m/s
    T typicalPressure = 0.5 * typicalFluidDensity * typicalVelocity * typicalVelocity;
    
    gsInfo << "典型动压（0.5*rho*v^2）: " << typicalPressure << " Pa\n";
    gsInfo << "压力量级比: " << globalMax / typicalPressure << "\n";
    
    if (globalMax > 100 * typicalPressure)
    {
        gsWarn << "警告：压力值可能过大！考虑检查：\n";
        gsWarn << "  1. 压力单位是否正确（Pa vs kPa）\n";
        gsWarn << "  2. 无量纲化是否正确\n";
        gsWarn << "  3. 边界条件是否合理\n";
    }
}

// 压力缩放建议
template<class T>
T suggestPressureScaling(const gsField<T>& pressureField,
                        T targetMaxPressure = 100.0)  // Pa
{
    T maxPressure = 0;
    
    for (size_t p = 0; p < pressureField.patches().nPatches(); ++p)
    {
        gsMatrix<T> pts = gsPointGrid<T>(pressureField.patches().patch(p).support(), 10);
        gsMatrix<T> vals;
        pressureField.function(p).eval_into(pts, vals);
        maxPressure = std::max(maxPressure, vals.cwiseAbs().maxCoeff());
    }
    
    T scaling = 1.0;
    if (maxPressure > targetMaxPressure)
    {
        scaling = targetMaxPressure / maxPressure;
        gsInfo << "\n建议的压力缩放因子: " << scaling << "\n";
        gsInfo << "这将把最大压力从 " << maxPressure << " 降低到 " << targetMaxPressure << "\n";
    }
    
    return scaling;
}