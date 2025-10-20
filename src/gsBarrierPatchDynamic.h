/** @file gsBarrierPatchDynamic.h
    @brief Helper for dynamic boundary mapping during mesh optimization
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsModeling/gsBarrierPatch.h>
#include <gsCore/gsBoundary.h>

namespace gismo
{

/**
 * \brief Helper class for dynamic boundary correspondence in rotating domains
 * 
 * This class provides utilities to handle changing boundary mappings
 * during rotation, particularly for ALE methods with rotating objects.
 * 
 * \tparam d domain dimension
 * \tparam T Coefficient type
 */
template<short_t d, typename T=real_t>
class gsBarrierPatchDynamic
{
public:
    /// Constructor
    gsBarrierPatchDynamic() : m_rotationAngle(0), m_dynamicBoundaryMapping(false) 
    { 
        m_rotationCenter.setZero(d);
    }
    
    /// Enable/disable dynamic boundary mapping
    void setDynamicBoundaryMapping(bool enable) { m_dynamicBoundaryMapping = enable; }
    
    /// Set rotation angle
    void setRotationAngle(T angle) { m_rotationAngle = angle; }
    
    /// Set rotation center
    void setRotationCenter(const gsVector<T>& center) { m_rotationCenter = center; }
    
    /// Compute optimized mesh with dynamic boundary handling
    gsMultiPatch<T> compute(const gsMultiPatch<T>& mp, const gsOptionList& options)
    {
        if (m_dynamicBoundaryMapping && std::abs(m_rotationAngle) > 0)
        {
            // For rotating domains, we need to handle boundaries differently
            // based on the rotation state
            
            // Determine which boundaries are currently close to each other
            std::vector<std::pair<patchSide, patchSide>> dynamicInterfaces;
            findDynamicInterfaces(mp, dynamicInterfaces);
            
            // Create gsBarrierPatch with appropriate settings
            // For now, use patch-wise optimization to avoid interface issues
            gsBarrierPatch<d, T> opt(mp, false); // true = patch-wise
            opt.options() = options;
            
            // Log dynamic interfaces for debugging
            if (!dynamicInterfaces.empty() && options.askInt("Verbose", 0) > 0)
            {
                gsInfo << "Dynamic interfaces detected after rotation of " 
                       << m_rotationAngle << " radians:\n";
                for (const auto& iface : dynamicInterfaces)
                {
                    gsInfo << "  Patch " << iface.first.patch << " side " << iface.first.side() 
                           << " <-> Patch " << iface.second.patch << " side " << iface.second.side() << "\n";
                }
            }
            
            opt.compute();
            return opt.result();
        }
        else
        {
            // Standard optimization without dynamic handling
            gsBarrierPatch<d, T> opt(mp, true); // false = global optimization
            opt.options() = options;
            opt.compute();
            return opt.result();
        }
    }
    
protected:
    /// Find boundaries that are close after rotation
    void findDynamicInterfaces(const gsMultiPatch<T>& mp, 
                               std::vector<std::pair<patchSide, patchSide>>& interfaces)
    {
        interfaces.clear();
        const T tol = 1e-6;
        
        // For each pair of patches, check if their boundaries are close
        for (size_t p1 = 0; p1 < mp.nPatches(); ++p1)
        {
            for (index_t s1 = 1; s1 <= 2*d; ++s1)
            {
                patchSide ps1(p1, boxSide(s1));
                
                // Sample points on this boundary
                gsMatrix<T> pts1;
                sampleBoundary(mp.patch(p1), boxSide(s1), pts1);
                
                // Check against other boundaries
                for (size_t p2 = p1; p2 < mp.nPatches(); ++p2)
                {
                    for (index_t s2 = (p1==p2 ? s1+1 : 1); s2 <= 2*d; ++s2)
                    {
                        patchSide ps2(p2, boxSide(s2));
                        
                        // Sample points on the other boundary
                        gsMatrix<T> pts2;
                        sampleBoundary(mp.patch(p2), boxSide(s2), pts2);
                        
                        // Check if boundaries are close
                        if (areBoundariesClose(pts1, pts2, tol))
                        {
                            interfaces.push_back({ps1, ps2});
                        }
                    }
                }
            }
        }
    }
    
    /// Sample points on a boundary
    void sampleBoundary(const gsGeometry<T>& patch, boxSide side, gsMatrix<T>& pts)
    {
        const index_t nSamples = 5;
        
        // Get parameter bounds for the boundary
        gsMatrix<T> params(d-1, nSamples);
        for (index_t i = 0; i < nSamples; ++i)
        {
            params(0, i) = i / T(nSamples - 1);
        }
        
        // Map to full parameter space
        gsMatrix<T> fullParams(d, nSamples);
        for (index_t i = 0; i < nSamples; ++i)
        {
            for (index_t j = 0; j < d; ++j)
            {
                if (side.direction() == j)
                {
                    fullParams(j, i) = (side.parameter() == 0) ? T(0) : T(1);
                }
                else
                {
                    fullParams(j, i) = params(j < side.direction() ? j : j-1, i);
                }
            }
        }
        
        // Evaluate patch at boundary points
        patch.eval_into(fullParams, pts);
    }
    
    /// Check if two sets of boundary points are close
    bool areBoundariesClose(const gsMatrix<T>& pts1, const gsMatrix<T>& pts2, T tol)
    {
        // For each point in pts1, find the closest point in pts2
        for (index_t i = 0; i < pts1.cols(); ++i)
        {
            T minDist = std::numeric_limits<T>::max();
            for (index_t j = 0; j < pts2.cols(); ++j)
            {
                T dist = (pts1.col(i) - pts2.col(j)).norm();
                minDist = std::min(minDist, dist);
            }
            
            if (minDist > tol)
                return false;
        }
        
        return true;
    }
    
private:
    T m_rotationAngle;
    gsVector<T> m_rotationCenter;
    bool m_dynamicBoundaryMapping;
};

} // namespace gismo