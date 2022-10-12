/** @file uwbSliceFieldCreator.h

Author(s): H. Hornikova
*/

#pragma once

#include <gismo.h>
#include <string>

#include <gsTensor/gsTensorTools.h>

#include "uwbINSBlockAssembler.h"

#include "uwbTurbineUtils.h"

namespace gismo
{

template<class T>
class uwbSliceFieldCreator
{
public:
    // swap: we have directions u, v, w in parametric space, the slice has two of them,
    //       for slice rotation, we need the radial direction to be first and the circumferential to be second,
    //       swap = true if the two directions of the slice have to be swapped to satisfy this

    uwbSliceFieldCreator(std::string inputFile, index_t knot, index_t numBlades, index_t dir = 2, bool swap = true)
        : m_knot(knot), m_numBlades(numBlades), m_dir(dir), m_swap(swap), m_compute(false)
    {
        initialize(inputFile);
    }

    uwbSliceFieldCreator(gsTensorBSpline<3, T> patch, gsTensorBSplineBasis<3, T> solBasis, gsMatrix<T> solCoefs, index_t knot, index_t numBlades, index_t dir = 2, bool swap = true)
        : m_knot(knot), m_numBlades(numBlades), m_dir(dir), m_swap(swap), m_compute(false)
    {
        initialize(patch, solBasis, solCoefs);
    }

    ~uwbSliceFieldCreator()
    { } 

protected:
    void initialize(std::string inputFile)
    {
        gsFileData<> fileData(inputFile);

        gsTensorBSpline<3, T> patch;
        gsTensorBSplineBasis<3, T> solBasis;
        gsMatrix<T> solCoefs;

        // read patch geometry
        if (fileData.has< gsTensorBSpline<3, T> >())
        {
            patch = *(fileData.getFirst< gsTensorBSpline<3, T> >());
        }
        else
        {
            GISMO_ERROR("Input file doesn't have a gsTensorBSpline<3, T> inside.\n");
        }

        // read discrete basis for the patch
        if (fileData.has< gsTensorBSplineBasis<3, T> >())
        {
            solBasis = *(fileData.getFirst< gsTensorBSplineBasis<3, T> >());
        }
        else
        {
            GISMO_ERROR("Input file doesn't have a gsTensorBSplineBasis<3, T> inside.\n");
        }

        // read solution coefficients
        if (fileData.has< gsMatrix<T> >())
        {
            solCoefs = *(fileData.getFirst< gsMatrix<T> >());
        }
        else
        {
            GISMO_ERROR("Input file doesn't have a gsMatrix inside.\n");
        }

        initialize(patch, solBasis, solCoefs);
    }

    void initialize(gsTensorBSpline<3, T> patch, gsTensorBSplineBasis<3, T> solBasis, gsMatrix<T> solCoefs)
    {
        // extract geometry slice 
        patch.slice(m_dir, patch.knots(m_dir).at(m_knot), m_origSlice);

        // extract solution basis slice
        m_origBasis = gsTensorBSplineBasis<2, T>(solBasis.component((m_dir + 1) % 3).knots(), solBasis.component((m_dir + 2) % 3).knots());

        // extract coefficients slice
        m_origCoefs = solCoefs.middleRows(solBasis.coefSlice(m_dir, m_knot).at(0), solBasis.coefSlice(m_dir, m_knot).rows());

        m_rotCoefs.setZero(m_numBlades * getNumSliceCoefs(), 3);
    }

    void rotatePatch()
    {
        gsVector<real_t, 3> axis;
        axis << 1, 0, 0;

        if(m_swap)
            m_origSlice.swapDirections(0, 1);

        gsKnotVector<T> kv1 = m_origSlice.knots(0);
        gsKnotVector<T> kv2 = m_origSlice.knots(1);

        size_t numCoefs = m_origSlice.coefs().rows();

        gsKnotVector<T> kv2new = kv2;
        gsMatrix<T> newCoefs(m_numBlades * numCoefs, 3);
        newCoefs.setZero();

        newCoefs.middleRows(0, numCoefs) = m_origSlice.coefs();

        gsTensorBSpline<2, T> other = *(m_origSlice.clone());

        for (int i = 1; i < m_numBlades; i++)
        {
            std::vector<real_t> tmp;
            for (int j = kv2.multLast(); j < kv2.size(); j++)
                tmp.push_back(kv2.at(j) + kv2new.last());

            other.rotate(-getPhi(), axis);

            newCoefs.middleRows(i*numCoefs, numCoefs) = other.coefs();

            gsKnotVector<T> tmpKV(tmp);
            kv2new = kv2new.knotUnion(tmpKV);
        }

        m_rotSlice = gsTensorBSpline<2, T>(kv1, kv2new, newCoefs);
    }

    void rotateBasis()
    {

        gsKnotVector<T> kv1 = m_origBasis.knots(m_swap ? 1 : 0);
        gsKnotVector<T> kv2 = m_origBasis.knots(m_swap ? 0 : 1);

        gsKnotVector<T> kv2new = kv2;

        for (int i = 1; i < m_numBlades; i++)
        {
            std::vector<real_t> tmp;
            for (int j = kv2.multLast(); j < kv2.size(); j++)
                tmp.push_back(kv2.at(j) + kv2new.last());

            gsKnotVector<T> tmpKV(tmp);
            kv2new = kv2new.knotUnion(tmpKV);
        }

        m_rotBasis = gsTensorBSplineBasis<2, T>(kv1, kv2new);
    }

    void rotateCoefs()
    {
        if (m_swap)
        {
            gsVector<int, 2> basisSizes;
            m_origBasis.size_cwise(basisSizes);
            swapTensorDirection(0, 1, basisSizes, m_origCoefs);
        }        

        gsMatrix<real_t> transformMatrix(3, 3);

        const real_t cos = math::cos(getPhi());
        const real_t sin = math::sin(getPhi());
        transformMatrix(0, 0) = 1;
        transformMatrix(0, 1) = 0;
        transformMatrix(0, 2) = 0;
        transformMatrix(1, 0) = 0;
        transformMatrix(1, 1) = cos;
        transformMatrix(1, 2) = sin;
        transformMatrix(2, 0) = 0;
        transformMatrix(2, 1) = -sin;
        transformMatrix(2, 2) = cos;

        m_rotCoefs.middleRows(0, getNumSliceCoefs()) = m_origCoefs;

        for (int i = 1; i < m_numBlades; i++)
        {
            for (int j = 0; j < getNumSliceCoefs(); j++)
            {
                m_rotCoefs.row(i*getNumSliceCoefs() + j) = transformMatrix * m_rotCoefs.row((i - 1)*getNumSliceCoefs() + j).transpose();
            }
        }
    }

    void areaAverage(int nPoints)
    {
        //gsInfo << "In area average ...\n";

        T area = 0;

        gsMatrix<T> Uaux(1,3);
        Uaux.setZero();
        gsMatrix<T> Uaveraged(nPoints+2, 3);
        Uaveraged.setZero();
        gsMatrix<T> parpoints(1, nPoints+2);
        parpoints.setZero();
        gsVector<T> flowRate(nPoints);
        flowRate.setZero();
        T flowRateTotal = 0;
        gsVector<T> flowRateAveraged(nPoints);
        flowRateAveraged.setZero();
        T flowRateAveragedTotal = 0;
        gsVector<T> flowRateinSlice(nPoints);
        flowRateinSlice.setZero();
        T flowRateinSliceTotal = 0;

        gsField<T> rotatedField = constructRotatedField();

        // AREA AVERAGING -------------------------------------------------------------------------------------------
        gsTensorBSplineBasis<2, T> basis = m_rotBasis;
        gsKnotVector<T> kv1 = basis.knots(0);
        gsKnotVector<T> kv2 = basis.knots(1);
        int maxUniqueKnots = kv2.breaks().size();

        gsVector<int> numQuadNodes(2);
        for (int i = 0; i < 2; ++i)
            numQuadNodes[i] = (2 * basis.degree(i) + 1);
        gsGaussRule<T> QuRule(numQuadNodes);// harmless slicing occurs here
        gsMatrix<T> quNodes; // Mapped nodes
        gsVector<T> quWeights; // Mapped weights
        unsigned evFlags = NEED_VALUE | NEED_MEASURE | NEED_OUTER_NORMAL;   // Evaluation flags for the Geometry map

        // Initialize geometry evaluator
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_rotSlice));

        T rr = 0;
        T areaTotal = 0.0;
        gsMatrix<T> stripeBoundaries(2, 2); stripeBoundaries.setZero();
        for (index_t i = 0; i < nPoints; i++)
        {
            // one circumferential area
            gsVector<T> unormal_transf(3); unormal_transf << 0, 0, 0;
            T div_low = i * (kv1.last()-kv1.first()) / (nPoints);
            T div_up = (i+1) * (kv1.last()-kv1.first()) / (nPoints);
            stripeBoundaries << div_low, div_up, 0, 0;
            gsMatrix<T> stripeBoundariesVals = rotatedField.point(stripeBoundaries, 0);
            gsMatrix<T> stripeBoundariesValsTransf(stripeBoundariesVals.rows(),stripeBoundariesVals.cols());
            for (index_t j = 0; j < maxUniqueKnots-1; j++)
            {
                gsVector<T> lowercorner(2);
                lowercorner << div_low, kv2.uValue(j);
                gsVector<T> uppercorner(2);
                uppercorner << div_up, kv2.uValue(j+1);
                QuRule.mapTo(lowercorner, uppercorner, quNodes, quWeights);

                // Compute image of Gauss nodes under geometry mapping as well as Jacobians
                geoEval->evaluateAt(quNodes);

                // Evaluate solution on element nodes
                gsMatrix<T> geometryVals = rotatedField.point(quNodes, 0);
                gsMatrix<T> solutionVals = rotatedField.value(quNodes, 0);

                gsVector<T> bvalsproject(2);
                bvalsproject << stripeBoundariesVals.col(0).row(1), stripeBoundariesVals.col(0).row(2);
                T phib = math::acos(bvalsproject(0) / bvalsproject.norm());
                if (bvalsproject(1) < 0)
                    phib = 2 * EIGEN_PI - phib;
                gsMatrix<real_t> transformMatrix2(3, 3);
                const real_t cos = math::cos(phib);
                const real_t sin = math::sin(phib);
                transformMatrix2(0, 0) = 1;
                transformMatrix2(0, 1) = 0;
                transformMatrix2(0, 2) = 0;
                transformMatrix2(1, 0) = 0;
                transformMatrix2(1, 1) = cos;
                transformMatrix2(1, 2) = sin;
                transformMatrix2(2, 0) = 0;
                transformMatrix2(2, 1) = -sin;
                transformMatrix2(2, 2) = cos;
                stripeBoundariesValsTransf.col(0) = transformMatrix2 * stripeBoundariesVals.col(0);
                stripeBoundariesValsTransf.col(1) = transformMatrix2 * stripeBoundariesVals.col(1);

                for (index_t k = 0; k < quWeights.rows(); k++)
                {
                    const T weight = quWeights(k) * geoEval->measure(k);

                    gsVector<T> yzproject(2);
                    yzproject << geometryVals.col(k).row(1), geometryVals.col(k).row(2);
                    T phi = math::acos(yzproject(0) / yzproject.norm());
                    if (yzproject(1) < 0)
                        phi = 2 * EIGEN_PI - phi;
                    gsMatrix<real_t> transformMatrix(3, 3);
                    const real_t cos = math::cos(phi);
                    const real_t sin = math::sin(phi);
                    transformMatrix(0, 0) = 1;
                    transformMatrix(0, 1) = 0;
                    transformMatrix(0, 2) = 0;
                    transformMatrix(1, 0) = 0;
                    transformMatrix(1, 1) = cos;
                    transformMatrix(1, 2) = sin;
                    transformMatrix(2, 0) = 0;
                    transformMatrix(2, 1) = -sin;
                    transformMatrix(2, 2) = cos;
                    gsVector<T> uvec = transformMatrix * solutionVals.col(k);

                    // Strip area computation
                    area += weight;

                    // Velocity integral on the strip
                    Uaux.col(0) += weight * uvec.row(0);
                    Uaux.col(1) += weight * uvec.row(1);
                    Uaux.col(2) += weight * uvec.row(2);

                    // Compute the outer normal vector and flowrate through the strip
                    gsVector<T> unormal, unormal2;
                    geoEval->normal(k, unormal);
                    unormal = unormal/unormal.norm();
                    unormal2 = transformMatrix * unormal;
                    unormal_transf += unormal_transf + unormal2;
                    unormal_transf = unormal_transf/unormal_transf.norm();

                    flowRate(i) += weight * unormal.dot(solutionVals.col(k));

                }

            }

            //Uaveraged.row(i+1) = Uaux.transpose() / area;
            Uaveraged.row(i+1) = Uaux / area;
            //Uaveraged.row(i+1) = Uaveraged.row(i+1)/Uaveraged.row(i+1).norm();

            parpoints(i+1) = (div_low + div_up) / 2;

            areaTotal += area;
            flowRateTotal += flowRate(i);

            T r1 = sqrt(pow(stripeBoundariesValsTransf(1,0), 2) + pow(stripeBoundariesValsTransf(2,0), 2));
            T r2 = sqrt(pow(stripeBoundariesValsTransf(1,1), 2) + pow(stripeBoundariesValsTransf(2,1), 2));
            flowRateinSlice(i) = flowRate(i) / (EIGEN_PI * (r1 + r2));
            flowRateinSliceTotal += flowRateinSlice(i);

            area = 0;
            rr = 0;
            Uaux.setZero();

        }
        parpoints(nPoints+1) = kv1.last();
        gsMatrix<T> boundaryParPoints(2,2);
        boundaryParPoints << kv1.first(), kv1.last(),
                             kv2.first(), kv2.first();
        gsMatrix<T> boundaryPoints = rotatedField.point(boundaryParPoints, 0);
        gsVector<T> vec = boundaryPoints.col(1) - boundaryPoints.col(0);
        T dist = vec.norm();

        gsMatrix<T> parpointssurf(2, parpoints.cols()); parpointssurf.setZero();
        //parpointssurf << parpoints;
        parpointssurf.row(0) = parpoints;
        geoEval->evaluateAt(parpointssurf);
        gsMatrix<T> geometryValsParPoints = rotatedField.point(parpointssurf, 0);

        for (index_t i = 0; i < parpoints.cols(); i++) {
            gsVector<T> yzproject(2);
            yzproject << geometryValsParPoints.col(i).row(1), geometryValsParPoints.col(i).row(2);
            T phi = math::acos(yzproject(0) / yzproject.norm());
            if (yzproject(1) < 0)
                phi = 2 * EIGEN_PI - phi;
            gsMatrix<real_t> transformMatrix(3, 3);
            const real_t cos = math::cos(phi);
            const real_t sin = math::sin(phi);
            transformMatrix(0, 0) = 1;
            transformMatrix(0, 1) = 0;
            transformMatrix(0, 2) = 0;
            transformMatrix(1, 0) = 0;
            transformMatrix(1, 1) = cos;
            transformMatrix(1, 2) = -sin;
            transformMatrix(2, 0) = 0;
            transformMatrix(2, 1) = sin;
            transformMatrix(2, 2) = cos;
            Uaveraged.row(i) = (transformMatrix * Uaveraged.row(i).transpose()).transpose();
        }

        // Data approximation via B-spline curve
        gsMatrix<T> dersStart(1,3);
        gsMatrix<T> dersEnd(1,3);
        dersStart.row(0) = Uaveraged.row(0);
        dersEnd.row(0) = Uaveraged.row(Uaveraged.rows()-1);
        //gsBSpline<T> UaveragedCurve = curveFittingWithEndDerivatives(Uaveraged, parpoints, kv1, dersStart, dersEnd);
        //gsBSpline<T> UaveragedCurve = curveFittingWithBoundary(Uaveraged, parpoints, kv1, true);
        //gsBSpline<T> UaveragedCurve = curveFittingWithBoundaryPreservingFlowrate(Uaveraged, parpoints, kv1, m_rotSlice, m_rotBasis, dist*flowRateTotal/areaTotal, true);
        gsBSpline<T> UaveragedCurve = curveFittingWithBoundaryPreservingFlowrate(Uaveraged, parpoints, kv1, m_rotSlice, m_rotBasis, flowRateinSlice, true);
        gsMatrix<T> UaveragedCurveCoefsNew = UaveragedCurve.coefs();

        // B-spline surface representing averaged velocity profile
        gsMatrix<real_t> transformMatrix(3, 3);
        real_t cos;
        real_t sin;
        gsKnotVector<T> kv2part = m_origBasis.knots(0);
        gsKnotVector<T> kv2partElevated = kv2part;
        gsKnotVector<T> kvpom(kv2part.first(), kv2part.last(), 0, 4);
        std::vector<T> kvdiff;
        kv2part.difference(kvpom, kvdiff);
        gsBSpline<T> bscurve;
        gsMatrix<T> axis(3, 1);
        axis << 1, 0, 0;
        gsMatrix<T> centre(3, 1);
        gsMatrix<T> start(3, 1);
        gsMatrix<T> result(4, 3);
        gsMatrix<T> result2(kv2part.size()-kv2part.degree()-1, 3);
        gsMatrix<T> Uxcoefs(UaveragedCurve.coefsSize(), kv2partElevated.size()-kv2partElevated.degree()-1);
        gsMatrix<T> Uycoefs(UaveragedCurve.coefsSize(), kv2partElevated.size()-kv2partElevated.degree()-1);
        gsMatrix<T> Uzcoefs(UaveragedCurve.coefsSize(), kv2partElevated.size()-kv2partElevated.degree()-1);
        for (index_t i = 0; i < UaveragedCurve.coefsSize(); i++)
        {
            centre << UaveragedCurve.coef(i, 0), 0, 0;
            start = UaveragedCurve.coef(i).transpose();
            if (start.norm() == 0)
                result.setZero();
            else
                result = computeCircleArc3D(start, centre, -getPhi(), axis);
            bscurve = gsBSpline<T> (kvpom, result);
            if (kv2part.degree() > kvpom.degree())
            {
                int div_degree = kv2part.degree() - kvpom.degree();
                bscurve.degreeElevate(div_degree);
            }
            for (index_t j = 0; j < kvdiff.size(); j++)
            {
                if ((kvdiff[j] > kv2part.first()) && (kvdiff[j] < kv2part.last()))
                    bscurve.insertKnot(kvdiff[j]);
            }
            result2 = bscurve.coefs();

            Uxcoefs.row(i) = result2.col(0).transpose();
            Uycoefs.row(i) = result2.col(1).transpose();
            Uzcoefs.row(i) = result2.col(2).transpose();
        }

        gsMatrix<T> UaveragedSurfaceCoefs(Uxcoefs.rows() * Uxcoefs.cols(), 3);
        UaveragedSurfaceCoefs.setZero();
        for (index_t i = 0; i < Uxcoefs.cols(); i++)
        {
            UaveragedSurfaceCoefs.block(i * Uxcoefs.rows(),0,Uxcoefs.rows(),1) = Uxcoefs.col(i);
            UaveragedSurfaceCoefs.block(i * Uxcoefs.rows(),1,Uxcoefs.rows(),1) = Uycoefs.col(i);
            UaveragedSurfaceCoefs.block(i * Uxcoefs.rows(),2,Uxcoefs.rows(),1) = Uzcoefs.col(i);
        }
        m_origAveragedCoefs = UaveragedSurfaceCoefs;

        // Rotation to include other periodic parts
        cos = math::cos(getPhi());
        sin = math::sin(getPhi());
        transformMatrix(0, 0) = 1;
        transformMatrix(0, 1) = 0;
        transformMatrix(0, 2) = 0;
        transformMatrix(1, 0) = 0;
        transformMatrix(1, 1) = cos;
        transformMatrix(1, 2) = sin;
        transformMatrix(2, 0) = 0;
        transformMatrix(2, 1) = -sin;
        transformMatrix(2, 2) = cos;

        gsMatrix<T> UaveragedSurfaceRotatedCoefs(UaveragedSurfaceCoefs.rows() * m_numBlades, 3);
        UaveragedSurfaceRotatedCoefs.middleRows(0, getNumSliceCoefs()) = UaveragedSurfaceCoefs;

        for (int i = 1; i < m_numBlades; i++)
        {
            for (int j = 0; j < getNumSliceCoefs(); j++)
            {
                UaveragedSurfaceRotatedCoefs.row(i*getNumSliceCoefs() + j) = transformMatrix * UaveragedSurfaceRotatedCoefs.row((i - 1)*getNumSliceCoefs() + j).transpose();
            }
        }
        m_rotAveragedCoefs = UaveragedSurfaceRotatedCoefs;

        // Total original flowrate through the area
        flowRateTotal = 0;
        for (index_t i = 0; i < nPoints; i++)
        {
            flowRateTotal += flowRate(i);
        }
        gsInfo << "Total original flowrate through the area is: " << flowRateTotal << " m3/s\n";

        // Computing flow rate after averaging
        gsField<T> rotatedAveragedField = constructRotatedAveragedField();
        for (index_t i = 0; i < nPoints; i++)
        {
            // one circumferential area
            T div_low = i * (kv1.last()-kv1.first()) / (nPoints);
            T div_up = (i+1) * (kv1.last()-kv1.first()) / (nPoints);
            for (index_t j = 0; j < maxUniqueKnots-1; j++)
            {
                gsVector<T> lowercorner(2);
                lowercorner << div_low, kv2.uValue(j);
                gsVector<T> uppercorner(2);
                uppercorner << div_up, kv2.uValue(j+1);
                QuRule.mapTo(lowercorner, uppercorner, quNodes, quWeights);

                // Compute image of Gauss nodes under geometry mapping as well as Jacobians
                geoEval->evaluateAt(quNodes);

                // Evaluate solution on element nodes
                gsMatrix<T> solutionVals = rotatedAveragedField.value(quNodes, 0);

                for (index_t k = 0; k < quWeights.rows(); k++)
                {
                    const T weight = quWeights(k) * geoEval->measure(k);

                    // Compute the outer normal vector and flowrate through the strip
                    gsVector<T> unormal;
                    geoEval->normal(k, unormal);
                    unormal = unormal/unormal.norm();

                    flowRateAveraged(i) += weight * unormal.dot(solutionVals.col(k));

                }

            }

        }

        // Total flowrate through the area after averaging
        flowRateAveragedTotal = 0;
        for (index_t i = 0; i < nPoints; i++)
        {
            flowRateAveragedTotal += flowRateAveraged(i);
        }
        gsInfo << "Total flowrate through the area after averaging is: " << flowRateAveragedTotal << " m3/s\n";

    }

public:
    void compute()
    {
        if (!m_compute) {
            rotatePatch();
            rotateBasis();
            rotateCoefs();
            m_compute = true;
        }
    }

    void computeWithAveraging(int nPoints = 30)
    {
        compute();
        areaAverage(nPoints);
    }

    gsTensorBSpline<2, T> getSlice()
    {
        return m_origSlice;
    }

    gsTensorBSpline<2, T> getRotatedSlice()
    {
        return m_rotSlice;
    }

    gsMatrix<T> getSliceCoefs()
    {
        return m_origCoefs;
    }

    gsMatrix<T> getRotatedCoefs()
    {
        return m_rotCoefs;
    }

    gsField<T> constructRotatedField()
    {
        m_rotatedSlice.clear();
        m_rotatedField.clear();

        m_rotatedSlice.addPatch(m_rotSlice);
        m_rotatedField.addPatch(m_rotBasis.makeGeometry(m_rotCoefs));
        
        return gsField<T>(m_rotatedSlice, m_rotatedField, true);
    }

    gsField<T> constructRotatedAveragedField()
    {
        m_rotatedSlice.clear();
        m_rotatedAveragedField.clear();

        m_rotatedSlice.addPatch(m_rotSlice);
        m_rotatedAveragedField.addPatch(m_rotBasis.makeGeometry(m_rotAveragedCoefs));

        return gsField<T>(m_rotatedSlice, m_rotatedAveragedField, true);
    }

    gsField<T> constructPartialField()
    {
        m_originalSlice.clear();
        m_originalField.clear();

        m_originalSlice.addPatch(m_origSlice);
        m_originalField.addPatch(m_origBasis.makeGeometry(m_origCoefs));

        return gsField<T>(m_originalSlice, m_originalField, true);
    }

    gsField<T> constructPartialAveragedField()
    {
        m_originalSlice.clear();
        m_originalAveragedField.clear();

        m_originalSlice.addPatch(m_origSlice);
        m_originalAveragedField.addPatch(m_origBasis.makeGeometry(m_origAveragedCoefs));

        return gsField<T>(m_originalSlice, m_originalAveragedField, true);
    }


    size_t getNumSliceCoefs() { return m_origCoefs.rows(); }

    real_t getPhi()
    {
        real_t phi = 2 * EIGEN_PI / m_numBlades;
        return phi;
    }

protected:
    index_t m_knot, m_numBlades, m_dir;
    bool m_swap, m_compute;
    gsTensorBSpline<2, T> m_origSlice, m_rotSlice;
    gsTensorBSplineBasis<2, T> m_origBasis, m_rotBasis;
    gsMatrix<T> m_origCoefs, m_rotCoefs, m_origAveragedCoefs, m_rotAveragedCoefs;
    gsMultiPatch<T> m_originalSlice, m_rotatedSlice, m_originalField, m_rotatedField, m_originalAveragedField, m_rotatedAveragedField;

}; // class uwbSliceFieldCreator

} // namespace gismo

