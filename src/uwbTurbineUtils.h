/** @file uwbTurbineUtils.h

Author(s): B. Bastl, K. Michalkova
*/

#ifndef UWBTURBINEUTILS
#define UWBTURBINEUTILS

#pragma once

#include <gismo.h>

#include "gsCore/gsGeometryEvaluator.h"

using namespace gismo;

template<class T>
gsMatrix<T> centripetalParameterization(gsMatrix<T> &points);

template<class T>
gsMatrix<T> chordalParameterization(gsMatrix<T> &points);

template<class T>
gsBSpline<T> curveFittingWithBoundary(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kvfit, bool print_info);

template<class T>
gsBSpline<T> curveFittingWithBoundaryPreservingFlowrate(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kvfit, gsTensorBSpline<2, T> m_rotSlice, gsTensorBSplineBasis<2, T> m_rotBasis, gsVector<T> FRconst, bool print_info);

template<class T>
gsBSpline<T> curveFittingWithBoundaryAndInputTangent(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kvfit, gsMatrix<T> T0);

template<class T>
gsBSpline<T> curveFittingWithEndDerivatives(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kvfit, gsMatrix<T> dersStart, gsMatrix<T> dersEnd);

template<class T>
int curveSphereIntersectionviaNR(gsTensorBSpline<2, T> & surface1, T par_surf1, T initial_solution, T radius, T & solution, T tol, int maxiter, int & num_iter, bool print_info);

template<class T>
int curveSurfaceIntersectionviaNR(gsTensorBSpline<2, T> & surface1, T par_surf1, gsTensorBSpline<2, T> & surface2, gsVector<T> initial_solution, gsVector<T> & solution, T tol, int maxiter, int & num_iter, bool print_info);

template<class T>
int computeLoftSurface(gsBSpline<T> defcurves[], gsKnotVector<T> kvcurves, int num_bladeprofiles, gsKnotVector<T> kvloft, gsMatrix<T> sec_parameters, gsTensorBSpline<2,T> & finalsurface);

template<class T>
int computeLoftSurface(std::vector<gsBSpline<T> > defcurves, gsKnotVector<T> kvcurves, int num_bladeprofiles, gsKnotVector<T> kvloft, gsMatrix<T> sec_parameters, gsTensorBSpline<2,T> & finalsurface);

template<class T>
int computeLoftFuction(std::vector<gsBSpline<T> > defcurves, gsKnotVector<T> kvcurves, int num_bladeprofiles, gsKnotVector<T> kvloft, gsMatrix<T> sec_parameters,gsTensorBSpline<2,T> & finalsurface, int dimension);

template<class T>
int fitData(gsVector<T> datax, gsVector<T> datay, int degree, gsMatrix<T> & polynomialcoefficients);

template<class T>
int fitError(gsVector<T> datax, gsVector<T> datay, gsMatrix<T> polynomialcoefficients, T & error);

template<class T>
int extrapolateData(gsVector<T> datax, gsVector<T> datay, T r, T & extra_data, T tol, int max_iter);

template<class T>
int extrapolateData(gsVector<T> datax, gsVector<T> datay, T r, T & extra_data, T tol, int max_iter);

template<class T>
gsMatrix<T> fitLeastSquares(gsMatrix<T> data);

template<class T>
int constructCylinderParametrization(T p0, T p1, T radius, gsTensorBSpline<2, T> & cylinder);

template<class T>
int minimizePointSurfaceDistanceviaNR(gsVector<T> point, gsTensorBSpline<2, T> & surface, gsVector<T> initial_solution, gsVector<T> & solution, T tol, int maxiter, int & num_iter, bool print_info);

template<class T>
int cross(gsVector<T> vec1, gsVector<T> vec2, gsVector<T> &result);

template<class T>
T euclideanNorm(gsMatrix<T> v);

template<class T>
T euclideanNorm(gsVector<T> v);

template<class T>
int discreteCoonsPatch(gsMatrix<T> a_cp, gsMatrix<T> b_cp, gsMatrix<T> c_cp, gsMatrix<T> d_cp, gsMatrix<T> & final_control_net, bool print_info);

template<class T>
int discreteCoonsPatch(gsMatrix<T> a_cp, gsMatrix<T> b_cp, gsMatrix<T> c_cp, gsMatrix<T> d_cp, T w1, T w2, T w3, gsMatrix<T> & final_control_net, bool print_info);

template<class T>
int springModelPatch(gsMatrix<T> a_cp, gsMatrix<T> b_cp, gsMatrix<T> c_cp, gsMatrix<T> d_cp, gsMatrix<T> & final_control_net, bool print_info);

template<class T>
int projectToCylinder(gsMatrix<T> point, T radius, gsMatrix<T> & new_point);

template<class T>
int discreteCoonsPatchModified(gsMatrix<T> a_cp, gsMatrix<T> b_cp, gsMatrix<T> c_cp, gsMatrix<T> d_cp, gsMatrix<T> & final_control_net);

template<class T>
gsTensorBSpline<2,T> surfaceFittingLSQWithBoundary(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kv1, gsKnotVector<T> kv2, gsMatrix<T> boundary_points, bool print_info);

template<class T>
int discreteCoonsPatch3D(gsTensorBSpline<2,T> a, gsTensorBSpline<2,T> b, gsTensorBSpline<2,T> c, gsTensorBSpline<2,T> d, gsTensorBSpline<2,T> e, gsTensorBSpline<2,T> f, gsTensorBSpline<3,T> & final_volume, bool print_info);

template<class T>
int springModelPatch3D(gsTensorBSpline<2,T> a, gsTensorBSpline<2,T> b, gsTensorBSpline<2,T> c, gsTensorBSpline<2,T> d, gsTensorBSpline<2,T> e, gsTensorBSpline<2,T> f, gsTensorBSpline<3,T> & final_volume, bool print_info);

template<class T>
gsMatrix<T> computeCircleArc2D(gsMatrix<T> start, T centreX, T centreY, T alpha);

template<class T>
gsMatrix<T> computeCircleArc3D(gsMatrix<T> start, gsMatrix<T> centre, T alpha, gsMatrix<T> axis);

template<class T>
gsMultiPatch<T> makeLinear(gsMultiPatch<T> mp,  int uR);

template<class T>
gsMatrix<T> LSQFergusonSmooth(gsMatrix<T> P0, gsMatrix<T> T0, gsMatrix<T> T1, gsMatrix<T> P3);

template<class T>
gsMatrix<T> AxisDirectionNormed(gsMatrix<T> T0, gsMatrix<T> T1);

// Change the knot vector of curve1 to be the same as the knot vector of curve2
template<class T>
void UnifyKnotVectors(gsBSpline<T> &curve1, gsBSpline<T> curve2);

template<class T>
gsMatrix<T> LSQFergusonShort(gsMatrix<T> P0, gsMatrix<T> T0, gsMatrix<T> T1, gsMatrix<T> P3);

template<class T>
std::vector<T> smartKnotIdentification(std::vector<T> kv_unique);

template<class T>
std::vector<T> smartKnotIdentification2(std::vector<T> kv_unique);

template<class T>
T smartKnotParameterInsert(T param, gsKnotVector<T> kv);




// Fits the given points witha B-spline curve where the first and the last control points correspond
// to the first and the last points of the given points, respectively.
template<class T>
gsBSpline<T> curveFittingWithBoundary(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kvfit, bool print_info)
{
    if (print_info) {
        gsInfo << "\nCurve approximation with least squares:\n";
        gsInfo << "=======================================\n";
    }
    // number of points
    int num_points=points.rows();
    // dimension of points will be also later dimension of coefficients
    int m_dimension=points.cols();
    // basis function definition
    //gsBSplineBasis<T> *curveBasis = new gsBSplineBasis<T>(kvfit);
    gsBSplineBasis<T> curveBasis = gsBSplineBasis<T>(kvfit);
    //number of basis functions
    index_t num_basis=curveBasis.size();


    //for computing the value of the basis function
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    //computing the values of the basis functions at some position
    curveBasis.eval_into(parameter_points,values);
    // which functions have been computed i.e. which are active
    curveBasis.active_into(parameter_points,actives);

    //how many rows and columns has the A matrix and how many rows has the b vector
    int num_rows=num_basis;
    //left side matrix
    gsMatrix<T> m_A(num_rows-2,num_rows-2);
    m_A.setZero(num_rows-2,num_rows-2); // ensure that all entries are zero in the beginning
    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_rows-2,m_dimension);
    m_B.setZero(num_rows-2,m_dimension); // enusure that all entris are zero in the beginning

    gsMatrix<T> R0values(1, num_points);
    gsBasisFun<T> R0 = curveBasis.function(0);
    R0.eval_into(parameter_points, R0values);
    gsMatrix<T> RNvalues(1, num_points);
    gsBasisFun<T> RN = curveBasis.function(num_basis-1);
    RN.eval_into(parameter_points, RNvalues);

    // building the matrix A and the vector b of the system of linear equations A*x==b(uses automatically the information of closed or not)
    for(index_t k=0;k<num_points;k++){
        for(index_t i=0;i<actives.rows();i++){
            if ((actives(i,k) > 0) && (actives(i,k) < num_basis-1)) {
                m_B.row(actives(i,k)-1) += values(i,k)*(points.row(k) - points.row(0)*R0values(k) - points.row(num_points-1)*RNvalues(k));
            }
            for(index_t j=0;j<actives.rows();j++){
                if ((actives(i,k) > 0) && (actives(j,k) > 0) && (actives(i,k) < num_basis-1) && (actives(j,k) < num_basis-1)) {
                    m_A(actives(i,k)-1,actives(j,k)-1) += values(i,k)*values(j,k);
                }

            }
        }
    }

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)
    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x=m_A.fullPivHouseholderQr().solve( m_B);

    // finally generate the B-spline curve
    gsMatrix<T> coefs(num_rows, m_dimension);
    coefs.row(0) = points.row(0);
    for (index_t i=0; i < num_rows-2; i++) {
        coefs.row(i+1) = x.row(i);
    }
    coefs.row(num_rows-1) = points.row(num_points-1);

    gsBSpline<T> curve = gsBSpline<T> (curveBasis, give(coefs));
    if (print_info) {
        gsInfo << "Final curve:\n" << curve << "\n";
    }

    return curve;
}

// Fits the given points witha B-spline curve where the first and the last control points correspond
// to the first and the last points of the given points, respectively.
template<class T>
gsBSpline<T> curveFittingWithBoundaryPreservingFlowrate(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kvfit, gsTensorBSpline<2, T> m_rotSlice, gsTensorBSplineBasis<2, T> m_rotBasis, gsVector<T> FRconst, bool print_info)
{
    if (print_info) {
        gsInfo << "\nCurve approximation with least squares:\n";
        gsInfo << "=======================================\n";
    }
    // number of points
    int num_points=points.rows();
    // dimension of points will be also later dimension of coefficients
    int m_dimension=points.cols();
    // basis function definition
    gsBSplineBasis<T> curveBasis = gsBSplineBasis<T>(kvfit);
    //number of basis functions
    unsigned int num_basis=curveBasis.size();

    index_t nPoints = FRconst.rows();
    gsVector<T> FRconst2 = FRconst;

    //for computing the value of the basis function
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    //computing the values of the basis functions at some position
    curveBasis.eval_into(parameter_points,values);
    // which functions have been computed i.e. which are active
    curveBasis.active_into(parameter_points,actives);

    int num_rows=num_basis; //how many rows and columns has the A matrix and how many rows has the b vector
    gsMatrix<T> m_NN(num_rows-2,num_rows-2); //left side matrix
    m_NN.setZero(); // ensure that all entries are zero in the beginning
    gsMatrix<T> m_Mx(num_rows-2, nPoints); //right side vector (more dimensional!)
    gsMatrix<T> m_My(num_rows-2, nPoints); //right side vector (more dimensional!)
    gsMatrix<T> m_Mz(num_rows-2, nPoints); //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_rows-2,m_dimension); //right side vector (more dimensional!)
    m_B.setZero(); // enusure that all entris are zero in the beginning
    m_Mx.setZero();
    m_My.setZero();
    m_Mz.setZero();

    gsMatrix<T> R0values(1, num_points);
    gsBasisFun<T> R0 = curveBasis.function(0);
    R0.eval_into(parameter_points, R0values);
    gsMatrix<T> RNvalues(1, num_points);
    gsBasisFun<T> RN = curveBasis.function(num_basis-1);
    RN.eval_into(parameter_points, RNvalues);

    // building the matrix A and the vector b of the system of linear equations A*x==b(uses automatically the information of closed or not)
    for(index_t k=0;k<num_points;k++){
        for(index_t i=0;i<actives.rows();i++){
            if ((actives(i,k) > 0) && (actives(i,k) < num_basis-1)) {
                m_B.row(actives(i,k)-1) += values(i,k)*(points.row(k) - points.row(0)*R0values(k) - points.row(num_points-1)*RNvalues(k));
            }
            for(index_t j=0;j<actives.rows();j++){
                if ((actives(i,k) > 0) && (actives(j,k) > 0) && (actives(i,k) < num_basis-1) && (actives(j,k) < num_basis-1)) {
                    m_NN(actives(i,k)-1,actives(j,k)-1) += values(i,k)*values(j,k);
                }

            }
        }
    }

    //gsTensorBSplineBasis<2, T> basis = m_rotSlice.basis();
    gsBSplineBasis<T> basis = m_rotBasis.component(0);
    gsKnotVector<T> kv1 = basis.knots(0);

    index_t numQuadNodes;
    numQuadNodes = (2 * basis.degree() + 1);
    gsGaussRule<T> QuRule(numQuadNodes);// harmless slicing occurs here
    gsMatrix<T> quNodes; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights
    unsigned evFlags = NEED_VALUE | NEED_MEASURE | NEED_OUTER_NORMAL;   // Evaluation flags for the Geometry map

    // Initialize geometry evaluator
    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_rotSlice));

    gsBSpline<T> rotSliceBeginning;
    m_rotSlice.slice(1, 0, rotSliceBeginning);
    typename gsGeometryEvaluator<T>::uPtr geoEvalSlice(getEvaluator(evFlags, rotSliceBeginning));

    T rr = 0;
    T length = 0;
    T lengthTotal = 0;
    gsMatrix<T> stripeBoundaries(2, 2); stripeBoundaries.setZero();
    for (index_t j = 0; j < nPoints; j++)
    {
        T div_low = j * (kv1.last()-kv1.first()) / (nPoints);
        T div_up = (j+1) * (kv1.last()-kv1.first()) / (nPoints);
        QuRule.mapTo(div_low, div_up, quNodes, quWeights);
        gsMatrix<T> quNodes2(2, quNodes.cols()); quNodes2.setZero();
        quNodes2.row(0) = quNodes.row(0);
        stripeBoundaries << div_low, div_up, 0, 0;
        gsMatrix<T> stripeBoundariesVals = m_rotSlice.eval(stripeBoundaries);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(quNodes2);
        geoEvalSlice->evaluateAt(quNodes);

        // Evaluate solution on element nodes
        gsMatrix<T> geometryVals = m_rotSlice.eval(quNodes2);
        gsMatrix<T> basisVals = basis.eval(quNodes);
        gsMatrix<index_t> activesVals = basis.active(quNodes);

        for (index_t k = 0; k < quWeights.rows(); k++)
        {
            //const T weight = quWeights(k) * geoEval->measure(k);
            const T weight = quWeights(k) * geoEvalSlice->measure(k);

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

            // Strip area computation
            length += weight;

            // Compute the outer normal vector and flowrate through the strip
            gsVector<T> unormal;
            geoEval->normal(k, unormal);
            unormal = unormal/unormal.norm();

            for (index_t i = 0; i < activesVals.rows(); i++) {
                if ((activesVals(i,k) > 0) && (activesVals(i,k) < num_basis-1)) {
                    m_Mx(activesVals(i,k)-1,j) += weight * basisVals(i,k) * unormal(0);
                    m_My(activesVals(i,k)-1,j) += weight * basisVals(i,k) * unormal(1);
                    m_Mz(activesVals(i,k)-1,j) += weight * basisVals(i,k) * unormal(2);
                }
            }

        }

        lengthTotal += length;
        length = 0;
        rr = 0;

    }

    gsMatrix<T> m_Zero(m_NN.rows(), m_NN.rows());
    m_Zero.setZero(); // ensure that all entries are zero in the beginning
    gsMatrix<T> m_Zero2(m_Mx.cols(), m_Mx.cols());
    m_Zero2.setZero(); // ensure that all entries are zero in the beginning
    gsMatrix<T> m_A(3 * m_NN.rows() + nPoints, 3 * m_NN.rows() + nPoints);
    m_A.setZero(); // ensure that all entries are zero in the beginning
    m_A << m_NN,                m_Zero,             m_Zero,             0.5 * m_Mx,
           m_Zero,              m_NN,               m_Zero,             0.5 * m_My,
           m_Zero,              m_Zero,             m_NN,               0.5 * m_Mz,
           m_Mx.transpose(),    m_My.transpose(),   m_Mz.transpose(),   m_Zero2;

    // right hand side of the linear system
    gsMatrix<T> m_B2(3 * m_NN.rows() + nPoints, 1);
    m_B2.setZero(); // enusure that all entris are zero in the beginning
    m_B2 << m_B.col(0), m_B.col(1), m_B.col(2), FRconst2;

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)
    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x=m_NN.fullPivHouseholderQr().solve( m_B);

    gsVector<T> xtest(nPoints);
    T xtesttotal = 0;
    for (index_t i = 0; i < nPoints; i++) {
        xtest(i) = m_Mx.col(i).dot(x.col(0)) + m_My.col(i).dot(x.col(1)) + m_Mz.col(i).dot(x.col(2));
        xtesttotal += xtest(i);
    }

    gsMatrix<T> x2 (m_B2.rows(), m_B2.cols());
    x2=m_A.fullPivHouseholderQr().solve( m_B2);
    gsMatrix<T> x3 (m_B.rows(), m_B.cols());
    for (index_t i = 0; i < x3.cols(); i++)
        for (index_t j = 0; j < x3.rows(); j++)
            x3(j, i) = x2(i*m_NN.rows() + j);

    gsVector<T> xtest2(nPoints);
    T xtest2total = 0;
    for (index_t i = 0; i < nPoints; i++) {
        xtest2(i) = m_Mx.col(i).dot(x3.col(0)) + m_My.col(i).dot(x3.col(1)) + m_Mz.col(i).dot(x3.col(2));
        xtest2total += xtest2(i);
    }

    // finally generate the B-spline curve
    gsMatrix<T> coefs(num_rows, m_dimension);
    coefs.row(0) = points.row(0);
    for (index_t i=0; i < num_rows-2; i++) {
        coefs.row(i+1) = x3.row(i);
    }
    coefs.row(num_rows-1) = points.row(num_points-1);

    gsBSpline<T> curve = gsBSpline<T> (curveBasis, give(coefs));
    if (print_info) {
        gsInfo << "Final curve:\n" << curve << "\n";
    }

    // test of the result
    gsVector<T> flowRateinLSQ(nPoints);
    flowRateinLSQ.setZero();
    T flowRateinLSQTotal = 0;
    for (index_t j = 0; j < nPoints; j++)
    {
        T div_low = j * (kv1.last()-kv1.first()) / (nPoints);
        T div_up = (j+1) * (kv1.last()-kv1.first()) / (nPoints);
        QuRule.mapTo(div_low, div_up, quNodes, quWeights);
        gsMatrix<T> quNodes2(2, quNodes.cols()); quNodes2.setZero();
        quNodes2.row(0) = quNodes.row(0);
        stripeBoundaries << div_low, div_up, 0, 0;
        gsMatrix<T> stripeBoundariesVals = m_rotSlice.eval(stripeBoundaries);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(quNodes2);
        geoEvalSlice->evaluateAt(quNodes);

        // Evaluate solution on element nodes
        gsMatrix<T> geometryVals = m_rotSlice.eval(quNodes2);
        gsMatrix<T> solutionVals = curve.eval(quNodes);

        for (index_t k = 0; k < quWeights.rows(); k++)
        {
            //const T weight = quWeights(k) * geoEval->measure(k);
            const T weight = quWeights(k) * geoEvalSlice->measure(k);


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

            // Strip area computation
            length += weight;


            // Compute the outer normal vector and flowrate through the strip
            gsVector<T> unormal;
            geoEval->normal(k, unormal);
            unormal = unormal/unormal.norm();

            flowRateinLSQ(j) += weight * unormal.dot(solutionVals.col(k));

        }
        flowRateinLSQTotal += flowRateinLSQ(j);

    }

    return curve;
}

/* OLD VERSION, NOT WORKING FOR 3D
// Fits the given points witha B-spline curve where the first and the last control points correspond
// to the first and the last points of the given points, respectively.
gsBSpline<T> curveFittingWithBoundaryAndInputTangent(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kvfit, gsMatrix<T> T0)
{
    // number of points
    int num_points=points.rows();
    // dimension of points will be also later dimension of coefficients
    int m_dimension=points.cols();
    // basis function definition
    gsBSplineBasis<T> *curveBasis = new gsBSplineBasis<T>(kvfit);
    //number of basis functions
    unsigned int num_basis=curveBasis->size();


    //for computing the value of the basis function
    gsMatrix<T> values;
    gsMatrix<unsigned> actives;

    //computing the values of the basis functions at some position
    curveBasis->eval_into(parameter_points,values);
    // which functions have been computed i.e. which are active
    curveBasis->active_into(parameter_points,actives);

    //gsInfo << "Vycislene bazove funkce:\n" << values << "\n\n";
    //gsInfo << "Aktivni bazove funkce:\n" << actives << "\n\n";

    //how many rows and columns has the A matrix and how many rows has the b vector
    int num_rows=num_basis;
    //left side matrix
    gsMatrix<T> m_A0(num_rows-2,num_rows-2);
    m_A0.setZero(); // ensure that all entries are zero in the beginning
    //right side vector (more dimensional!)
    gsMatrix<T> m_B0(num_rows-2,m_dimension);
    m_B0.setZero(); // enusure that all entris are zero in the beginning

    gsMatrix<T> R0values(1, num_points);
    gsBasisFun<T> R0 = curveBasis->function(0);
    R0.eval_into(parameter_points, R0values);
    //gsInfo << "R0:\n" << R0values << "\n";

    gsMatrix<T> R1values(1, num_points);
    gsBasisFun<T> R1 = curveBasis->function(1);
    R1.eval_into(parameter_points, R1values);
    //gsInfo << "R1:\n" << R1values << "\n";

    gsMatrix<T> RNvalues(1, num_points);
    gsBasisFun<T> RN = curveBasis->function(num_basis-1);
    RN.eval_into(parameter_points, RNvalues);
    //gsInfo << "RN:\n" << RNvalues << "\n";

    // building auxiliary matrix A and vector b for consequent building of the final matrix A and the vector b of the system of linear equations A*x==b
    real_t m_Bp = 0;
    for(index_t k=0;k<num_points;k++){
        m_Bp += R1values(k)*(T0(0)*points(k,0) + T0(1)*points(k,1) - (points(0,0)*T0(0) + points(0,1)*T0(1))*(R0values(k)+R1values(k)) - RNvalues(k)*(T0(0)*points(num_basis-1,0) + T0(1)*points(num_basis-1,1)));
        for(index_t i=0;i<actives.rows();i++){
            if ((actives(i,k) > 0) && (actives(i,k) < num_basis-1)) {
                //- points.row(0)*R1values(k)
                m_B0.row(actives(i,k)-1) += values(i,k)*(points.row(k) - points.row(0)*R0values(k) - points.row(num_points-1)*RNvalues(k));
            }
            for(index_t j=0;j<actives.rows();j++){
                if ((actives(i,k) > 0) && (actives(j,k) > 0) && (actives(i,k) < num_basis-1) && (actives(j,k) < num_basis-1)) {
                    m_A0(actives(i,k)-1,actives(j,k)-1) += values(i,k)*values(j,k);
                }

            }
        }
    }

    //gsInfo << T0 << "\n";
    //gsInfo << m_A0 << "\n";
    //gsInfo << m_B0 << "\n";
    //gsInfo << m_B0.block(1,0,num_rows-3,1) << "\n";

    //left side matrix of the linear system
    gsMatrix<T> m_Zero(num_rows-3,num_rows-3);
    m_Zero.setZero(); // ensure that all entries are zero in the beginning
    gsMatrix<T> m_A(2*(num_rows-3)+1,2*(num_rows-3)+1);
    m_A.setZero(); // ensure that all entries are zero in the beginning
    m_A << (T0(0)*T0(0) + T0(1)*T0(1))*m_A0(0,0),       T0(0)*m_A0.block(0,1,1,num_rows-3),         T0(1)*m_A0.block(0,1,1,num_rows-3),
           T0(0)*m_A0.block(1,0,num_rows-3,1),          m_A0.block(1,1,num_rows-3,num_rows-3),      m_Zero,
           T0(1)*m_A0.block(1,0,num_rows-3,1),          m_Zero,                                     m_A0.block(1,1,num_rows-3,num_rows-3);

    // right hand side of the linear system
    gsMatrix<T> m_B(2*(num_rows-3)+1,1);
    m_B.setZero(); // enusure that all entris are zero in the beginning
    m_B << m_Bp, m_B0.block(1,0,num_rows-3,1), m_B0.block(1,1,num_rows-3,1);

    gsInfo << "Soustava sestavena. \n";

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)
    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x = m_A.fullPivHouseholderQr().solve( m_B);

    gsInfo << "Soustava vyresena.\n" << x << "\n" << x.size() << "\n";

    // finally generate the B-spline curve
    gsMatrix<T> coefs(num_rows, m_dimension);
    coefs.setZero();

    coefs.row(0) = points.row(0);
    coefs.row(1) = coefs.row(0) + x(0)*T0.transpose();
    for (index_t i=0; i < num_rows-3; i++) {
        coefs(i+2,0) = x(i+1);
        coefs(i+2,1) = x(i+1+num_rows-3);
    }
    coefs.row(num_rows-1) = points.row(num_points-1);

    gsInfo << "Ridici body: \n" << coefs << "\n";

    gsBSpline<T> curve = gsBSpline<T> (*curveBasis, give(coefs));
    delete curveBasis;

    gsInfo << "Vysledna krivka:\n" << curve << "\n";

    return curve;
}
*/

template<class T>
gsBSpline<T> curveFittingWithBoundaryAndInputTangent(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kvfit, gsMatrix<T> T0)
{
    // number of points
    int num_points=points.rows();
    // dimension of points will be also later dimension of coefficients
    int m_dimension=points.cols();
    //gsInfo << "Dimenze: " << m_dimension << "\n";
    // basis function definition
    gsBSplineBasis<T> *curveBasis = new gsBSplineBasis<T>(kvfit);
    //number of basis functions
    unsigned int num_basis=curveBasis->size();


    //for computing the value of the basis function
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    //computing the values of the basis functions at some position
    curveBasis->eval_into(parameter_points,values);
    // which functions have been computed i.e. which are active
    curveBasis->active_into(parameter_points,actives);

    //gsInfo << "Vycislene bazove funkce:\n" << values << "\n\n";
    //gsInfo << "Aktivni bazove funkce:\n" << actives << "\n\n";

    //how many rows and columns has the A matrix and how many rows has the b vector
    int num_rows=num_basis;
    //left side matrix
    gsMatrix<T> m_A0(num_rows-2,num_rows-2);
    m_A0.setZero(); // ensure that all entries are zero in the beginning
    //right side vector (more dimensional!)
    gsMatrix<T> m_B0(num_rows-2,m_dimension);
    m_B0.setZero(); // enusure that all entris are zero in the beginning

    gsMatrix<T> R0values(1, num_points);
    gsBasisFun<T> R0 = curveBasis->function(0);
    R0.eval_into(parameter_points, R0values);
    //gsInfo << "R0:\n" << R0values << "\n";

    gsMatrix<T> R1values(1, num_points);
    gsBasisFun<T> R1 = curveBasis->function(1);
    R1.eval_into(parameter_points, R1values);
    //gsInfo << "R1:\n" << R1values << "\n";

    gsMatrix<T> RNvalues(1, num_points);
    gsBasisFun<T> RN = curveBasis->function(num_basis-1);
    RN.eval_into(parameter_points, RNvalues);
    //gsInfo << "RN:\n" << RNvalues << "\n";

    // building auxiliary matrix A and vector b for consequent building of the final matrix A and the vector b of the system of linear equations A*x==b
    T m_Bp = 0;
    T s1, s2, s3;
    for(index_t k=0;k<num_points;k++){
        if (m_dimension == 2) {
            s1 = T0(0)*points(k,0) + T0(1)*points(k,1);
            s2 = points(0,0)*T0(0) + points(0,1)*T0(1);
            s3 = T0(0)*points(num_basis-1,0) + T0(1)*points(num_basis-1,1);
        }
        else {
            s1 = T0(0)*points(k,0) + T0(1)*points(k,1) + T0(2)*points(k,2);
            s2 = points(0,0)*T0(0) + points(0,1)*T0(1) + points(0,2)*T0(2);
            s3 = T0(0)*points(num_basis-1,0) + T0(1)*points(num_basis-1,1) + T0(2)*points(num_basis-1,2);
        }
        m_Bp += R1values(k)*(s1 - s2*(R0values(k)+R1values(k)) - RNvalues(k)*s3);
        for(index_t i=0;i<actives.rows();i++){
            if ((actives(i,k) > 1) && ((unsigned) actives(i,k) < num_basis-1)) {
                m_B0.row(actives(i,k)-1) += values(i,k)*(points.row(k) - points.row(0)*R0values(k) - points.row(0)*R1values(k) - points.row(num_points-1)*RNvalues(k));
            }
            for(index_t j=0;j<actives.rows();j++){
                if ((actives(i,k) > 0) && (actives(j,k) > 0) && ((unsigned) actives(i,k) < num_basis-1) && ((unsigned) actives(j,k) < num_basis-1)) {
                    m_A0(actives(i,k)-1,actives(j,k)-1) += values(i,k)*values(j,k);
                }

            }
        }
    }

    //gsInfo << T0 << "\n";
    //gsInfo << m_A0 << "\n";
    //gsInfo << m_B0 << "\n";
    //gsInfo << m_B0.block(1,0,num_rows-3,1) << "\n";

    //left side matrix of the linear system
    gsMatrix<T> m_Zero(num_rows-3,num_rows-3);
    m_Zero.setZero(); // ensure that all entries are zero in the beginning
    gsMatrix<T> m_A(m_dimension*(num_rows-3)+1,m_dimension*(num_rows-3)+1);
    m_A.setZero(); // ensure that all entries are zero in the beginning
    if (m_dimension == 2) {
        m_A << (T0(0)*T0(0) + T0(1)*T0(1))*m_A0(0,0),       T0(0)*m_A0.block(0,1,1,num_rows-3),         T0(1)*m_A0.block(0,1,1,num_rows-3),
               T0(0)*m_A0.block(1,0,num_rows-3,1),          m_A0.block(1,1,num_rows-3,num_rows-3),      m_Zero,
               T0(1)*m_A0.block(1,0,num_rows-3,1),          m_Zero,                                     m_A0.block(1,1,num_rows-3,num_rows-3);
    }
    else {
        m_A << (T0(0)*T0(0) + T0(1)*T0(1) + T0(2)*T0(2))*m_A0(0,0),     T0(0)*m_A0.block(0,1,1,num_rows-3),         T0(1)*m_A0.block(0,1,1,num_rows-3),     T0(2)*m_A0.block(0,1,1,num_rows-3),
               T0(0)*m_A0.block(1,0,num_rows-3,1),                      m_A0.block(1,1,num_rows-3,num_rows-3),      m_Zero,                                 m_Zero,
               T0(1)*m_A0.block(1,0,num_rows-3,1),                      m_Zero,                                     m_A0.block(1,1,num_rows-3,num_rows-3),  m_Zero,
               T0(2)*m_A0.block(1,0,num_rows-3,1),                      m_Zero,                                     m_Zero,                                 m_A0.block(1,1,num_rows-3,num_rows-3);
    }

    // right hand side of the linear system
    gsMatrix<T> m_B(m_dimension*(num_rows-3)+1,1);
    m_B.setZero(); // enusure that all entris are zero in the beginning
    if (m_dimension == 2) {
        m_B << m_Bp, m_B0.block(1,0,num_rows-3,1), m_B0.block(1,1,num_rows-3,1);
    }
    else {
        m_B << m_Bp, m_B0.block(1,0,num_rows-3,1), m_B0.block(1,1,num_rows-3,1), m_B0.block(1,2,num_rows-3,1);
    }



    //gsInfo << "Soustava sestavena. \n";

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)
    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x = m_A.fullPivHouseholderQr().solve( m_B);

    //gsInfo << "Soustava vyresena.\n" << x << "\n" << x.size() << "\n";

    // finally generate the B-spline curve
    gsMatrix<T> coefs(num_rows, m_dimension);
    coefs.setZero();

    coefs.row(0) = points.row(0);
    coefs.row(1) = coefs.row(0) + x(0)*T0.transpose();
    for (index_t i=0; i < num_rows-3; i++) {
        if (m_dimension == 2) {
            coefs(i+2,0) = x(i+1);
            coefs(i+2,1) = x(i+1+num_rows-3);
        }
        else {
            coefs(i+2,0) = x(i+1);
            coefs(i+2,1) = x(i+1+num_rows-3);
            coefs(i+2,2) = x(i+1+2*(num_rows-3));
        }
    }
    coefs.row(num_rows-1) = points.row(num_points-1);

    //gsInfo << "Ridici body: \n" << coefs << "\n";

    gsBSpline<T> curve = gsBSpline<T> (*curveBasis, give(coefs));
    delete curveBasis;

    //gsInfo << "Vysledna krivka:\n" << curve << "\n";

    return curve;
}

template<class T>
gsBSpline<T> curveFittingWithEndDerivatives(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kvfit, gsMatrix<T> dersStart, gsMatrix<T> dersEnd) {

    //if (print_info) {
        gsInfo << "\nB-spline approximation respecting boundary derivatives (up to second derivative) with least squares (Piegl, Tiller):\n";
        gsInfo << "=====================================================================================================================\n";
    //}
    // number of points
    int num_points=points.rows();
    // dimension of points will be also later dimension of coefficients
    int m_dimension=points.cols();
    // basis function definition
    gsBSplineBasis<T> curveBasis = gsBSplineBasis<T>(kvfit);
    //number of basis functions
    unsigned int num_basis=curveBasis.size();
    int num_rows = num_basis - dersStart.rows() - dersEnd.rows();
    int num_nonZero = kvfit.degree() + 1;

    gsInfo << num_basis << "\n";

    //for computing the value of the basis function
    gsMatrix<T> values;
    gsMatrix<T> dersMatrixStart;
    gsMatrix<T> ders2MatrixStart;
    gsMatrix<T> dersMatrixEnd;
    gsMatrix<T> ders2MatrixEnd;
    gsMatrix<index_t> actives;
    gsMatrix<T> parStart(1,1); parStart << kvfit.first();
    gsMatrix<T> parEnd(1,1); parEnd << kvfit.last();

    //computing the values of the basis functions at some position
    curveBasis.eval_into(parameter_points,values);
    // which functions have been computed i.e. which are active
    curveBasis.active_into(parameter_points,actives);
    // evaluating derivatives of basis functions at the start and end
    curveBasis.deriv_into(parStart, dersMatrixStart);
    curveBasis.deriv2_into(parStart, ders2MatrixStart);
    curveBasis.deriv_into(parEnd, dersMatrixEnd);
    curveBasis.deriv2_into(parEnd, ders2MatrixEnd);

    gsInfo << values.rows() << ", " << values.cols() << "\n";
    gsInfo << dersMatrixEnd << "\n\n";
    gsInfo << ders2MatrixEnd << "\n\n";
    gsInfo << values << "\n\n";
    gsInfo << actives << "\n\n";

    //
    gsMatrix<T> final_curve_cp(num_basis, m_dimension);
    final_curve_cp.setZero(num_basis, m_dimension);

    // respecting boundary conditions at the start
    for (index_t i = 0; i < dersStart.rows(); i++) {
        if (i == 0) {
            final_curve_cp.row(0) = dersStart.row(0);
        }
        else if (i == 1) {
            final_curve_cp.row(1) = (dersStart.row(1) - dersMatrixStart(0) * final_curve_cp.row(0)) / dersMatrixStart(1);
        }
        else if (i == 2) {
            final_curve_cp.row(2) = (dersStart.row(2) - ders2MatrixStart(0) * final_curve_cp.row(0) - ders2MatrixStart(1) * final_curve_cp.row(1)) / ders2MatrixStart(2);
        }
        else {
            GISMO_ASSERT(dersStart.rows() > 2, "Boundary conditions only up to second derivative can be entered!");
        }
    }

    // respecting boundary conditions at the end
    for (index_t i = 0; i < dersEnd.rows(); i++) {
        if (i == 0) {
            final_curve_cp.row(num_basis-1) = dersEnd.row(0);
        }
        else if (i == 1) {
            //final_curve_cp(num_basis-2,0) = (dersEnd(1,0) - dersMatrixEnd(num_basis-1,0) * final_curve_cp(num_basis-1,0)) / dersMatrixEnd(num_basis-2,0);
            //final_curve_cp(num_basis-2,1) = (dersEnd(1,1) - dersMatrixEnd(num_basis-1,1) * final_curve_cp(num_basis-1,1)) / dersMatrixEnd(num_basis-2,1);
            //if (m_dimension > 2) final_curve_cp(num_basis-2,2) = (dersEnd(1,2) - dersMatrixEnd(num_basis-1,2) * final_curve_cp(num_basis-1,2)) / dersMatrixEnd(num_basis-2,2);
            //final_curve_cp.row(num_basis-2) = (dersEnd.row(1) - dersMatrixEnd(num_rows-1) * final_curve_cp.row(num_basis-1)) / dersMatrixEnd(num_rows-2);
            //final_curve_cp.row(num_basis-2) = (dersEnd.row(1) - dersMatrixEnd(num_basis-1) * final_curve_cp.row(num_basis-1)) / dersMatrixEnd(num_basis-2);
            final_curve_cp.row(num_basis-2) = (dersEnd.row(1) - dersMatrixEnd(num_nonZero-1) * final_curve_cp.row(num_basis-1)) / dersMatrixEnd(num_nonZero-2);
        }
        else if (i == 2) {
            //final_curve_cp(num_basis-3,0) = (dersEnd(2,0) - ders2MatrixEnd(num_basis-1,0) * final_curve_cp(num_basis-1,0) - ders2MatrixEnd(num_basis-2,0) * final_curve_cp(num_basis-2,0)) / ders2MatrixEnd(num_basis-3,0);
            //final_curve_cp(num_basis-3,1) = (dersEnd(2,1) - ders2MatrixEnd(num_basis-1,1) * final_curve_cp(num_basis-1,1) - ders2MatrixEnd(num_basis-2,1) * final_curve_cp(num_basis-2,1)) / ders2MatrixEnd(num_basis-3,1);
            //if (m_dimension > 2) final_curve_cp(num_basis-3,2) = (dersEnd(2,2) - ders2MatrixEnd(num_basis-1,2) * final_curve_cp(num_basis-1,2) - ders2MatrixEnd(num_basis-2,2) * final_curve_cp(num_basis-2,2)) / ders2MatrixEnd(num_basis-3,2);
            final_curve_cp.row(num_basis-3) = (dersEnd.row(2) - ders2MatrixEnd(num_nonZero-1) * final_curve_cp.row(num_basis-1) - ders2MatrixEnd(num_nonZero-2) * final_curve_cp.row(num_basis-2)) / ders2MatrixEnd(num_nonZero-3);
        }
        else {
            GISMO_ERROR("Boundary conditions only up to second derivative can be entered!");
        }
    }

    gsInfo << final_curve_cp << "\n\n";

    //how many rows and columns has the A matrix and how many rows has the b vector
    int indStart = dersStart.rows();
    int indEnd = num_basis - 1 - dersEnd.rows();
    //left side matrix
    gsMatrix<T> m_N(num_points-2,num_rows);
    m_N.setZero(num_points-2,num_rows);
    for (index_t i = 1; i < num_points-1; i++)
        for (index_t j = 0; j < values.rows(); j++) {
            if ((actives(j,i) >= indStart) && (actives(j,i) <= indEnd))
                m_N(i-1,actives(j,i)-dersStart.rows()) = values(j,i);
        }
    gsInfo << m_N.rows() << ", " << m_N.cols() << "\n";
    gsInfo << m_N << "\n\n";
    gsMatrix<T> m_A(num_rows, num_rows);
    m_A = (m_N.transpose()) * m_N;
    gsInfo << m_A.rows() << ", " << m_A.cols() << "\n";

    //right side vector (more dimensional!)
    gsMatrix<T> sum1(1, m_dimension), sum2(1, m_dimension);
    gsMatrix<T> m_R(num_points-2, m_dimension);
    m_R.setZero(num_points-2, m_dimension);
    for (index_t r = 1; r < num_points-1; r++) {
        sum1.setZero(1, m_dimension); sum2.setZero(1, m_dimension);
        for (index_t j = 0; j < dersStart.rows(); j++) {
            if (actives(j,r) < indStart)
                sum1 = sum1 + values(j,r) * final_curve_cp.row(actives(j,r));
        }
        gsInfo << sum1 << "\n";
        for (index_t j = 0; j < dersEnd.rows(); j++) {
            //if (actives(num_rows-1-j,r) > indEnd) sum2 = sum2 + values(num_rows-1-j,r) * final_curve_cp.row(actives(num_rows-1-j,r));
            //if (actives(num_basis-1-j,r) > indEnd)
            if (actives(actives.rows()-1-j,r) > indEnd)
                //sum2 = sum2 + values(num_basis-1-j,r) * final_curve_cp.row(actives(num_basis-1-j,r));
                sum2 = sum2 + values(values.rows()-1-j,r) * final_curve_cp.row(actives(actives.rows()-1-j,r));
        }
        gsInfo << sum2 << "\n";
        m_R.row(r-1) = points.row(r) - sum1 - sum2;
    }
    gsInfo << m_R.rows() << ", " << m_R.cols() << "\n";
    gsInfo << m_R << "\n\n";
    gsMatrix<T> m_B(num_rows,m_dimension);
    m_B.setZero();
    m_B = (m_N.transpose()) * m_R;
    gsInfo << m_B.rows() << ", " << m_B.cols() << "\n";

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)
    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x=m_A.fullPivHouseholderQr().solve( m_B);
    gsInfo << x << "\n";

    // finally generate the B-spline curve
    //gsMatrix<T> coefs(num_rows, m_dimension);
    //coefs.row(0) = points.row(0);
    for (index_t i=0; i < num_rows; i++) {
        final_curve_cp.row(indStart+i) = x.row(i);
    }
    //coefs.row(num_rows-1) = points.row(num_points-1);

    gsInfo << final_curve_cp << "\n\n";

    gsBSpline<T> curve = gsBSpline<T> (curveBasis, give(final_curve_cp));
    //if (print_info) {
        gsInfo << "Final curve:\n" << curve << "\n";
    //}

    return curve;
}



// Returns chordal parametrization for given point in plane or space.
// Requires "points" as gsMatrix of size nxd, where n is the number of given points and d is the dimension
template<class T>
gsMatrix<T> chordalParameterization(gsMatrix<T> &points) {

    T l;
    int num_points=points.rows();
    int m_dimension = points.cols();
    gsMatrix<T> pars(1, num_points);
    pars.setZero();


    for (int i = 0; i < num_points-1; i++) {
        l = 0;
        for (int j = 0; j< m_dimension; j++) {
            l +=(points(i+1, j) - points(i, j))*(points(i+1, j) - points(i, j));
        }
        //l = sqrt((points(i+1, 0)-points(i, 0))*(points(i+1, 0)-points(i, 0)) + (points(i+1, 1)-points(i, 1))*(points(i+1, 1)-points(i, 1)));
        l = sqrt(l);
        pars(i+1) = pars(i) + l;
    }

    return pars/pars(num_points-1);

}

// Returns centripetal parametrization for given point in plane or space.
// Requires "points" as gsMatrix of size nxd, where n is the number of given points and d is the dimension
template<class T>
gsMatrix<T> centripetalParameterization(gsMatrix<T> &points) {

    T l;
    int num_points = points.rows();
    int m_dimension = points.cols();
    gsMatrix<T> pars(1, num_points);
    pars.setZero();


    for (int i = 0; i < num_points-1; i++) {
        l = 0;
        for (int j = 0; j< m_dimension; j++) {
            l +=(points(i+1, j) - points(i, j))*(points(i+1, j) - points(i, j));
        }
        //l = sqrt(sqrt((points(i+1, 0)-points(i, 0))*(points(i+1, 0)-points(i, 0)) + (points(i+1, 1)-points(i, 1))*(points(i+1, 1)-points(i, 1))));
        l = sqrt(sqrt(l));
        pars(i+1) = pars(i) + l;
    }

    return pars/pars(num_points-1);

}

// Finding intersection of a B-spline curve with B-spline surface via Newton-Raphson iteration
// par_surf1 defines parameter value of second parameter of surface1, for which the curve on surface1 is selected to compute intersection with surface2
template<class T>
int curveSurfaceIntersectionviaNR(gsTensorBSpline<2, T> & surface1, T par_surf1, gsTensorBSpline<2, T> & surface2, gsVector<T> initial_solution, gsVector<T> & solution, T tol, int maxiter, int & num_iter, bool print_info) {

    gsMatrix<T> pars_for_surface2(2, 1);
    gsMatrix<T> pars_for_surface1(2, 1);
    gsMatrix<T> surface2_evaluated(3,1);
    gsMatrix<T> surface1_evaluated(3,1);
    gsMatrix<T> surface2_ders(6, 1);
    gsMatrix<T> surface1_ders(6, 1);
    gsMatrix<T> curve_ders(3, 1);
    gsMatrix<T> m_B(3, 1);
    gsMatrix<T> m_J(3, 3);
    T error;

    pars_for_surface2(0) = initial_solution(0);
    pars_for_surface2(1) = initial_solution(1);
    pars_for_surface1(0) = initial_solution(2);
    pars_for_surface1(1) = par_surf1;

    surface2.eval_into(pars_for_surface2, surface2_evaluated);
    surface1.eval_into(pars_for_surface1, surface1_evaluated);
    //curve.eval_into(par_for_curve, curve_evaluated);
    m_B = surface2_evaluated - surface1_evaluated;
    error = sqrt(pow(m_B(0), 2) + pow(m_B(1), 2) + pow(m_B(2), 2));

    if (print_info)
    {
        gsInfo << "==============================\n";
        gsInfo << "Newton-Raphson iteration: 0\n";
        gsInfo << "==============================\n";
        gsInfo << "Point on surface: " << surface2_evaluated(0) << ", " << surface2_evaluated(1) << ", " << surface2_evaluated(2) << "\n";
        gsInfo << "Point on curve: " << surface1_evaluated(0) << ", " << surface1_evaluated(1)<< ", " << surface1_evaluated(2) << "\n";
    }

    //if (surface1_evaluated(1) > 0) {
    //    pars_for_surface2(1) = 1-pars_for_surface2(1);
    //}

    for (index_t i = 0; i < maxiter; i++) {

        num_iter = i;

        if (error < tol) {
            if (print_info)
            {
                gsInfo << "===================================\n";
                gsInfo << "Newton-Raphson iteration converged!\n";
                gsInfo << "===================================\n";
                gsInfo << "Solution: " << solution << "\n";
                gsInfo << "Error of the solution: " << error << "\n";
            }
            return 0;
        }
        else {
            surface2.deriv_into(pars_for_surface2, surface2_ders);
            //gsInfo << "Vycisleni derivaci plochy v parametrech: " << surface2_ders << "\n";
            surface1.deriv_into(pars_for_surface1, surface1_ders);
            //curve.deriv_into(par_for_curve, curve_ders);
            //gsInfo << "Vycisleni derivace druhe plochy v parametru: " << surface1_ders << "\n";
            curve_ders(0) = surface1_ders(0);
            curve_ders(1) = surface1_ders(2);
            curve_ders(2) = surface1_ders(4);
            //gsInfo << "Vycisleni derivace krivky v parametru: " << curve_ders << "\n";
            m_J(0,0) = surface2_ders(0,0);
            m_J(0,1) = surface2_ders(1,0);
            m_J(0,2) = -curve_ders(0);
            m_J(1,0) = surface2_ders(2,0);
            m_J(1,1) = surface2_ders(3,0);
            m_J(1,2) = -curve_ders(1);
            m_J(2,0) = surface2_ders(4,0);
            m_J(2,1) = surface2_ders(5,0);
            m_J(2,2) = -curve_ders(2);

            //gsInfo << "m_J = " << m_J << "\n";
            //gsInfo << "m_B = " << m_B << "\n";

            gsMatrix<T> x (m_B.rows(), m_B.cols());
            x=m_J.fullPivHouseholderQr().solve( m_B);

            pars_for_surface2(0) -= x(0);
            pars_for_surface2(1) -= x(1);
            pars_for_surface1(0) -= x(2);
            //par_for_curve(0) += x(2);
            pars_for_surface2(0) = (pars_for_surface2(0) < 0) ? 0.0 : pars_for_surface2(0);
            pars_for_surface2(0) = (pars_for_surface2(0) > 1) ? 1.0 : pars_for_surface2(0);
            pars_for_surface2(1) = (pars_for_surface2(1) < 0) ? 0.0 : pars_for_surface2(1);
            pars_for_surface2(1) = (pars_for_surface2(1) > 1) ? 1.0 : pars_for_surface2(1);
            pars_for_surface1(0) = (pars_for_surface1(0) < 0) ? 0.0 : pars_for_surface1(0);
            pars_for_surface1(0) = (pars_for_surface1(0) > 1) ? 1.0 : pars_for_surface1(0);

            solution(0) = pars_for_surface2(0);
            solution(1) = pars_for_surface2(1);
            solution(2) = pars_for_surface1(0);

            if (print_info)
            {
                gsInfo << "x = " << x << "\n";
                gsInfo << "solution = " << solution << "\n";
                gsInfo << "pars_for_surface1 = " << pars_for_surface1 << "\n";
                gsInfo << "pars_for_surface2 = " << pars_for_surface2 << "\n";
            }

            surface2.eval_into(pars_for_surface2, surface2_evaluated);
            //gsInfo << "Vycisleni plochy v parametrech: " << surface2_evaluated << "\n";
            surface1.eval_into(pars_for_surface1, surface1_evaluated);
            //curve.eval_into(par_for_curve, curve_evaluated);
            //gsInfo << "Vycisleni krivky v parametru: " << surface1_evaluated << "\n";
            m_B = surface2_evaluated - surface1_evaluated;

            error = sqrt(pow(m_B(0), 2) + pow(m_B(1), 2) + pow(m_B(2), 2));

            if (print_info)
            {
                gsInfo << "==============================\n";
                gsInfo << "Newton_Raphson iteration: " << i+1 << "\n";
                gsInfo << "==============================\n";
                gsInfo << "Point on surface: " << surface2_evaluated(0) << ", " << surface2_evaluated(1) << ", " << surface2_evaluated(2) << "\n";
                gsInfo << "Point on curve: " << surface1_evaluated(0) << ", " << surface1_evaluated(1)<< ", " << surface1_evaluated(2) << "\n";
                gsInfo << "Error of the current solution: " << error << "\n";
            }
        }
    }

    if (print_info)
    {
        gsInfo << "Solution does not converged in " << maxiter << "iterations.\n Current solution: " << solution << "\n Its error: " << error << "\n";
    }
    return 0;
}

// Finding intersection of a B-spline curve with sphere via Newton-Raphson iteration
// par_surf1 defines parameter value of second parameter of surface1, for which the curve on surface1 is selected to compute intersection with surface2
template<class T>
int curveSphereIntersectionviaNR(gsTensorBSpline<2, T> & surface1, T par_surf1, T initial_solution, T radius, T & solution, T tol, int maxiter, int & num_iter, bool print_info) {

    gsMatrix<T> pars_for_surface1(2, 1);
    gsMatrix<T> surface1_evaluated(3,1);
    gsMatrix<T> surface1_ders(6, 1);
    gsMatrix<T> curve_ders(3, 1);
    T ft;
    T derft;
    T error;

    pars_for_surface1(0) = initial_solution;
    pars_for_surface1(1) = par_surf1;
    solution = initial_solution;

    surface1.eval_into(pars_for_surface1, surface1_evaluated);
    ft = pow(surface1_evaluated(0), 2) + pow(surface1_evaluated(1), 2) + pow(surface1_evaluated(2), 2) - pow(radius, 2);

    error = fabs(pow(surface1_evaluated(0), 2) + pow(surface1_evaluated(1), 2) + pow(surface1_evaluated(2), 2) - pow(radius, 2));

    if (print_info)
    {
        gsInfo << "==============================\n";
        gsInfo << "Newton-Raphson iteration: 0\n";
        gsInfo << "==============================\n";
        gsInfo << "Point on curve: " << surface1_evaluated(0) << ", " << surface1_evaluated(1)<< ", " << surface1_evaluated(2) << "\n";
    }

    //if (surface1_evaluated(1) > 0) {
    //    pars_for_surface2(1) = 1-pars_for_surface2(1);
    //}

    for (index_t i = 0; i < maxiter; i++) {

        num_iter = i;

        if (error < tol) {
            if (print_info)
            {
                gsInfo << "===================================\n";
                gsInfo << "Newton-Raphson iteration converged!\n";
                gsInfo << "===================================\n";
                gsInfo << "Solution: " << solution << "\n";
                gsInfo << "Error of the solution: " << error << "\n";
            }
            return 0;
        }
        else {
            surface1.deriv_into(pars_for_surface1, surface1_ders);
            //curve.deriv_into(par_for_curve, curve_ders);
            //gsInfo << "Vycisleni derivace druhe plochy v parametru: " << surface1_ders << "\n";
            curve_ders(0) = surface1_ders(0);
            curve_ders(1) = surface1_ders(2);
            curve_ders(2) = surface1_ders(4);
            //gsInfo << "Vycisleni derivace krivky v parametru: " << curve_ders << "\n";
            derft = 2 * surface1_evaluated(0) * curve_ders(0) + 2 * surface1_evaluated(1) * curve_ders(1) + 2 * surface1_evaluated(2) * curve_ders(2);

            //gsInfo << "m_J = " << m_J << "\n";
            //gsInfo << "m_B = " << m_B << "\n";

            //gsMatrix<T> x (m_B.rows(), m_B.cols());
            //x=m_J.fullPivHouseholderQr().solve( m_B);
            T x;
            x = -ft/derft;

            pars_for_surface1(0) += x;
            //par_for_curve(0) += x(2);
            pars_for_surface1(0) = (pars_for_surface1(0) < 0) ? 0.0 : pars_for_surface1(0);
            pars_for_surface1(0) = (pars_for_surface1(0) > 1) ? 1.0 : pars_for_surface1(0);

            solution = pars_for_surface1(0);

            if (print_info)
            {
                gsInfo << "x = " << x << "\n";
                gsInfo << "solution = " << solution << "\n";
                gsInfo << "pars_for_surface1 = " << pars_for_surface1 << "\n";
            }

            surface1.eval_into(pars_for_surface1, surface1_evaluated);
            ft = pow(surface1_evaluated(0), 2) + pow(surface1_evaluated(1), 2) + pow(surface1_evaluated(2), 2) - pow(radius, 2);
            //curve.eval_into(par_for_curve, curve_evaluated);
            //gsInfo << "Vycisleni krivky v parametru: " << surface1_evaluated << "\n";

            error = fabs(pow(surface1_evaluated(0), 2) + pow(surface1_evaluated(1), 2) + pow(surface1_evaluated(2), 2) - pow(radius, 2));

            if (print_info)
            {
                gsInfo << "==============================\n";
                gsInfo << "Newton_Raphson iteration: " << i+1 << "\n";
                gsInfo << "==============================\n";
                gsInfo << "Point on curve: " << surface1_evaluated(0) << ", " << surface1_evaluated(1)<< ", " << surface1_evaluated(2) << "\n";
                gsInfo << "Error of the current solution: " << error << "\n";
            }
        }
    }

    if (print_info)
    {
        gsInfo << "Solution does not converged in " << maxiter << "iterations.\n Current solution: " << solution << "\n Its error: " << error << "\n";
    }
    return 0;
}

// Computes lift surface from a given set of spatial B-spline curves
template<class T>
int computeLoftSurface(gsBSpline<T> defcurves[], gsKnotVector<T> kvcurves, int num_bladeprofiles, gsKnotVector<T> kvloft, gsMatrix<T> sec_parameters, gsTensorBSpline<2,T> & finalsurface) {

    int num_points_section;
    gsBSplineBasis<T> loftBasis;
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    num_points_section = defcurves[0].coefs().rows();
    loftBasis = gsBSplineBasis<T>(kvloft);
    int num_loftbasis=loftBasis.size();
    //gsInfo << "Pocet bazovych funkci pro loft: " << num_loftbasis << "\n";

    loftBasis.eval_into(sec_parameters, values);
    loftBasis.active_into(sec_parameters, actives);

    //gsInfo << "Hodnoty bazovych funkci:\n" << values << "\n";
    //gsInfo << "Aktivni bazove funkce:\n" << actives << "\n";

    gsMatrix<T> m_A(num_loftbasis, num_loftbasis);
    m_A.setZero();
    gsMatrix<T> m_B(num_loftbasis, 3);
    m_B.setZero();

    for (index_t i = 0; i < values.cols(); i++) {
        for (index_t k = 0; k < actives.rows(); k++) {
            m_A(i, actives(k, i)) = values(k, i);
        }
    }

    //gsInfo << "Matice soustavy:\n" << m_A << "\n";

    //gsMatrix<T> controlnet (num_bladeprofiles*num_points_section, 3);
    gsMatrix<T> controlnet (num_loftbasis*num_points_section, 3);

    gsMatrix<T> x (m_B.rows(), m_B.cols());
    for (index_t j = 0; j < num_points_section; j ++) {
        //for (index_t i = 0; i < num_bladeprofiles; i++) {
        for (index_t i = 0; i < num_loftbasis; i++) {
            m_B.row(i) = defcurves[i].coefs().row(j);
        }
        //gsInfo << "Prava strana:\n" << m_B << "\n";

        x = m_A.fullPivHouseholderQr().solve( m_B);

        //gsInfo << "Reseni soustavy:\n" << x << "\n";

        //for (index_t i = 0; i < num_bladeprofiles; i++) {
        for (index_t i = 0; i < num_loftbasis; i++) {
            //controlnet.row(j*num_bladeprofiles+i) = x.row(i);
            controlnet.row(j*num_loftbasis+i) = x.row(i);
        }
    }

    //gsInfo << "Ridici sit:\n" << controlnet << "\n";

    gsTensorBSplineBasis<2, T> basis(kvloft, kvcurves);
    finalsurface = gsTensorBSpline<2, T> (basis, controlnet);

    return 0;
}
template<class T>
int computeLoftSurface(std::vector<gsBSpline<T> > defcurves, gsKnotVector<T> kvcurves, int num_bladeprofiles, gsKnotVector<T> kvloft, gsMatrix<T> sec_parameters, gsTensorBSpline<2,T> & finalsurface) {

    int num_points_section;
    gsBSplineBasis<T> loftBasis;
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    num_points_section = defcurves[0].coefs().rows();
    loftBasis = gsBSplineBasis<T>(kvloft);
    //gsInfo << "loftBasis" << loftBasis << "\n";
    int num_loftbasis=loftBasis.size();
    //gsInfo << "Pocet bazovych funkci pro loft: " << num_loftbasis << "\n";

    loftBasis.eval_into(sec_parameters, values);
    loftBasis.active_into(sec_parameters, actives);

    //gsInfo << "Hodnoty bazovych funkci:\n" << values << "\n";
    //gsInfo << "Aktivni bazove funkce:\n" << actives << "\n";

    gsMatrix<T> m_A(num_loftbasis, num_loftbasis);
    m_A.setZero();
    gsMatrix<T> m_B(num_loftbasis, 3);
    m_B.setZero();

    for (index_t i = 0; i < values.cols(); i++) {
        for (index_t k = 0; k < actives.rows(); k++) {
            m_A(i, actives(k, i)) = values(k, i);
        }
    }

    //gsInfo << "Matice soustavy:\n" << m_A << "\n";

    //gsMatrix<T> controlnet (num_bladeprofiles*num_points_section, 3);
    gsMatrix<T> controlnet (num_loftbasis*num_points_section, 3);

    gsMatrix<T> x (m_B.rows(), m_B.cols());
    for (index_t j = 0; j < num_points_section; j ++) {
        //for (index_t i = 0; i < num_bladeprofiles; i++) {
        for (index_t i = 0; i < num_loftbasis; i++) {
            m_B.row(i) = defcurves[i].coefs().row(j);
        }
        //gsInfo << "Prava strana:\n" << m_B << "\n";

        x = m_A.fullPivHouseholderQr().solve( m_B);

        //gsInfo << "Reseni soustavy:\n" << x << "\n";

        //for (index_t i = 0; i < num_bladeprofiles; i++) {
        for (index_t i = 0; i < num_loftbasis; i++) {
            //controlnet.row(j*num_bladeprofiles+i) = x.row(i);
            controlnet.row(j*num_loftbasis+i) = x.row(i);
        }
    }

    //gsInfo << "Ridici sit:\n" << controlnet << "\n";

    gsTensorBSplineBasis<2, T> basis(kvloft, kvcurves);
    finalsurface = gsTensorBSpline<2, T> (basis, controlnet);

    return 0;
}

template<class T>
int computeLoftFuction(std::vector<gsBSpline<T> > defcurves, gsKnotVector<T> kvcurves, int num_bladeprofiles, gsKnotVector<T> kvloft, gsMatrix<T> sec_parameters, gsTensorBSpline<2,T> & finalsurface, int dimension) {

    int num_points_section;
    gsBSplineBasis<T> loftBasis;
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    num_points_section = defcurves[0].coefs().rows();
    loftBasis = gsBSplineBasis<T>(kvloft);
    //gsInfo << "loftBasis" << loftBasis << "\n";
    int num_loftbasis=loftBasis.size();
    //gsInfo << "Pocet bazovych funkci pro loft: " << num_loftbasis << "\n";

    loftBasis.eval_into(sec_parameters, values);
    loftBasis.active_into(sec_parameters, actives);

    //gsInfo << "Hodnoty bazovych funkci:\n" << values << "\n";
    //gsInfo << "Aktivni bazove funkce:\n" << actives << "\n";

    gsMatrix<T> m_A(num_loftbasis, num_loftbasis);
    m_A.setZero(num_loftbasis, num_loftbasis);
    gsMatrix<T> m_B(num_loftbasis, dimension);
    m_B.setZero(num_loftbasis, dimension);

    for (index_t i = 0; i < values.cols(); i++) {
        for (index_t k = 0; k < actives.rows(); k++) {
            m_A(i, actives(k, i)) = values(k, i);
        }
    }

    //gsInfo << "Matice soustavy:\n" << m_A << "\n";

    //gsMatrix<T> controlnet (num_bladeprofiles*num_points_section, 3);
    gsMatrix<T> controlnet (num_loftbasis*num_points_section, dimension);

    gsMatrix<T> x (m_B.rows(), m_B.cols());
    for (index_t j = 0; j < num_points_section; j ++) {
        //for (index_t i = 0; i < num_bladeprofiles; i++) {
        for (index_t i = 0; i < num_loftbasis; i++) {
            m_B.row(i) = defcurves[i].coefs().row(j);
        }
        //gsInfo << "Prava strana:\n" << m_B << "\n";

        x = m_A.fullPivHouseholderQr().solve( m_B);

        //gsInfo << "Reseni soustavy:\n" << x << "\n";

        //for (index_t i = 0; i < num_bladeprofiles; i++) {
        for (index_t i = 0; i < num_loftbasis; i++) {
            //controlnet.row(j*num_bladeprofiles+i) = x.row(i);
            controlnet.row(j*num_loftbasis+i) = x.row(i);
        }
    }

    //gsInfo << "Ridici sit:\n" << controlnet << "\n";

    gsTensorBSplineBasis<2, T> basis(kvloft, kvcurves);
    finalsurface = gsTensorBSpline<2,T> (basis, controlnet);

    return 0;
}

// Fits polynomial of a given degree to given data (by least squares), return coefficients of this polynomial
template<class T>
int fitData(gsVector<T> datax, gsVector<T> datay, int degree, gsMatrix<T> & polynomialcoefficients) {

    int num_points;

    gsMatrix<T> m_A(degree+1, degree+1);
    m_A.setZero();
    gsMatrix<T> m_B(degree+1, 1);
    m_B.setZero();

    num_points = datax.size();
    //gsInfo << datax.size() << "\n";
    //gsInfo << pow(datax(1), 0) << "\n";
    for (index_t i = 0; i < num_points; i++) {
        for (index_t j = 0; j < degree+1; j++) {
            for (index_t k = j; k < degree+1; k++) {
                m_A(j,k) += pow(datax(i), j+k);
                if (j != k) {
                    m_A(k,j) = m_A(j,k);
                }
            }
            m_B(j) += pow(datax(i), j) * datay(i);
        }
    }
    //gsInfo << "Matice soustavy:\n" << m_A << "\n";
    //gsInfo << "Prava strana:\n" << m_B << "\n";

    gsMatrix<T> x(m_B.rows(), m_B.cols());
    x = m_A.fullPivHouseholderQr().solve( m_B);
    //gsInfo << x << "\n";
    //gsInfo << x.transpose() << "\n";

    polynomialcoefficients = x.transpose();

    return 0;
}

// Computes error of fitting polynomial to given data, requires coefficients of this polynomial computed by fitData()
template<class T>
int fitError(gsVector<T> datax, gsVector<T> datay, gsMatrix<T> polynomialcoefficients, T & error) {

    int num_points;
    int order;
    T x;
    T pom;

    error = 0;
    num_points = datax.size();
    order = polynomialcoefficients.size();
    for (index_t i = 0; i < num_points; i++) {
        x = datax(i);
        pom = polynomialcoefficients(0);
        for (index_t j = 1; j < order; j++) {
            pom += polynomialcoefficients(j) * x;
            x = x*datax(i);
        }
        //gsInfo << pom << "\n";
        error += pow((pom - datay(i)), 2);
    }

    return 0;
}

// Extrapolates given data - used for extrapolating given data for blade profiles
// Computes fitting to given data with increasing degree of polynomial, until the prescribed error is under tolerance,
// or maximal number of iteration is reached. Then, for obtained polynomial fit computes the extrapoled value.
template<class T>
int extrapolateData(gsVector<T> datax, gsVector<T> datay, T r, T & extra_data, T tol, int max_iter) {

    gsMatrix<T> polynomialcoefficients(1, max_iter+1);
    T error = 1000;
    int degree=0;
    T x;

    gsInfo << "Extrapolation of blade input data by fitting:\n";
    gsInfo << "=============================================\n";
    if (max_iter > datax.size() - 1) {
        max_iter = datax.size() - 1;
    }
    for (index_t i = 0; i < max_iter; i++) {
        degree = i+1;
        if (error < tol) {
            gsInfo << "Final fitting error for degree " << degree-1 << " polynomial is " << error << "\n\n";
            break;
        }
        else {
            fitData(datax, datay, degree, polynomialcoefficients);
            fitError(datax, datay, polynomialcoefficients, error);
            gsInfo << "Error for degree " << degree << " polynomial: " << error << "\n";
        }
    }
    x = r;
    extra_data = polynomialcoefficients(0);
    for (index_t j = 1; j < degree; j++) {
        extra_data += polynomialcoefficients(j) * x;
        x = x*r;
    }
    gsInfo << "Extrapoled function value for variable " << r << " is " << extra_data << "\n";

    return 0;
}

//least squares fitting for given data
template<class T>
gsMatrix<T> fitLeastSquares(gsMatrix<T> data)
 {

    int i,j,k;
    int degree = data.rows();
    gsMatrix<T> sigmaX(2*degree + 1,1);
    for (i=0;i<2*degree+1;i++)
    {
        sigmaX(i,0)=0;
        for (j=0;j<data.rows();j++)
            sigmaX(i,0)=sigmaX(i,0)+math::pow(data(j,0),i);
    }
    gsMatrix<T> normal(degree+1,degree+2);
    gsMatrix<T> coef(degree+1,1);

    for (i=0;i<=degree;i++)
        for (j=0;j<=degree;j++)
            normal(i,j)=sigmaX(i+j,0);
    gsMatrix<T> sigmaY(degree+1,1);
    for (i=0;i<degree+1;i++)
    {
        sigmaY(i,0)=0;
        for (j=0;j<data.rows();j++)
        sigmaY(i,0)=sigmaY(i,0)+math::pow(data(j,0),i)*data(j,1);
    }
    for (i=0;i<=degree;i++)
        normal(i,degree+1)=sigmaY(i,0);
    degree=degree+1;


    for (i=0;i<degree;i++)
        for (k=i+1;k<degree;k++)
            if (normal(i,i)<normal(k,i))
                for (j=0;j<=degree;j++)
                {
                    T temp=normal(i,j);
                    normal(i,j)=normal(k,j);
                    normal(k,j)=temp;
                }

    for (i=0;i<degree-1;i++)
        for (k=i+1;k<degree;k++)
            {
                T t=normal(k,i)/normal(i,i);
                for (j=0;j<=degree;j++)
                    normal(k,j)=normal(k,j)-t*normal(i,j);
            }

    for (i=degree-1;i>=0;i--)
    {

        coef(i,0)=normal(i,degree);
        for (j=0;j<degree;j++)
            if (j!=i)
                coef(i,0)=coef(i,0)-normal(i,j)*coef(j,0);
        coef(i,0)=coef(i,0)/normal(i,i);

    }

    return coef;
    }

//Construction of (approximate) cylinder (represented as B-spline surface) with axis coincident with z-axis
// from p0 to p1 on z-axis with given radius
template<class T>
int constructCylinderParametrization(T p0, T p1, T radius, gsTensorBSpline<2, T> & cylinder) {

    T const k = 0.55228475; //constant for circle section

    gsMatrix<T> coefsclp(26, 3);

    coefsclp << radius, 0, p0,
                radius, 0, p1,
                radius, k * radius, p0,
                radius, k * radius, p1,
                k * radius, radius, p0,
                k * radius, radius, p1,
                0, radius, p0,
                0, radius, p1,
                -k * radius, radius, p0,
                -k * radius, radius, p1,
                -radius, k * radius, p0,
                -radius, k * radius, p1,
                -radius, 0, p0,
                -radius, 0, p1,
                -radius, -k * radius, p0,
                -radius, -k * radius, p1,
                -k * radius, -radius, p0,
                -k * radius, -radius, p1,
                0, -radius, p0,
                0, -radius, p1,
                k * radius, - radius, p0,
                k * radius, - radius, p1,
                radius, -k * radius, p0,
                radius, -k * radius, p1,
                radius, 0, p0,
                radius, 0, p1;
    gsKnotVector<T> kv1(0,1,3,4,3);
    gsKnotVector<T> kv2(0,1,0,2);
    gsTensorBSplineBasis<2,T> basis(kv2, kv1);
    cylinder =  gsTensorBSpline<2,T>  (basis, coefsclp);

    return 0;
}

// Finding the closest poit on the surface "surface" to the point "point" with the help of Newton method. The result
// of the method outputs in "solution" and represents the pair of parameters on "surface" corresponding to the closest
// point.
template<class T>
int minimizePointSurfaceDistanceviaNR(gsVector<T> point, gsTensorBSpline<2, T> & surface, gsVector<T> initial_solution, gsVector<T> & solution, T tol, int maxiter, int & num_iter, bool print_info) {

    gsMatrix<T> pars_for_surface(2, 1);
    gsMatrix<T> surface_evaluated(3,1);
    gsMatrix<T> surface_ders(6, 1);
    gsMatrix<T> surface_ders2(9, 1);
    gsMatrix<T> m_B(2, 1);
    gsMatrix<T> m_J(2, 2);
    T error;

    pars_for_surface(0) = initial_solution(0);
    pars_for_surface(1) = initial_solution(1);

    surface.eval_into(pars_for_surface, surface_evaluated);
    surface.deriv_into(pars_for_surface, surface_ders);

    m_B(0) = 2 * ((surface_evaluated(0) - point(0)) * surface_ders(0)) + 2 * ((surface_evaluated(1) - point(1)) * surface_ders(2)) + 2 * ((surface_evaluated(2) - point(2)) * surface_ders(4));
    m_B(1) = 2 * ((surface_evaluated(0) - point(0)) * surface_ders(1)) + 2 * ((surface_evaluated(1) - point(1)) * surface_ders(3)) + 2 * ((surface_evaluated(2) - point(2)) * surface_ders(5));
    error = sqrt(pow(m_B(0), 2) + pow(m_B(1), 2));

    if (print_info)
    {
        gsInfo << "==============================\n";
        gsInfo << "Newton-Raphson iteration: 0\n";
        gsInfo << "==============================\n";
        gsInfo << "Point on surface: " << surface_evaluated(0) << ", " << surface_evaluated(1) << ", " << surface_evaluated(2) << "\n";
        gsInfo << "Point: " << point(0) << ", " << point(1)<< ", " << point(2) << "\n";
    }

    //if (point(1) > 0) {
    //    pars_for_surface(1) = 1-pars_for_surface(1);
    //}

    for (index_t i = 0; i < maxiter; i++) {

        num_iter = i;

        if (error < tol) {
            if (print_info)
            {
                gsInfo << "===================================\n";
                gsInfo << "Newton-Raphson iteration converged!\n";
                gsInfo << "===================================\n";
                gsInfo << "Solution: " << solution << "\n";
                gsInfo << "Error of the solution: " << error << "\n";
            }
            return 0;
        }
        else {
            surface.eval_into(pars_for_surface, surface_evaluated);
            surface.deriv_into(pars_for_surface, surface_ders);
            surface.deriv2_into(pars_for_surface, surface_ders2);
            //gsInfo << "Vycisleni derivaci plochy v parametrech: " << surface2_ders << "\n";
            //gsInfo << "Vycisleni derivace druhe plochy v parametru: " << surface1_ders << "\n";
            //gsInfo << "Vycisleni derivace krivky v parametru: " << curve_ders << "\n";
            m_J(0,0) = 2 * (pow(surface_ders(0), 2) + (surface_evaluated(0) - point(0))*surface_ders2(0) + pow(surface_ders(2), 2) + (surface_evaluated(1) - point(1))*surface_ders2(3) + pow(surface_ders(4), 2) + (surface_evaluated(2) - point(2))*surface_ders2(6));
            m_J(0,1) = 2 * (surface_ders(0)*surface_ders(1) + (surface_evaluated(0) - point(0))*surface_ders2(2) + surface_ders(2)*surface_ders(3) + (surface_evaluated(1) - point(1))*surface_ders2(5) + surface_ders(4)*surface_ders(5) + (surface_evaluated(2) - point(2))*surface_ders2(8));
            m_J(1,0) = m_J(0,1);
            m_J(1,1) = 2 * (pow(surface_ders(1), 2) + (surface_evaluated(0) - point(0))*surface_ders2(1) + pow(surface_ders(3), 2) + (surface_evaluated(1) - point(1))*surface_ders2(4) + pow(surface_ders(5), 2) + (surface_evaluated(2) - point(2))*surface_ders2(7));;

            //gsInfo << "m_J = " << m_J << "\n";
            //gsInfo << "m_B = " << m_B << "\n";

            gsMatrix<T> x (m_B.rows(), m_B.cols());
            x=m_J.fullPivHouseholderQr().solve( m_B);

            pars_for_surface(0) -= x(0);
            pars_for_surface(1) -= x(1);

            pars_for_surface(0) = (pars_for_surface(0) < 0) ? 0.0 : pars_for_surface(0);
            pars_for_surface(0) = (pars_for_surface(0) > 1) ? 1.0 : pars_for_surface(0);
            pars_for_surface(1) = (pars_for_surface(1) < 0) ? 0.0 : pars_for_surface(1);
            pars_for_surface(1) = (pars_for_surface(1) > 1) ? 1.0 : pars_for_surface(1);

            solution(0) = pars_for_surface(0);
            solution(1) = pars_for_surface(1);

            if (print_info)
            {
                gsInfo << "x = " << x << "\n";
                gsInfo << "solution = " << solution << "\n";
                gsInfo << "pars_for_surface = " << pars_for_surface << "\n";
            }

            surface.eval_into(pars_for_surface, surface_evaluated);
            //gsInfo << "Vycisleni plochy v parametrech: " << surface2_evaluated << "\n";
            //gsInfo << "Vycisleni krivky v parametru: " << surface1_evaluated << "\n";
            m_B(0) = 2 * ((surface_evaluated(0) - point(0)) * surface_ders(0)) + 2 * ((surface_evaluated(1) - point(1)) * surface_ders(2)) + 2 * ((surface_evaluated(2) - point(2)) * surface_ders(4));
            m_B(1) = 2 * ((surface_evaluated(0) - point(0)) * surface_ders(1)) + 2 * ((surface_evaluated(1) - point(1)) * surface_ders(3)) + 2 * ((surface_evaluated(2) - point(2)) * surface_ders(5));

            error = sqrt(pow(m_B(0), 2) + pow(m_B(1), 2));

            if (print_info)
            {
                gsInfo << "==============================\n";
                gsInfo << "Newton_Raphson iteration: " << i+1 << "\n";
                gsInfo << "==============================\n";
                gsInfo << "Point on surface: " << surface_evaluated(0) << ", " << surface_evaluated(1) << ", " << surface_evaluated(2) << "\n";
                gsInfo << "Point: " << point(0) << ", " << point(1)<< ", " << point(2) << "\n";
                gsInfo << "Error of the current solution: " << error << "\n";
            }
        }
    }

    if (print_info)
    {
        gsInfo << "Solution does not converged in " << maxiter << "iterations.\n Current solution: " << solution << "\n Its error: " << error << "\n";
    }
    return 0;
}

// Cross product of two vectors in R^3
template<class T>
int cross(gsVector<T> vec1, gsVector<T> vec2, gsVector<T> &result) {

    if ((vec1.rows() != 3) || (vec2.rows() != 3)) {
        GISMO_ERROR("For cross product computation, two vecotrs in R^3 are required. Exiting.\n");
        return 1;
    }
    else {
        result(0) = vec1(1)*vec2(2) - vec2(1)*vec1(2);
        result(1) = vec2(0)*vec1(2) - vec2(2)*vec1(0);
        result(2) = vec1(0)*vec2(1) - vec2(0)*vec1(1);

        return 0;
    }
}

// Euclidean norm of a vector
template<class T>
T euclideanNorm(gsMatrix<T> v) {

    T norm = 0.0;

    if (v.cols() > 1) {
        for (index_t i = 0; i < v.cols(); i++) {
            norm += pow(v(i), 2);
        }
    }
    else {
        for (index_t i = 0; i < v.rows(); i++) {
            norm += pow(v(i), 2);
        }
    }

    return sqrt(norm);

}

template<class T>
T euclideanNorm(gsVector<T> v) {

    T norm = 0.0;

    for (index_t i = 0; i < v.size(); i++) {
        norm += pow(v(i), 2);
    }

    return sqrt(norm);

}

template<class T>
int discreteCoonsPatch(gsMatrix<T> a_cp, gsMatrix<T> b_cp, gsMatrix<T> c_cp, gsMatrix<T> d_cp, gsMatrix<T> & final_control_net, bool print_info) {

    int num1 = a_cp.rows();
    int num2 = c_cp.rows();
    T alfa = 0.0;
    T beta = 0.0;
    if (print_info) {
        gsInfo << "\nDiscrete Coons patch computation:\n";
        gsInfo << "=================================\n";
        gsInfo << "Dimensions of point net: " << num1 << ", " << num2 << "\n";
    }
    //gsMatrix<T> final_control_net(num1*num2, 3);
    //final_control_net.setZero();

    gsVector<T> vec1 = a_cp.row(0) - c_cp.row(0);
    gsVector<T> vec2 = b_cp.row(0) - c_cp.row(num2-1);
    gsVector<T> vec3 = a_cp.row(num1-1) - d_cp.row(0);
    gsVector<T> vec4 = b_cp.row(num1-1) - d_cp.row(num2-1);
    if (print_info) {
        gsInfo << "Check of compotatibility conditions: " << euclideanNorm(vec1) << ", " << euclideanNorm(vec2) << ", " << euclideanNorm(vec3) << ", " << euclideanNorm(vec4) << "\n";
    }

    if ((a_cp.rows() != b_cp.rows()) || (c_cp.rows() != d_cp.rows())) {
        GISMO_ERROR("1: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }
    //if ((a_cp.row(0) != c_cp.row(0)) || (b_cp.row(0) != c_cp.row(num2-1)) || (a_cp.row(num1-1) != d_cp.row(0)) || (b_cp.row(num1-1) != d_cp.row(num2-1))) {
    //    GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
    //    return 1;
    //}
    if ((euclideanNorm(vec1) > static_cast<T>(1e-8)) || (euclideanNorm(vec2) > static_cast<T>(1e-8)) || (euclideanNorm(vec3) > static_cast<T>(1e-8)) || (euclideanNorm(vec4) > static_cast<T>(1e-8))) {
        GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }

    // Initialization - filling boundary curves control points to the final control net
    for (index_t j = 0; j < num1; j++) {
        final_control_net.row(j) = a_cp.row(j);
    }
    for (index_t j = 0; j < num1; j++) {
        final_control_net.row((num2-1)*num1+j) = b_cp.row(j);
    }
    for (index_t i = 0; i < num2; i++) {
        final_control_net.row(i*num1) = c_cp.row(i);
    }
    for (index_t i = 0; i < num2; i++) {
        final_control_net.row(i*num1+num1-1) = d_cp.row(i);
    }
    //gsInfo << "Po inicializaci ...";
    //gsInfo << final_control_net << "\n";

    // surface 1
    for (index_t i = 1; i < num2-1; i++) {
        for (index_t j = 1; j < num1-1; j++) {
            alfa = (T) j/(num1-1);
            final_control_net.row(i*num1+j) = (1-alfa) * c_cp.row(i) + alfa * d_cp.row(i);
        }
    }
    //gsInfo << "Po doplneni prechodove plochy mezi c d ...\n";
    //gsInfo << final_control_net << "\n";

    // surface 2
    for (index_t j = 1; j < num1-1; j++) {
        for (index_t i = 1; i < num2-1; i++) {
            alfa = (T) i/(num2-1);
            final_control_net.row(i*num1+j) += (1-alfa) * a_cp.row(j) + alfa * b_cp.row(j);
        }
    }
    //gsInfo << "Po doplneni prechodove plochy mezi a b ...";
    //gsInfo << final_control_net << "\n";

    // surface 3
    for (index_t i = 1; i < num2-1; i++) {
        for (index_t j = 1; j < num1-1; j++) {
            alfa = (T) j/(num1-1);
            beta = (T) i/(num2-1);
            final_control_net.row(i*num1+j) -= (1-beta) * ((1-alfa) * a_cp.row(0) + alfa * a_cp.row(num1-1)) + beta * ((1-alfa) * b_cp.row(0) + alfa * b_cp.row(num1-1));
        }
    }
    gsInfo << "Po doplneni hyperbolickeho paraboloidu ...";
    gsInfo << final_control_net << "\n";

    return 0;

}

template<class T>
int discreteCoonsPatch(gsMatrix<T> a_cp, gsMatrix<T> b_cp, gsMatrix<T> c_cp, gsMatrix<T> d_cp, T w1, T w2, T w3, gsMatrix<T> & final_control_net, bool print_info) {

    int num1 = a_cp.rows();
    int num2 = c_cp.rows();
    T alfa = 0.0;
    T beta = 0.0;
    if (print_info) {
        gsInfo << "\nDiscrete Coons patch computation:\n";
        gsInfo << "=================================\n";
        gsInfo << "Dimensions of point net: " << num1 << ", " << num2 << "\n";
    }
    //gsMatrix<T> final_control_net(num1*num2, 3);
    //final_control_net.setZero();

    gsVector<T> vec1 = a_cp.row(0) - c_cp.row(0);
    gsVector<T> vec2 = b_cp.row(0) - c_cp.row(num2-1);
    gsVector<T> vec3 = a_cp.row(num1-1) - d_cp.row(0);
    gsVector<T> vec4 = b_cp.row(num1-1) - d_cp.row(num2-1);
    if (print_info) {
        gsInfo << "Check of compotatibility conditions: " << euclideanNorm(vec1) << ", " << euclideanNorm(vec2) << ", " << euclideanNorm(vec3) << ", " << euclideanNorm(vec4) << "\n";
    }

    if ((a_cp.rows() != b_cp.rows()) || (c_cp.rows() != d_cp.rows())) {
        GISMO_ERROR("1: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }
    //if ((a_cp.row(0) != c_cp.row(0)) || (b_cp.row(0) != c_cp.row(num2-1)) || (a_cp.row(num1-1) != d_cp.row(0)) || (b_cp.row(num1-1) != d_cp.row(num2-1))) {
    //    GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
    //    return 1;
    //}
    if ((euclideanNorm(vec1) > static_cast<T>(1e-8)) || (euclideanNorm(vec2) > static_cast<T>(1e-8)) || (euclideanNorm(vec3) > static_cast<T>(1e-8)) || (euclideanNorm(vec4) > static_cast<T>(1e-8))) {
        GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }

    // Initialization - filling boundary curves control points to the final control net
    for (index_t j = 0; j < num1; j++) {
        final_control_net.row(j) = a_cp.row(j);
    }
    for (index_t j = 0; j < num1; j++) {
        final_control_net.row((num2-1)*num1+j) = b_cp.row(j);
    }
    for (index_t i = 0; i < num2; i++) {
        final_control_net.row(i*num1) = c_cp.row(i);
    }
    for (index_t i = 0; i < num2; i++) {
        final_control_net.row(i*num1+num1-1) = d_cp.row(i);
    }
    //gsInfo << "Po inicializaci ...";
    //gsInfo << final_control_net << "\n";

    // surface 1
    for (index_t i = 1; i < num2-1; i++) {
        for (index_t j = 1; j < num1-1; j++) {
            alfa = (T) j/(num1-1);
            final_control_net.row(i*num1+j) = w2 * ((1-alfa) * c_cp.row(i) + alfa * d_cp.row(i));
        }
    }
    //gsInfo << "Po doplneni prechodove plochy mezi c d ...\n";
    //gsInfo << final_control_net << "\n";

    // surface 2
    for (index_t j = 1; j < num1-1; j++) {
        for (index_t i = 1; i < num2-1; i++) {
            alfa = (T) i/(num2-1);
            final_control_net.row(i*num1+j) += w1 * ((1-alfa) * a_cp.row(j) + alfa * b_cp.row(j));
        }
    }
    //gsInfo << "Po doplneni prechodove plochy mezi a b ...";
    //gsInfo << final_control_net << "\n";

    // surface 3
    for (index_t i = 1; i < num2-1; i++) {
        for (index_t j = 1; j < num1-1; j++) {
            alfa = (T) j/(num1-1);
            beta = (T) i/(num2-1);
            final_control_net.row(i*num1+j) -= -w3 * ((1-beta) * ((1-alfa) * a_cp.row(0) + alfa * a_cp.row(num1-1)) + beta * ((1-alfa) * b_cp.row(0) + alfa * b_cp.row(num1-1)));
        }
    }
    gsInfo << "Po doplneni hyperbolickeho paraboloidu ...";
    gsInfo << final_control_net << "\n";

    return 0;

}


template<class T>
int springModelPatch(gsMatrix<T> a_cp, gsMatrix<T> b_cp, gsMatrix<T> c_cp, gsMatrix<T> d_cp, gsMatrix<T> & final_control_net, bool print_info) {

    int num1 = a_cp.rows();
    int num2 = c_cp.rows();
    int m_dim = a_cp.cols();
    if (print_info) {
        gsInfo << "\nSpring model patch computation:\n";
        gsInfo << "=================================\n";
        gsInfo << "Dimensions of point net: " << num1 << ", " << num2 << "\n";
    }
    //gsMatrix<T> final_control_net(num1*num2, 3);
    //final_control_net.setZero();

    gsVector<T> vec1 = a_cp.row(0) - c_cp.row(0);
    gsVector<T> vec2 = b_cp.row(0) - c_cp.row(num2-1);
    gsVector<T> vec3 = a_cp.row(num1-1) - d_cp.row(0);
    gsVector<T> vec4 = b_cp.row(num1-1) - d_cp.row(num2-1);
    if (print_info) {
        gsInfo << "Check of compotatibility conditions: " << euclideanNorm(vec1) << ", " << euclideanNorm(vec2) << ", " << euclideanNorm(vec3) << ", " << euclideanNorm(vec4) << "\n";
    }

    if ((a_cp.rows() != b_cp.rows()) || (c_cp.rows() != d_cp.rows())) {
        GISMO_ERROR("1: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }
    //if ((a_cp.row(0) != c_cp.row(0)) || (b_cp.row(0) != c_cp.row(num2-1)) || (a_cp.row(num1-1) != d_cp.row(0)) || (b_cp.row(num1-1) != d_cp.row(num2-1))) {
    //    GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
    //    return 1;
    //}
    if ((euclideanNorm(vec1) > static_cast<T>(1e-8)) || (euclideanNorm(vec2) > static_cast<T>(1e-8)) || (euclideanNorm(vec3) > static_cast<T>(1e-8)) || (euclideanNorm(vec4) > static_cast<T>(1e-8))) {
        GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }

    // Initialization - filling boundary curves control points to the final control net
    for (index_t j = 0; j < num1; j++) {
        final_control_net.row(j) = a_cp.row(j);
    }
    for (index_t j = 0; j < num1; j++) {
        final_control_net.row((num2-1)*num1+j) = b_cp.row(j);
    }
    for (index_t i = 0; i < num2; i++) {
        final_control_net.row(i*num1) = c_cp.row(i);
    }
    for (index_t i = 0; i < num2; i++) {
        final_control_net.row(i*num1+num1-1) = d_cp.row(i);
    }

    // Setting up linear system for inner control points
    gsMatrix<T> m_A((num1-2)*(num2-2), (num1-2)*(num2-2));
    m_A.setZero();
    gsMatrix<T> m_B((num1-2)*(num2-2), m_dim);
    m_B.setZero();
    for (index_t i = 1; i < num2-1; i++) {
        for (index_t j = 1; j < num1-1; j++) {
            //gsInfo << i << ", " << j << "\n";
            m_A((i-1)*(num1-2)+j-1,(i-1)*(num1-2)+j-1) = 4;
            if ((i-1) == 0) {
                m_B.row((i-1)*(num1-2)+j-1) += final_control_net.row(j);
            }
            else {
                m_A((i-1)*(num1-2)+j-1,(i-2)*(num1-2)+j-1) = -1;
            }
            if ((i+1) == num2-1) {
                m_B.row((i-1)*(num1-2)+j-1) += final_control_net.row((num2-1)*num1+j);
            }
            else {
                m_A((i-1)*(num1-2)+j-1,i*(num1-2)+j-1) = -1;
            }
            if ((j-1) == 0) {
                m_B.row((i-1)*(num1-2)+j-1) += final_control_net.row(i*num1);
            }
            else {
                m_A((i-1)*(num1-2)+j-1,(i-1)*(num1-2)+j-2) = -1;
            }
            if ((j+1) == num1-1) {
                m_B.row((i-1)*(num1-2)+j-1) += final_control_net.row(i*num1+num1-1);
            }
            else {
                m_A((i-1)*(num1-2)+j-1,(i-1)*(num1-2)+j) = -1;
            }
        }
    }
    //gsInfo << m_A << "\n";
    //gsInfo << m_B << "\n";

    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x=m_A.fullPivHouseholderQr().solve( m_B);

    //gsInfo << x << "\n";

    for (index_t i = 1; i < num2-1; i++) {
        for (index_t j = 1; j < num1-1; j++) {
            final_control_net.row(i*num1+j) = x.row((i-1)*(num1-2)+j-1);
        }
    }

    return 0;
}

template<class T>
int projectToCylinder(gsMatrix<T> point, T radius, gsMatrix<T> & new_point) {

    T norm = sqrt(pow(point(0),2) + pow(point(1),2));
    T koef = radius/norm;
    //gsInfo << "koef r/sqrt = " << koef << "\n";
    new_point(0) = koef * point(0);
    new_point(1) = koef * point(1);
    new_point(2) = point(2);

    return 0;
}

template<class T>
int discreteCoonsPatchModified(gsMatrix<T> a_cp, gsMatrix<T> b_cp, gsMatrix<T> c_cp, gsMatrix<T> d_cp, gsMatrix<T> & final_control_net) {

    int num1 = a_cp.rows();
    int num2 = c_cp.rows();
    T alfa = 0.0;
    T beta = 0.0;
    T r1 = 0.0;
    T r2 = 0.0;
    T r = 0.0;
    gsMatrix<T> new_point(1,3);
    gsMatrix<T> pom;
    gsInfo << num1 << ", " << num2 << "\n";
    //gsMatrix<T> final_control_net(num1*num2, 3);
    //final_control_net.setZero();

    gsVector<T> vec1 = a_cp.row(0) - c_cp.row(0);
    gsVector<T> vec2 = b_cp.row(0) - c_cp.row(num2-1);
    gsVector<T> vec3 = a_cp.row(num1-1) - d_cp.row(0);
    gsVector<T> vec4 = b_cp.row(num1-1) - d_cp.row(num2-1);
    gsInfo << euclideanNorm(vec1) << "\n";
    gsInfo << euclideanNorm(vec2) << "\n";
    gsInfo << euclideanNorm(vec3) << "\n";
    gsInfo << euclideanNorm(vec4) << "\n";

    if ((a_cp.rows() != b_cp.rows()) || (c_cp.rows() != d_cp.rows())) {
        GISMO_ERROR("1: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }
    //if ((a_cp.row(0) != c_cp.row(0)) || (b_cp.row(0) != c_cp.row(num2-1)) || (a_cp.row(num1-1) != d_cp.row(0)) || (b_cp.row(num1-1) != d_cp.row(num2-1))) {
    //    GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
    //    return 1;
    //}
    if ((euclideanNorm(vec1) > static_cast<T>(1e-8)) || (euclideanNorm(vec2) > static_cast<T>(1e-8)) || (euclideanNorm(vec3) > static_cast<T>(1e-8)) || (euclideanNorm(vec4) > static_cast<T>(1e-8))) {
        GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }

    // Initialization - filling boundary curves control points to the final control net
    for (index_t j = 0; j < num1; j++) {
        final_control_net.row(j) = a_cp.row(j);
    }
    for (index_t j = 0; j < num1; j++) {
        final_control_net.row((num2-1)*num1+j) = b_cp.row(j);
    }
    for (index_t i = 0; i < num2; i++) {
        final_control_net.row(i*num1) = c_cp.row(i);
    }
    for (index_t i = 0; i < num2; i++) {
        final_control_net.row(i*num1+num1-1) = d_cp.row(i);
    }
    gsInfo << "Po inicializaci ...";
    gsInfo << final_control_net << "\n";

    // surface 1
    for (index_t i = 1; i < num2-1; i++) {
        for (index_t j = 1; j < num1-1; j++) {
            alfa = (T) j/(num1-1);
            final_control_net.row(i*num1+j) = (1-alfa) * c_cp.row(i) + alfa * d_cp.row(i);
        }
    }
    gsInfo << "Po doplneni prechodove plochy mezi c d ...\n";
    gsInfo << final_control_net << "\n";

    // surface 2
    for (index_t j = 1; j < num1-1; j++) {
        r1 = sqrt(pow(a_cp(j,0),2) + pow(a_cp(j,1),2));
        r2 = sqrt(pow(b_cp(j,0),2) + pow(b_cp(j,1),2));
        gsInfo << "r1 = " << r1 << "\n";
        gsInfo << "r2 = " << r2 << "\n";
        for (index_t i = 1; i < num2-1; i++) {
            alfa = (T) i/(num2-1);
            r = (1-alfa) * r1 + alfa * r2;
            gsInfo << "r = " << r << "\n";
            pom = (1-alfa) * a_cp.row(j) + alfa * b_cp.row(j);
            projectToCylinder(pom, r, new_point);
            final_control_net.row(i*num1+j) += new_point;
        }
    }
    gsInfo << "Po doplneni prechodove plochy mezi a b ...";
    gsInfo << final_control_net << "\n";

    // surface 3
    for (index_t i = 1; i < num2-1; i++) {
        for (index_t j = 1; j < num1-1; j++) {
            alfa = (T) j/(num1-1);
            beta = (T) i/(num2-1);
            final_control_net.row(i*num1+j) -= (1-beta) * ((1-alfa) * a_cp.row(0) + alfa * a_cp.row(num1-1)) + beta * ((1-alfa) * b_cp.row(0) + alfa * b_cp.row(num1-1));
        }
    }
    gsInfo << "Po doplneni hyperbolickeho paraboloidu ...";
    gsInfo << final_control_net << "\n";

    return 0;

}


// Fits the given points with a B-spline surface where the boundary control points are prescribed and preserved.
template<class T>
gsTensorBSpline<2,T> surfaceFittingLSQWithBoundary(gsMatrix<T> &points, gsMatrix<T> &parameter_points, gsKnotVector<T> kv1, gsKnotVector<T> kv2, gsMatrix<T> boundary_points, bool print_info)
{
    if (print_info) {
        gsInfo << "\nSurface fitting with least squares:\n";
        gsInfo << "===================================\n";
    }
    // number of points
    int num_points=points.rows();
    // dimension of points will be also later dimension of coefficients
    int m_dimension=points.cols();
    // basis function definition
    gsTensorBSplineBasis<2, T> surfaceBasis = gsTensorBSplineBasis<2, T>(kv1, kv2);
    //number of basis functions
    unsigned int num_basis=surfaceBasis.size();

    //for computing the value of the basis function
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    //computing the values of the basis functions at some position
    surfaceBasis.eval_into(parameter_points,values);
    // which functions have been computed i.e. which are active
    surfaceBasis.active_into(parameter_points,actives);

    gsMatrix<index_t> boundaryindices = surfaceBasis.allBoundary();
    //gsInfo << boundaryindices << "\n";
    int num_boundaryindices = boundaryindices.size();

    //how many rows and columns has the A matrix and how many rows has the b vector
    unsigned int num_rows=num_basis;
    //left side matrix
    gsMatrix<T> m_A(num_rows-num_boundaryindices,num_rows-num_boundaryindices);
    m_A.setZero(); // ensure that all entries are zero in the beginning
    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_rows-num_boundaryindices,m_dimension);
    m_B.setZero(); // enusure that all entris are zero in the beginning

    gsMatrix<T> boundaryvalues(num_boundaryindices, num_points);
    boundaryvalues.setZero();
    //gsBasisFun<T> basisfunction;
    gsMatrix<T> val(1, num_points);
    for (index_t i = 0; i < num_boundaryindices; i++) {
        gsBasisFun<T> basisfunction = surfaceBasis.function(boundaryindices(i));
        val.setZero();
        basisfunction.eval_into(parameter_points, val);
        boundaryvalues.row(i) = val;
    }

    // building the matrix A and the vector b of the system of linear equations A*x==b(uses automatically the information of closed or not)
    int q = 0;
    bool inner = true;
    gsMatrix<int> innerindices(2, num_basis);
    innerindices.setZero();
    for (unsigned i = 0; i < num_basis; i++) {
        inner = true;
        for (index_t k = 0; k < num_boundaryindices; k++) {
            if (i == boundaryindices(k)) {
                inner = false;
            }
        }
        if (inner) {
            innerindices(0,i) = i;
            innerindices(1,i) = q;
            q++;
        }
    }

    bool inner1 = true;
    bool inner2 = true;
    for(index_t k=0;k<num_points;k++){
        for(index_t i=0;i<actives.rows();i++){
            inner1 = true;
            inner2 = true;
            for (index_t l = 0; l < num_boundaryindices; l++) {
                if (actives(i,k) == boundaryindices(l)) {
                    inner1 = false;
                }
            }
            if (inner1) {
                m_B.row(innerindices(1, actives(i,k))) += values(i,k)*points.row(k);
                for (index_t l = 0; l < num_boundaryindices; l++) {
                    m_B.row(innerindices(1, actives(i,k))) -= values(i,k)*boundary_points.row(l)*boundaryvalues(l,k);
                }
            }
            for(index_t j=0;j<actives.rows();j++){
                inner2 = true;
                for (index_t l = 0; l < num_boundaryindices; l++) {
                    if (actives(i,k) == boundaryindices(l)) {
                        inner1 = false;
                    }
                    if (actives(j,k) == boundaryindices(l)) {
                        inner2 = false;
                    }
                }
                if (inner1 && inner2) {
                    m_A(innerindices(1,actives(i,k)),innerindices(1,actives(j,k))) += values(i,k)*values(j,k);
                }

            }
        }
    }

    //gsInfo << "Soustava sestavena. \n";

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)
    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x=m_A.fullPivHouseholderQr().solve( m_B);

    //gsInfo << "Soustava vyresena.\n" << x << "\n";
    //gsInfo << "m_A:\n" << m_A << "\n";
    //gsInfo << "m_B:\n" << m_B << "\n";

    // finally generate the B-spline curve
    gsMatrix<T> coefs(num_rows, m_dimension);
    int inner_index = -1;
    for (unsigned i=0; i < num_rows; i++) {
        inner = true;
        for (index_t l = 0; l < num_boundaryindices; l++) {
            if (i == boundaryindices(l)) {
                inner = false;
                inner_index = l;
            }
        }
        if (inner) {
            coefs.row(i) = x.row(innerindices(1,i));
        }
        else {
            coefs.row(i) = boundary_points.row(inner_index);
        }
    }

    gsTensorBSpline<2,T> surface = gsTensorBSpline<2, T> (surfaceBasis, coefs);
    //delete surfaceBasis;

    if (print_info) {
        gsInfo << "Final approximation surface (surfaceFittingLSQWithBoundary):\n" << surface << "\n";
    }

    return surface;
}


template<class T>
int discreteCoonsPatch3D(gsTensorBSpline<2,T> a, gsTensorBSpline<2,T> b, gsTensorBSpline<2,T> c, gsTensorBSpline<2,T> d, gsTensorBSpline<2,T> e, gsTensorBSpline<2,T> f, gsTensorBSpline<3,T> & final_volume, bool print_info) {

    // a, b determine knot vectors U, V
    gsKnotVector<T> kv1 = a.knots(0);
    gsKnotVector<T> kv2 = a.knots(1);
    gsKnotVector<T> kv3 = c.knots(1);

    int num1 = kv1.size() - kv1.degree() - 1;
    int num2 = kv2.size() - kv2.degree() - 1;
    int num3 = kv3.size() - kv3.degree() - 1;

    T alfa = 0.0;
    T beta = 0.0;
    T gamma = 0.0;
    if (print_info) {
        gsInfo << "\nDiscrete Coons patch 3D computation:\n";
        gsInfo << "=================================\n";
        gsInfo << "Dimensions of point net: " << num1 << ", " << num2 << ", " <<num3 << "\n";
    }

    gsMatrix<T> final_control_net(num1*num2*num3, 3);
    final_control_net.setZero();

    gsMatrix<T> a_cp = a.coefs();
    gsMatrix<T> b_cp = b.coefs();
    gsMatrix<T> c_cp = c.coefs();
    gsMatrix<T> d_cp = d.coefs();
    gsMatrix<T> e_cp = e.coefs();
    gsMatrix<T> f_cp = f.coefs();

    if ((a_cp.rows() != b_cp.rows()) || (c_cp.rows() != d_cp.rows()) || (e_cp.rows() != f_cp.rows())) {
        GISMO_ERROR("1: Incompatible input curves for computation of discrete Coons patch 3D! Exiting function ...");
        return 1;
    }

    T euclnorm1 = 0.0;
    gsVector<T> vec(3);
    gsMatrix<T> vec1(num1, 3);
    gsMatrix<T> vec2(num1, 3);
    gsMatrix<T> vec3(num1, 3);
    gsMatrix<T> vec4(num1, 3);

    for (index_t i = 0; i < num1; i++) {
        vec1.row(i) = a_cp.row(i) - c_cp.row(i); // a, c edge
        vec = vec1.row(i); euclnorm1 += euclideanNorm(vec);
        vec2.row(i) = a_cp.row((num2-1)*num1+i) - d_cp.row(i); // a, d edge
        vec = vec2.row(i); euclnorm1 += euclideanNorm(vec);
        vec3.row(i) = b_cp.row(i) - c_cp.row((num3-1)*num1+i); // b, c edge
        vec = vec3.row(i); euclnorm1 += euclideanNorm(vec);
        vec4.row(i) = b_cp.row((num2-1)*num1+i) - d_cp.row((num3-1)*num1+i); // b, d edge
        vec = vec4.row(i); euclnorm1 += euclideanNorm(vec);
    }

    T euclnorm2 = 0.0;
    gsMatrix<T> vec5(num2, 3);
    gsMatrix<T> vec6(num2, 3);
    gsMatrix<T> vec7(num2, 3);
    gsMatrix<T> vec8(num2, 3);

    for (index_t i = 0; i < num2; i++) {
        vec5.row(i) = a_cp.row(i*num1) - e_cp.row(i); // a, e edge
        vec = vec5.row(i); euclnorm2 += euclideanNorm(vec);
        vec6.row(i) = a_cp.row((i+1)*num1-1) - f_cp.row(i); // a, f edge
        vec = vec6.row(i); euclnorm2 += euclideanNorm(vec);
        vec7.row(i) = b_cp.row(i*num1) - e_cp.row((num3-1)*num2+i); // b, e edge
        vec = vec7.row(i); euclnorm2 += euclideanNorm(vec);
        vec8.row(i) = b_cp.row((i+1)*num1-1) - f_cp.row((num3-1)*num2+i); // b, f edge
        vec = vec8.row(i); euclnorm2 += euclideanNorm(vec);
    }

    T euclnorm3 = 0.0;
    gsMatrix<T> vec9(num3, 3);
    gsMatrix<T> vec10(num3, 3);
    gsMatrix<T> vec11(num3, 3);
    gsMatrix<T> vec12(num3, 3);

    for (index_t i = 0; i < num3; i++) {
        vec9.row(i) = c_cp.row(i*num1) - e_cp.row(i*num2); // c, e edge
        vec = vec9.row(i); euclnorm3 += euclideanNorm(vec);
        vec10.row(i) = c_cp.row((i+1)*num1-1) - f_cp.row(i*num2); // c, f edge
        vec = vec10.row(i); euclnorm3 += euclideanNorm(vec);
        vec11.row(i) = d_cp.row(i*num1) - e_cp.row((i+1)*num2-1); // d, e edge
        vec = vec11.row(i); euclnorm3 += euclideanNorm(vec);
        vec12.row(i) = d_cp.row((i+1)*num1-1) - f_cp.row((i+1)*num2-1); // d, f edge
        vec = vec12.row(i); euclnorm3 += euclideanNorm(vec);
    }

    if (print_info) {
        gsInfo << "Check of compotatibility conditions: " << euclnorm1 << " (U), " << euclnorm2 << " (V), " << euclnorm3 << " (W)\n";
    }

    if ((euclnorm1 > static_cast<T>(1e-8)) || (euclnorm2 > static_cast<T>(1e-8)) || (euclnorm3 > static_cast<T>(1e-8))) {
        GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }

    // Initialization - filling boundary surfaces control points to the final control net
    for (index_t j = 0; j < num2; j++) { // adding a, b
        for (index_t i = 0; i < num1; i++) {
            //gsInfo << j*num1+i << "\n";
            final_control_net.row(j*num1+i) = a_cp.row(j*num1+i);
            final_control_net.row((num3-1)*num1*num2+j*num1+i) = b_cp.row(j*num1+i);
        }
    }
    for (index_t k = 0; k < num3; k++) { // adding c, d
        for (index_t i = 0; i < num1; i++) {
            //gsInfo << k*num1*num2+i << "\n";
            final_control_net.row(k*num1*num2+i) = c_cp.row(k*num1+i);
            final_control_net.row(k*num1*num2+(num2-1)*num1+i) = d_cp.row(k*num1+i);
        }
    }
    for (index_t k = 0; k < num3; k++) { // adding e, f
        for (index_t j = 0; j < num2; j++) {
            //gsInfo << k*num1*num2+j*num1 << "\n";
            final_control_net.row(k*num1*num2+j*num1) = e_cp.row(k*num2+j);
            final_control_net.row(k*num1*num2+j*num1+num1-1) = f_cp.row(k*num2+j);
        }
    }
    //gsInfo << "Po inicializaci ...\n";
    //gsInfo << final_control_net << "\n";

    //gsInfo << num1*num2 << "\n";
    for (index_t k = 0; k < num3; k++) {
        gamma = (T) k/(num3-1);
        for (index_t j = 0; j < num2; j++) {
            beta = (T) j/(num2-1);
            for (index_t i = 0; i < num1; i++) {
                alfa = (T) i/(num1-1);
                final_control_net.row(k*num1*num2+j*num1+i) = (1-alfa) * final_control_net.row(k*num1*num2+j*num1) + alfa * final_control_net.row(k*num1*num2+j*num1+num1-1) +
                        (1-beta) * final_control_net.row(k*num1*num2+i) + beta * final_control_net.row(k*num1*num2+(num2-1)*num1+i) +
                        (1-gamma) * final_control_net.row(j*num1+i) + gamma * final_control_net.row((num3-1)*num1*num2+j*num1+i)
                        - ((1-alfa) * (1-beta) * final_control_net.row(k*num1*num2) + alfa * (1-beta) * final_control_net.row(k*num1*num2+num1-1) + (1-alfa) * beta * final_control_net.row(k*num1*num2+(num2-1)*num1) + alfa * beta * final_control_net.row(k*num1*num2+(num2-1)*num1+num1-1))
                        - ((1-beta) * (1-gamma) * final_control_net.row(i) + beta * (1-gamma) * final_control_net.row((num2-1)*num1+i) + (1-beta) * gamma * final_control_net.row((num3-1)*num1*num2+i) + beta * gamma * final_control_net.row((num3-1)*num1*num2+(num2-1)*num1+i))
                        - ((1-gamma) * (1-alfa) * final_control_net.row(j*num1) + gamma * (1-alfa) * final_control_net.row((num3-1)*num1*num2+j*num1) + (1-gamma) * alfa * final_control_net.row(j*num1+num1-1) + gamma * alfa * final_control_net.row((num3-1)*num1*num2+j*num1+num1-1))
                        + (1-gamma) * ((1-alfa) * (1-beta) * final_control_net.row(0) + alfa * (1-beta) * final_control_net.row(num1-1) + (1-alfa) * beta * final_control_net.row((num2-1)*num1) + alfa * beta * final_control_net.row((num2-1)*num1+num1-1))
                        + gamma * ((1-alfa) * (1-beta) * final_control_net.row((num3-1)*num1*num2) + alfa * (1-beta) * final_control_net.row((num3-1)*num1*num2+num1-1) + (1-alfa) * beta * final_control_net.row((num3-1)*num1*num2+(num2-1)*num1) + alfa * beta * final_control_net.row((num3-1)*num1*num2+(num2-1)*num1+num1-1));
            }
        }
    }
    //gsInfo << "Konecna ...\n";
    //gsInfo << final_control_net << "\n";

    /*gsMatrix<T> control_net_test(num1*num2,3);
    control_net_test.setZero();
    int k = 6;
    for (index_t j = 0; j < num2; j++) { // adding a, b
        for (index_t i = 0; i < num1; i++) {
            //gsInfo << j*num1+i << ", " << k*num1*num2+j*num1+i << "\n";
            control_net_test.row(j*num1+i) = final_control_net.row(k*num1*num2+j*num1+i);
        }
    }
    //gsInfo << control_net_test.rows() << "\n";
    //gsInfo << control_net_test << "\n";

    gsTensorBSpline<2,T> test_surface (kv1, kv2, control_net_test);
    gsWriteParaview(test_surface, "test_surface", 10000);
    gsMesh<T> test_surface_mesh;
    test_surface.controlNet(test_surface_mesh);
    gsWriteParaview( test_surface_mesh, "test_surface_mesh");*/

    gsTensorBSplineBasis<3, T> basis (kv1, kv2, kv3);
    final_volume = gsTensorBSpline<3, T> (basis, final_control_net);

    return 0;

}


template<class T>
int springModelPatch3D(gsTensorBSpline<2,T> a, gsTensorBSpline<2,T> b, gsTensorBSpline<2,T> c, gsTensorBSpline<2,T> d, gsTensorBSpline<2,T> e, gsTensorBSpline<2,T> f, gsTensorBSpline<3,T> & final_volume, bool print_info) {

    // a, b determine knot vectors U, V
    gsKnotVector<T> kv1 = a.knots(0);
    gsKnotVector<T> kv2 = a.knots(1);
    gsKnotVector<T> kv3 = c.knots(1);

    int num1 = kv1.size() - kv1.degree() - 1;
    int num2 = kv2.size() - kv2.degree() - 1;
    int num3 = kv3.size() - kv3.degree() - 1;

    if (print_info) {
        gsInfo << "\nSpring model patch 3D computation:\n";
        gsInfo << "=================================\n";
        gsInfo << "Dimensions of point net: " << num1 << ", " << num2 << ", " <<num3 << "\n";
    }

    gsMatrix<T> final_control_net(num1*num2*num3, 3);
    final_control_net.setZero();

    gsMatrix<T> a_cp = a.coefs();
    gsMatrix<T> b_cp = b.coefs();
    gsMatrix<T> c_cp = c.coefs();
    gsMatrix<T> d_cp = d.coefs();
    gsMatrix<T> e_cp = e.coefs();
    gsMatrix<T> f_cp = f.coefs();

    if ((a_cp.rows() != b_cp.rows()) || (c_cp.rows() != d_cp.rows()) || (e_cp.rows() != f_cp.rows())) {
        GISMO_ERROR("1: Incompatible input curves for computation of discrete Coons patch 3D! Exiting function ...");
        return 1;
    }

    T euclnorm1 = 0.0;
    gsVector<T> vec(3);
    gsMatrix<T> vec1(num1, 3);
    gsMatrix<T> vec2(num1, 3);
    gsMatrix<T> vec3(num1, 3);
    gsMatrix<T> vec4(num1, 3);

    for (index_t i = 0; i < num1; i++) {
        vec1.row(i) = a_cp.row(i) - c_cp.row(i); // a, c edge
        vec = vec1.row(i); euclnorm1 += euclideanNorm(vec);
        vec2.row(i) = a_cp.row((num2-1)*num1+i) - d_cp.row(i); // a, d edge
        vec = vec2.row(i); euclnorm1 += euclideanNorm(vec);
        vec3.row(i) = b_cp.row(i) - c_cp.row((num3-1)*num1+i); // b, c edge
        vec = vec3.row(i); euclnorm1 += euclideanNorm(vec);
        vec4.row(i) = b_cp.row((num2-1)*num1+i) - d_cp.row((num3-1)*num1+i); // b, d edge
        vec = vec4.row(i); euclnorm1 += euclideanNorm(vec);
    }

    T euclnorm2 = 0.0;
    gsMatrix<T> vec5(num2, 3);
    gsMatrix<T> vec6(num2, 3);
    gsMatrix<T> vec7(num2, 3);
    gsMatrix<T> vec8(num2, 3);

    for (index_t i = 0; i < num2; i++) {
        vec5.row(i) = a_cp.row(i*num1) - e_cp.row(i); // a, e edge
        vec = vec5.row(i); euclnorm2 += euclideanNorm(vec);
        vec6.row(i) = a_cp.row((i+1)*num1-1) - f_cp.row(i); // a, f edge
        vec = vec6.row(i); euclnorm2 += euclideanNorm(vec);
        vec7.row(i) = b_cp.row(i*num1) - e_cp.row((num3-1)*num2+i); // b, e edge
        vec = vec7.row(i); euclnorm2 += euclideanNorm(vec);
        vec8.row(i) = b_cp.row((i+1)*num1-1) - f_cp.row((num3-1)*num2+i); // b, f edge
        vec = vec8.row(i); euclnorm2 += euclideanNorm(vec);
    }

    T euclnorm3 = 0.0;
    gsMatrix<T> vec9(num3, 3);
    gsMatrix<T> vec10(num3, 3);
    gsMatrix<T> vec11(num3, 3);
    gsMatrix<T> vec12(num3, 3);

    for (index_t i = 0; i < num3; i++) {
        vec9.row(i) = c_cp.row(i*num1) - e_cp.row(i*num2); // c, e edge
        vec = vec9.row(i); euclnorm3 += euclideanNorm(vec);
        vec10.row(i) = c_cp.row((i+1)*num1-1) - f_cp.row(i*num2); // c, f edge
        vec = vec10.row(i); euclnorm3 += euclideanNorm(vec);
        vec11.row(i) = d_cp.row(i*num1) - e_cp.row((i+1)*num2-1); // d, e edge
        vec = vec11.row(i); euclnorm3 += euclideanNorm(vec);
        vec12.row(i) = d_cp.row((i+1)*num1-1) - f_cp.row((i+1)*num2-1); // d, f edge
        vec = vec12.row(i); euclnorm3 += euclideanNorm(vec);
    }

    if (print_info) {
        gsInfo << "Check of compotatibility conditions: " << euclnorm1 << " (U), " << euclnorm2 << " (V), " << euclnorm3 << " (W)\n";
    }

    if ((euclnorm1 > static_cast<T>(1e-8)) || (euclnorm2 > static_cast<T>(1e-8)) || (euclnorm3 > static_cast<T>(1e-8))) {
        GISMO_ERROR("2: Incompatible input curves for computation of discrete Coons patch! Exiting function ...");
        return 1;
    }

    // Initialization - filling boundary surfaces control points to the final control net
    for (index_t j = 0; j < num2; j++) { // adding a, b
        for (index_t i = 0; i < num1; i++) {
            //gsInfo << j*num1+i << "\n";
            final_control_net.row(j*num1+i) = a_cp.row(j*num1+i);
            final_control_net.row((num3-1)*num1*num2+j*num1+i) = b_cp.row(j*num1+i);
        }
    }
    for (index_t k = 0; k < num3; k++) { // adding c, d
        for (index_t i = 0; i < num1; i++) {
            //gsInfo << k*num1*num2+i << "\n";
            final_control_net.row(k*num1*num2+i) = c_cp.row(k*num1+i);
            final_control_net.row(k*num1*num2+(num2-1)*num1+i) = d_cp.row(k*num1+i);
        }
    }
    for (index_t k = 0; k < num3; k++) { // adding e, f
        for (index_t j = 0; j < num2; j++) {
            //gsInfo << k*num1*num2+j*num1 << "\n";
            final_control_net.row(k*num1*num2+j*num1) = e_cp.row(k*num2+j);
            final_control_net.row(k*num1*num2+j*num1+num1-1) = f_cp.row(k*num2+j);
        }
    }
    //gsInfo << "Po inicializaci ...\n";
    //gsInfo << final_control_net << "\n";

    // Setting up linear system for inner control points
    gsMatrix<T> m_A((num1-2)*(num2-2)*(num3-2), (num1-2)*(num2-2)*(num3-2));
    m_A.setZero();
    gsMatrix<T> m_B((num1-2)*(num2-2)*(num3-2), 3);
    m_B.setZero();
    int index;
    for (index_t k = 1; k < num3-1; k++) {
        for (index_t j = 1; j < num2-1; j++) {
            for (index_t i = 1; i < num1-1; i++) {
                //gsInfo << i << ", " << j << ", " << k << "\n";
                index = (k-1)*(num1-2)*(num2-2)+(j-1)*(num1-2)+i-1;
                m_A(index, index) = 6;
                if ((i-1) == 0) {
                    m_B.row(index) += final_control_net.row(k*num1*num2+j*num1);
                }
                else {
                    m_A(index, (k-1)*(num1-2)*(num2-2)+(j-1)*(num1-2)+i-2) = -1;
                }
                if ((i+1) == num1-1) {
                    m_B.row(index) += final_control_net.row(k*num1*num2+j*num1+num1-1);
                }
                else {
                    m_A(index, (k-1)*(num1-2)*(num2-2)+(j-1)*(num1-2)+i) = -1;
                }
                if ((j-1) == 0) {
                    m_B.row(index) += final_control_net.row(k*num1*num2+i);
                }
                else {
                    m_A(index, (k-1)*(num1-2)*(num2-2)+(j-2)*(num1-2)+i-1) = -1;
                }
                if ((j+1) == num2-1) {
                    m_B.row(index) += final_control_net.row(k*num1*num2+(num2-1)*num1+i);
                }
                else {
                    m_A(index, (k-1)*(num1-2)*(num2-2)+(j)*(num1-2)+i-1) = -1;
                }
                if ((k-1) == 0) {
                    m_B.row(index) += final_control_net.row(j*num1+i);
                }
                else {
                    m_A(index, (k-2)*(num1-2)*(num2-2)+(j-1)*(num1-2)+i-1) = -1;
                }
                if ((k+1) == num3-1) {
                    m_B.row(index) += final_control_net.row((num3-1)*num1*num2+j*num1+i);
                }
                else {
                    m_A(index, (k)*(num1-2)*(num2-2)+(j-1)*(num1-2)+i-1) = -1;
                }
            }
        }
    }

    //gsInfo << m_A << "\n";
    //gsInfo << m_B << "\n";

    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x=m_A.fullPivHouseholderQr().solve( m_B);

    //gsInfo << x << "\n";

    for (index_t k = 1; k < num3-1; k++) {
        for (index_t j = 1; j < num2-1; j++) {
            for (index_t i = 1; i < num1-1; i++) {
                final_control_net.row(k*num1*num2+j*num1+i) = x.row((k-1)*(num1-2)*(num2-2)+(j-1)*(num1-2)+i-1);
            }
        }
    }

    gsTensorBSplineBasis<3, T> basis (kv1, kv2, kv3);
    final_volume = gsTensorBSpline<3, T> (basis, final_control_net);

    return 0;

}

template<class T>
 gsMatrix<T> computeCircleArc2D(gsMatrix<T> start, T centreX, T centreY, T alpha)
{
     //constant for circle section

    T radius = math::sqrt(math::pow(start(0)-centreX,2)+math::pow(start(1)-centreY,2));
    T norm_koef =   radius * (4.0/3.0)*(math::tan(alpha/4));

    gsMatrix<T> end(2,1);
    end << centreX + (-centreX + start(0))*( math::cos(alpha)) - (-centreY + start(1)) * (math:: sin(alpha)),
           centreY + (-centreY + start(1))*(math::cos(alpha)) + (-centreX + start(0)) *(math::sin(alpha));

    gsMatrix<T> coefsclp(4, 2);
    coefsclp << start(0), start(1),
                start(0) - (norm_koef*(start(1)-centreY))/(math::sqrt(math::pow(start(0)-centreX,2)+math::pow(start(1)-centreY,2))), start(1) + (norm_koef*(start(0)-centreX))/(math::sqrt(math::pow(-start(0)+centreX,2)+math::pow(start(1)-centreY,2))),
                end(0) + (norm_koef*(end(1)-centreY))/(math::sqrt(math::pow(-end(1)+centreY,2)+math::pow(end(0)-centreX,2))), end(1) - (norm_koef*(end(0)-centreX))/(math::sqrt(math::pow(-end(1)+centreY,2)+math::pow(end(0)-centreX,2))),
                end(0), end(1);

//    gsInfo<<"\n--------------------------------------\n";
//    gsInfo<<"arc \n";
//    gsInfo<< coefsclp;
//    gsInfo<<"\n--------------------------------------\n";

    return coefsclp;
 }

 template<class T>
 gsMatrix<T> computeCircleArc3D(gsMatrix<T> start, gsMatrix<T> centre, T alpha, gsMatrix<T> axis)
{
     gsMatrix<T> circleArc3Dpom(4,2);
     gsMatrix<T> circleArc3D(4,3);
     gsMatrix<T> startpom(2,1);


     if (axis(0,0) == 1) //axis x
     {


         startpom(0,0) = start(1,0);
        startpom(1,0) = start(2,0);
        circleArc3Dpom = computeCircleArc2D(startpom,  centre(1,0), centre(2,0), alpha);

        for(int i = 0; i < 4; i++)
        {
             circleArc3D(i,0) = centre(0,0);
             circleArc3D(i,1) = circleArc3Dpom(i,0);
             circleArc3D(i,2) = circleArc3Dpom(i,1);

        }
     }

     if (axis(1,0) == 1) //axis y
     {
        startpom(0,0) = start(0,0);
        startpom(1,0) = start(2,0);
        circleArc3Dpom = computeCircleArc2D(startpom,  centre(0,0), centre(2,0), alpha);

        for(int i = 0; i < 4; i++)
        {
             circleArc3D(i,0) = circleArc3Dpom(i,0);
             circleArc3D(i,1) = centre(1,0);
             circleArc3D(i,2) = circleArc3Dpom(i,1);

        }
     }

     if (axis(2,0) == 1) //axis z
     {
        startpom(0) = start(0);
        startpom(1) = start(1);
        circleArc3Dpom = computeCircleArc2D(startpom,  centre(0,0), centre(1,0), alpha);

        for(int i = 0; i < 4; i++)
        {
             circleArc3D(i,0) = circleArc3Dpom(i,0);
             circleArc3D(i,1) = circleArc3Dpom(i,1);
             circleArc3D(i,2) = centre(2,0);

        }
     }

    return circleArc3D;
}

 template<class T>
 gsMultiPatch<T> makeLinear(gsMultiPatch<T> mp, int uR)
 {
     gsMultiPatch<T> mp_new;
     int dim = mp.dim();
     gsMatrix<int> degs(mp.nPatches(), dim);

     gsInfo << "Dimension: " << dim << "\n";

     for (index_t i = 0; i < mp.nPatches(); i++)
         for (int j = 0; j < dim; j++) {
             degs(i, j) = mp[i].degree(j);
         }
     gsInfo << "Degrees of patches of the multipatch:\n" << degs << "\n";

     for (index_t i = 0; i < mp.nPatches(); i++) {
         if (dim == 2) {
             gsTensorBSpline<2, T>* tbspline = dynamic_cast<gsTensorBSpline<2, T>*>(&(mp.patch(i)));
             /*for (int j = 0; j < 2; j++) {
                 tbspline->degreeReduce(degs(i,j)-1, j);
             }*/

             gsInfo << "Input B-spline:\n" << *tbspline << "\n";

             for (int j = 0; j < uR; j++) {
                 tbspline->uniformRefine();
             }

             std::vector<T> knots1 = (tbspline->knots(0)).breaks();
             std::vector<T> knots2 = (tbspline->knots(1)).breaks();
             gsMatrix<T> parpoints(2, (knots1.size()) * (knots2.size()));
             gsMatrix<T> points(2, (knots1.size()) * (knots2.size()));

             for (index_t b = 0; b < knots2.size(); b++)
                 for (index_t a = 0; a < knots1.size(); a++) {
                     parpoints(0, b * (knots1.size()) + a) = knots1[a];
                     parpoints(1, b * (knots1.size()) + a) = knots2[b];
                 }

             tbspline->eval_into(parpoints, points);

             gsWriteParaviewPoints<real_t>( points, "points_linear");

             gsKnotVector<T> kv1((tbspline->knots(0)).unique(), 1, 0);
             gsKnotVector<T> kv2((tbspline->knots(1)).unique(), 1, 0);

             gsInfo << kv1 << "\n";
             gsInfo << kv2 << "\n";
             gsInfo << points.rows() << ", " << points.cols() << "\n";

             gsTensorBSpline<2, T> tbspline_new(kv1, kv2, points.transpose());

             gsInfo << "Final B-spline:\n" << tbspline_new << "\n";

             mp_new.addPatch(tbspline_new);
         }
         else {
             gsTensorBSpline<3, T>* tbspline = dynamic_cast<gsTensorBSpline<3, T>*>(&(mp.patch(i)));
             /*for (int j = 0; j < dim; j++) {
                 tbspline->degreeReduce(degs(i,j)-1, j);
             }*/

             gsInfo << "Input B-spline:\n" << *tbspline << "\n";

             for (int j = 0; j < uR; j++) {
                 tbspline->uniformRefine();
             }

             std::vector<T> knots1 = (tbspline->knots(0)).breaks();
             std::vector<T> knots2 = (tbspline->knots(1)).breaks();
             std::vector<T> knots3 = (tbspline->knots(2)).breaks();
             gsMatrix<T> parpoints(3, (knots1.size()) * (knots2.size()) * (knots3.size()));
             gsMatrix<T> points(3, (knots1.size()) * (knots2.size()) * (knots3.size()));

             for (index_t c = 0; c < knots3.size(); c++)
                 for (index_t b = 0; b < knots2.size(); b++)
                     for (index_t a = 0; a < knots1.size(); a++) {
                         parpoints(0, c * (knots1.size()) * (knots2.size()) + b * (knots1.size()) + a) = knots1[a];
                         parpoints(1, c * (knots1.size()) * (knots2.size()) + b * (knots1.size()) + a) = knots2[b];
                         parpoints(2, c * (knots1.size()) * (knots2.size()) + b * (knots1.size()) + a) = knots3[c];
                     }

             tbspline->eval_into(parpoints, points);

             gsWriteParaviewPoints<real_t>( points, "points_linear");

             gsKnotVector<T> kv1((tbspline->knots(0)).unique(), 1, 0);
             gsKnotVector<T> kv2((tbspline->knots(1)).unique(), 1, 0);
             gsKnotVector<T> kv3((tbspline->knots(2)).unique(), 1, 0);

             gsInfo << kv1 << "\n";
             gsInfo << kv2 << "\n";
             gsInfo << kv3 << "\n";
             gsInfo << points.rows() << ", " << points.cols() << "\n";

             gsTensorBSpline<3, T> tbspline_new(kv1, kv2, kv3, points.transpose());

             gsInfo << "Final B-spline:\n" << tbspline_new << "\n";

             mp_new.addPatch(tbspline_new);
         }
     }

     for (index_t i = 0; i < mp.nInterfaces(); i++) {
         mp_new.addInterface(mp.bInterface(i));
     }

     mp_new.addAutoBoundaries();

     return mp_new;
 }



 template<class T>
 gsMultiPatch<T> makeLinearBasis(gsMultiPatch<T> mp, gsMultiBasis<T> tbasis, int uR)
 {
     int dim = mp.dim();
     gsMatrix<int> degs(mp.nPatches(), dim);


     gsMultiPatch<T> mp_new;

     for (size_t i = 0; i < mp.nPatches(); i++)
         for (int j = 0; j < dim; j++) {
             degs(i, j) = mp[i].degree(j);
         }

     gsInfo << "Dimension: " << dim << "\n";

     for (size_t i = 0; i < mp.nPatches(); i++)
         for (int j = 0; j < dim; j++) {
             degs(i, j) = mp[i].degree(j);
         }
     gsInfo << "Degrees of patches of the multipatch:\n" << degs << "\n";

     for (int j = 0; j < uR; j++) {
         tbasis.uniformRefine();
     }

     for (size_t i = 0; i < mp.nPatches(); i++) {
         if (dim == 2) {

             const gsTensorBSplineBasis<2, T>* tbsplinebasis = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&(tbasis.piece(i)));
             /*for (int j = 0; j < 2; j++) {
                 tbspline->degreeReduce(degs(i,j)-1, j);
             }*/
              gsTensorBSpline<2, T>* tbspline = dynamic_cast<gsTensorBSpline<2, T>*>(&(mp.patch(i)));

             gsKnotVector<T> kv1((tbsplinebasis->knots(0)).unique(), 1, 0);
             gsKnotVector<T> kv2((tbsplinebasis->knots(1)).unique(), 1, 0);

             std::vector<T> knots1 = (tbsplinebasis->knots(0)).breaks();
             std::vector<T> knots2 = (tbsplinebasis->knots(1)).breaks();
             gsMatrix<T> parpoints(2, (knots1.size()) * (knots2.size()));
             gsMatrix<T> points(2, (knots1.size()) * (knots2.size()));


             for (size_t b = 0; b < knots2.size(); b++)
                 for (size_t a = 0; a < knots1.size(); a++) {
                     parpoints(0, b * (knots1.size()) + a) = knots1[a];
                     parpoints(1, b * (knots1.size()) + a) = knots2[b];
                 }

             tbspline->eval_into(parpoints, points);

             gsTensorBSpline<2, T> tbspline_new(kv1, kv2, points.transpose());

             mp_new.addPatch(tbspline_new);


         }
         else {
             const gsTensorBSplineBasis<3, T>* tbsplinebasis = dynamic_cast<const gsTensorBSplineBasis<3, T>*>(&(tbasis.piece(i)));
             /*for (int j = 0; j < dim; j++) {
                 tbspline->degreeReduce(degs(i,j)-1, j);
             }*/
             gsTensorBSpline<3, T>* tbspline = dynamic_cast<gsTensorBSpline<3, T>*>(&(mp.patch(i)));

             gsKnotVector<T> kv1((tbsplinebasis->knots(0)).unique(), 1, 0);
             gsKnotVector<T> kv2((tbsplinebasis->knots(1)).unique(), 1, 0);
             gsKnotVector<T> kv3((tbsplinebasis->knots(2)).unique(), 1, 0);

             gsInfo << kv1 << "\n";
             gsInfo << kv2 << "\n";
             gsInfo << kv3 << "\n";

             std::vector<T> knots1 = (tbsplinebasis->knots(0)).breaks();
             std::vector<T> knots2 = (tbsplinebasis->knots(1)).breaks();
             std::vector<T> knots3 = (tbsplinebasis->knots(2)).breaks();
             gsMatrix<T> parpoints(3, (knots1.size()) * (knots2.size()) * (knots3.size()));
             gsMatrix<T> points(3, (knots1.size()) * (knots2.size()) * (knots3.size()));

             for (size_t c = 0; c < knots3.size(); c++)
                 for (size_t b = 0; b < knots2.size(); b++)
                     for (size_t a = 0; a < knots1.size(); a++) {
                         parpoints(0, c * (knots1.size()) * (knots2.size()) + b * (knots1.size()) + a) = knots1[a];
                         parpoints(1, c * (knots1.size()) * (knots2.size()) + b * (knots1.size()) + a) = knots2[b];
                         parpoints(2, c * (knots1.size()) * (knots2.size()) + b * (knots1.size()) + a) = knots3[c];
                     }

             tbspline->eval_into(parpoints, points);

             gsWriteParaviewPoints<real_t>( points, "points_linear");


             gsTensorBSpline<3, T> tbspline_new(kv1, kv2, kv3, points.transpose());

             gsInfo << "Final B-spline:\n" << tbspline_new << "\n";

             mp_new.addPatch(tbspline_new);
         }


     }


     for (size_t i = 0; i < mp.nInterfaces(); i++) {
         mp_new.addInterface(mp.bInterface(i));
     }

     mp_new.addAutoBoundaries();
     //gsMultiBasis<T> tbasisNew(mp_new);


     return mp_new;
 }


  template<class T>
 gsMatrix<T> LSQFergusonSmooth(gsMatrix<T> P0, gsMatrix<T> T0, gsMatrix<T> T1, gsMatrix<T> P3) {

     real_t alfa = (-3*P0(0,0)*(T0(0,0) - 4*T0(0,1)*T1(0,0)*T1(0,1) + 4*T0(0,0)*pow(T1(0,1),2)) + 3*P3(0,0)*(T0(0,0) - 4*T0(0,1)*T1(0,0)*T1(0,1) + 4*T0(0,0)*pow(T1(0,1),2)) + 3*(P0(0,1) - P3(0,1))*(4*T0(0,0)*T1(0,0)*T1(0,1) + T0(0,1)*(-5 + 4*pow(T1(0,1),2))))/
             (9 + 16*(-2*T0(0,0)*T0(0,1)*T1(0,0)*T1(0,1) + pow(T1(0,1),2) + pow(T0(0,1),2)*(1 - 2*pow(T1(0,1),2))));
     real_t beta = (3*(-(P0(0,0)*(T1(0,0) + 4*pow(T0(0,1),2)*T1(0,0) - 4*T0(0,0)*T0(0,1)*T1(0,1))) + P3(0,0)*(T1(0,0) + 4*pow(T0(0,1),2)*T1(0,0) - 4*T0(0,0)*T0(0,1)*T1(0,1)) + (P0(0,1) - P3(0,1))*(-5*T1(0,1) + 4*T0(0,1)*(T0(0,0)*T1(0,0) + T0(0,1)*T1(0,1)))))/
             (-9 + 32*T0(0,0)*T0(0,1)*T1(0,0)*T1(0,1) - 16*pow(T1(0,1),2) + 16*pow(T0(0,1),2)*(-1 + 2*pow(T1(0,1),2)));
     gsInfo << alfa << "\n";
     gsInfo << beta << "\n";

     gsMatrix<real_t> res(4,2);
     res << P0(0,0), P0(0,1),
            P0(0,0) + alfa * T0(0,0), P0(0,1) + alfa * T0(0,1),
            P3(0,0) + beta * T1(0,0), P3(0,1) + beta * T1(0,1),
            P3(0,0), P3(0,1);

     return res;
 }

   template<class T>
 gsMatrix<T> AxisDirectionNormed(gsMatrix<T> T0, gsMatrix<T> T1) {

     gsMatrix<real_t> aux = (T0/T0.norm() + T1/T1.norm());

     return aux/aux.norm();
 }

 template<class T>
 void UnifyKnotVectors(gsBSpline<T> &curve1, gsBSpline<T> curve2) {

     std::vector<real_t> knotsdiff;

     gsKnotVector<real_t> kv1 = curve1.knots();
     gsKnotVector<real_t> kv2 = curve2.knots();
     kv2.difference(kv1, knotsdiff);
     for (unsigned i = 0; i < knotsdiff.size(); i++) {
         curve1.insertKnot(knotsdiff[i]);
     }

 }

 template<class T>
 gsMatrix<T> LSQFergusonShort(gsMatrix<T> P0, gsMatrix<T> T0, gsMatrix<T> T1, gsMatrix<T> P3) {

     real_t alfa = (-(P0(0,0)*(T0(0,0) - T0(0,1)*T1(0,0)*T1(0,1) + T0(0,0)*pow(T1(0,1),2))) + P3(0,0)*(T0(0,0) - T0(0,1)*T1(0,0)*T1(0,1) + T0(0,0)*pow(T1(0,1),2)) + (P0(0,1) - P3(0,1))*(T0(0,0)*T1(0,0)*T1(0,1) + T0(0,1)*(-2 + pow(T1(0,1),2))))/
             (3 - 2*T0(0,0)*T0(0,1)*T1(0,0)*T1(0,1) + pow(T1(0,1),2) + pow(T0(0,1),2)*(1 - 2*pow(T1(0,1),2)));
     real_t beta = -((-(P0(0,0)*(1 + pow(T0(0,1),2))*T1(0,0)) + P0(0,0)*T0(0,0)*T0(0,1)*T1(0,1) + P3(0,0)*(T1(0,0) + pow(T0(0,1),2)*T1(0,0) - T0(0,0)*T0(0,1)*T1(0,1)) + (P0(0,1) - P3(0,1))*(T0(0,0)*T0(0,1)*T1(0,0) + (-2 + pow(T0(0,1),2))*T1(0,1)))/
                     (3 - 2*T0(0,0)*T0(0,1)*T1(0,0)*T1(0,1) + pow(T1(0,1),2) + pow(T0(0,1),2)*(1 - 2*pow(T1(0,1),2))));
     gsInfo << alfa << "\n";
     gsInfo << beta << "\n";
     alfa = abs(alfa);
     beta = abs(beta);

     gsMatrix<real_t> res(4,2);
     res << P0(0,0), P0(0,1),
            P0(0,0) + alfa * T0(0,0), P0(0,1) + alfa * T0(0,1),
            P3(0,0) + beta * T1(0,0), P3(0,1) + beta * T1(0,1),
            P3(0,0), P3(0,1);

     return res;
 }

 template<class T>
 std::vector<T> smartKnotIdentification(std::vector<T> kv_unique) {
     std::vector<real_t> kv_spans(kv_unique.size()-1);
     for (size_t i = 0; i < kv_unique.size()-1; i++)
         kv_spans[i] = kv_unique[i+1] - kv_unique[i];
     for (size_t i = 0; i < kv_spans.size(); i++)
         gsInfo << kv_spans[i] << ",";
     gsInfo << "\n";

     bool changed = false;
     gsVector<real_t> fr1(2), fr2(3);
     fr1 << 0.35, 0.65;
     fr2 << 0.25, 0.33, 0.42;
     std::vector<real_t> inserted_knots;
     index_t iter = 0;
     auto it1 = kv_spans.begin();
     do {
         changed = false;
         for (size_t i = 0; i < kv_spans.size()-1; i++) {
             gsInfo << i << "\n";
             if (kv_spans[i+1] <= fr1[1] * kv_spans[i]) {
                 gsInfo << "a\n";
                 inserted_knots.push_back(kv_unique[i] + fr2[0] * kv_spans[i]);
                 inserted_knots.push_back(kv_unique[i] + (fr2[0]+fr2[1]) * kv_spans[i]);
                 it1 = kv_spans.insert(it1 + i + 1, fr2[2] * kv_spans[i]);
                 kv_spans.insert(it1, fr2[1] * kv_spans[i]);
                 kv_spans[i] = fr2[0] * kv_spans[i];
                 it1 = kv_spans.begin();
                 for (size_t i = 0; i < kv_spans.size(); i++)
                     gsInfo << kv_spans[i] << ",";
                 gsInfo << "\n";
                 changed = true;
             }
             else if (kv_spans[i+1] < (0.99 * kv_spans[i])) {
                 gsInfo << "b\n";
                 inserted_knots.push_back(kv_unique[i] + fr1[0] * kv_spans[i]);
                 it1 = kv_spans.insert(it1 + i + 1, fr1[1] * kv_spans[i]);
                 kv_spans[i] = fr1[0] * kv_spans[i];
                 it1 = kv_spans.begin();
                 for (size_t i = 0; i < kv_spans.size(); i++)
                     gsInfo << kv_spans[i] << ",";
                 gsInfo << "\n";
                 changed = true;
             }
         }
         iter++;
         gsInfo << "iter = " << iter << "\n";
     } while (changed && (iter < 3));
     gsInfo << "inserted_knots: ";
     for (size_t i = 0; i < inserted_knots.size(); i++)
         gsInfo << inserted_knots[i] << ",";
     gsInfo << "\n";
     do {
         changed = false;
         for (size_t i = 0; i < kv_spans.size()-1; i++) {
             gsInfo << i << "\n";
             if ((3 * kv_spans[i]) < kv_spans[i+1]) {
                 gsInfo << "c\n";
                 inserted_knots.push_back(kv_unique[i+1] + fr1[0] * kv_spans[i+1]);
                 T span_old = kv_spans[i+1];
                 kv_spans[i+1] = fr1[1] * span_old;
                 it1 = kv_spans.insert(it1 + i + 1, fr1[0] * span_old);
                 it1 = kv_spans.begin();
                 for (size_t i = 0; i < kv_spans.size(); i++)
                     gsInfo << kv_spans[i] << ",";
                 gsInfo << "\n";
                 changed = true;
             }
         }
         iter++;
         gsInfo << "iter = " << iter << "\n";
     } while (changed);
     gsInfo << "inserted_knots: ";
     for (size_t i = 0; i < inserted_knots.size(); i++)
         gsInfo << inserted_knots[i] << ",";
     gsInfo << "\n";

     return inserted_knots;

 }

 template<class T>
 std::vector<T> smartKnotIdentification2(std::vector<T> kv_unique) {

     std::vector<real_t> inserted_knots;

     std::vector<real_t> kv_spans(kv_unique.size()-1);
     for (size_t i = 0; i < kv_unique.size()-1; i++)
         kv_spans[i] = kv_unique[i+1] - kv_unique[i];
     for (size_t i = 0; i < kv_spans.size(); i++)
         gsInfo << kv_spans[i] << ",";
     gsInfo << "\n";

     T min_span_size = 10;
     for (size_t i = 0; i < kv_spans.size(); i++)
         if (kv_spans[i] < min_span_size)
             min_span_size = kv_spans[i];
     min_span_size = 2 * min_span_size;
     gsInfo << " Min span size: " << min_span_size << "\n";

     auto it1 = kv_spans.begin();
     for (size_t i = 0; i < kv_spans.size(); i++) {
         if (0.5 * kv_spans[i] > min_span_size)
         {
             int ratio = kv_spans[i] / min_span_size;
             T new_span_length = kv_spans[i] / ratio;
             for (int j = 1; j < ratio; j++) {
                 inserted_knots.push_back(kv_unique[i] + j * new_span_length);
                 //kv_spans.insert(it1, fr2[1] * new_span_length);
             }
             //kv_spans[i] = new_span_length;
         }
         else
             it1++;
     }

     gsInfo << "inserted_knots: ";
     for (size_t i = 0; i < inserted_knots.size(); i++)
         gsInfo << inserted_knots[i] << ",";
     gsInfo << "\n";

     /*
     bool changed = false;
     gsVector<real_t> fr1(2), fr2(3);
     fr1 << 0.5, 0.5;
     fr2 << 0.33, 0.34, 0.33;
     std::vector<real_t> inserted_knots;
     index_t iter = 0;
     auto it1 = kv_spans.begin();
     do {
         changed = false;
         for (size_t i = 0; i < kv_spans.size()-1; i++) {
             gsInfo << i << "\n";
             if (kv_spans[i+1] <= fr1[1] * kv_spans[i]) {
                 gsInfo << "a\n";
                 inserted_knots.push_back(kv_unique[i] + fr2[0] * kv_spans[i]);
                 inserted_knots.push_back(kv_unique[i] + (fr2[0]+fr2[1]) * kv_spans[i]);
                 it1 = kv_spans.insert(it1 + i + 1, fr2[2] * kv_spans[i]);
                 kv_spans.insert(it1, fr2[1] * kv_spans[i]);
                 kv_spans[i] = fr2[0] * kv_spans[i];
                 it1 = kv_spans.begin();
                 for (size_t i = 0; i < kv_spans.size(); i++)
                     gsInfo << kv_spans[i] << ",";
                 gsInfo << "\n";
                 changed = true;
             }
             else if (kv_spans[i+1] < (0.99 * kv_spans[i])) {
                 gsInfo << "b\n";
                 inserted_knots.push_back(kv_unique[i] + fr1[0] * kv_spans[i]);
                 it1 = kv_spans.insert(it1 + i + 1, fr1[1] * kv_spans[i]);
                 kv_spans[i] = fr1[0] * kv_spans[i];
                 it1 = kv_spans.begin();
                 for (size_t i = 0; i < kv_spans.size(); i++)
                     gsInfo << kv_spans[i] << ",";
                 gsInfo << "\n";
                 changed = true;
             }
         }
         iter++;
         gsInfo << "iter = " << iter << "\n";
     } while (changed && (iter < 3));
     gsInfo << "inserted_knots: ";
     for (size_t i = 0; i < inserted_knots.size(); i++)
         gsInfo << inserted_knots[i] << ",";
     gsInfo << "\n";
     do {
         changed = false;
         for (size_t i = 0; i < kv_spans.size()-1; i++) {
             gsInfo << i << "\n";
             if ((3 * kv_spans[i]) < kv_spans[i+1]) {
                 gsInfo << "c\n";
                 inserted_knots.push_back(kv_unique[i+1] + fr1[0] * kv_spans[i+1]);
                 T span_old = kv_spans[i+1];
                 kv_spans[i+1] = fr1[1] * span_old;
                 it1 = kv_spans.insert(it1 + i + 1, fr1[0] * span_old);
                 it1 = kv_spans.begin();
                 for (size_t i = 0; i < kv_spans.size(); i++)
                     gsInfo << kv_spans[i] << ",";
                 gsInfo << "\n";
                 changed = true;
             }
         }
         iter++;
         gsInfo << "iter = " << iter << "\n";
     } while (changed);
     gsInfo << "inserted_knots: ";
     for (size_t i = 0; i < inserted_knots.size(); i++)
         gsInfo << inserted_knots[i] << ",";
     gsInfo << "\n";

     */

     return inserted_knots;

 }

 template<class T>
 T smartKnotParameterInsert(T param, gsKnotVector<T> kv) {
     T new_param;



     GISMO_ASSERT(kv.inDomain(param), "param is not in the parametric domain");

     std::vector<real_t>::const_iterator span;
     span = kv.iFind(param);

     T start_knot = kv.at(0);
     T start_end = kv.at(kv.size()-1);

     //insert param as the boundary of span or in the middle
     T middle_span = (span[0]+span[1])/2.;


     std::vector<T> diff_span(3);
     diff_span[0] = math::abs(span[0]-param);
     diff_span[1] = math::abs(span[1]-param);
     diff_span[2] = math::abs(middle_span-param);



     auto min_diff_span = std::min_element( diff_span.begin(), diff_span.end());
     bool case0 = ((span[0] - start_knot) < std::numeric_limits<T>::epsilon()) && (*min_diff_span == diff_span[0]);
     bool case1 = (math::abs(span[1] - start_end) < std::numeric_limits<T>::epsilon()) && (*min_diff_span == diff_span[1]);

     if (*min_diff_span == diff_span[2] || (case0 || case1))
     {
        new_param = middle_span;
     }
     else
     {
       new_param = span[std::distance(diff_span.begin(), min_diff_span)];
     }

     return new_param;

 }

#endif // UWBTURBINEUTILS
