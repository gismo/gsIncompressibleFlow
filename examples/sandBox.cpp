#include <gismo.h>
#include <bitset>

#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowFieldCreator.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>

using namespace gismo;

int main(int argc, char *argv[])
{ 
    int deg = 2;
    int numRefine = 2;
    int maxIt = 5;
    real_t viscosity = 0.1;
    real_t tol = 1e-5;

    real_t a = 8;
    real_t b = 2;
    real_t a_in = 1;
    gsMultiPatch<> patches = BSplineStep2D<real_t>(deg, a, b, a_in);
    patches.patch(2).coef(0,0) -= 0.3;
    patches.patch(2).coef(3,0) -= 0.1;

    gsBoundaryConditions<> bcInfo;
    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall; // containers of patch sides corresponding to inflow, outflow and wall boundaries
    gsFunctionExpr<> f("0", "0", 2); // external force

    gsMultiBasis<> basis(patches);
    defineBCs_step(bcInfo, bndIn, bndOut, bndWall, 2); 
    refineBasis_step(basis, numRefine, 0, 0, 0, 0, 2, a, b);

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases);
    gsINSSolverSteady<real_t> NSsolver(params);

    gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";
    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();
    NSsolver.solve(maxIt, tol);

    gsField<> velocity = NSsolver.constructSolution(0);
    gsField<> pressure = NSsolver.constructSolution(1);

    gsWriteParaview<>(velocity, "step_velocity", 30000, true);
    gsWriteParaview<>(pressure, "step_pressure", 30000);

    // -------------------------------
    // test of new boundary condition:

    real_t pTarget = 1.8;
    gsDiffScaledOuterNormalField<real_t> Uin(2, boundary::west, pressure, pTarget);
    gsBoundaryConditions<> bcInfo1;
    gsFunctionExpr<>Uwall("0", "0", 2);

    for (size_t i = 0; i < bndWall.size(); i++)
        bcInfo1.addCondition(bndWall[i].first, bndWall[i].second, condition_type::dirichlet, Uwall, 0);

    for (size_t i = 0; i < bndIn.size(); i++)
        bcInfo1.addCondition(bndIn[i].first, bndIn[i].second, condition_type::dirichlet, Uin, 0, true);

    gsNavStokesPde<real_t> NSpde1(patches, bcInfo1, &f, viscosity);
    gsFlowSolverParams<real_t> params1(NSpde1, discreteBases);
    gsINSSolverSteady<real_t> NSsolver1(params1);

    gsInfo << "numDofs: " << NSsolver1.numDofs() << "\n";
    gsInfo << "\ninitialization...\n";
    NSsolver1.initialize();
    NSsolver1.solve(1, tol);

    gsField<> velocity1 = NSsolver1.constructSolution(0);

    gsWriteParaview<>(velocity1, "step_velocity_newBC", 30000, true);

    // ----------------------------------------------------------------------------

    // gsTensorBSpline<2, real_t> geo = BSplineRectangle<real_t>(1, 0, 0, 1, 1);
    // geo.rotate(EIGEN_PI / 4);

    // gsWriteParaview<>(geo, "geo");

    // gsMultiPatch<> mp, field;
    // mp.addPatch(geo);

    // gsMatrix<> coeffs(geo.basis().size(), 1);
    // coeffs.setOnes();
    // coeffs *= 2*math::sqrt(2);
    // field.addPatch(geo.basis().makeGeometry(coeffs));

    // gsField<> pressure(mp, field, true);

    // gsMatrix<> result;

    // gsMatrix<> pts(1,2);
    // pts << 0.5, 0;
    // gsScaledOuterNormalField<real_t> PnormalField(0, boundary::south, pressure);
    // PnormalField.eval_into(pts.transpose(), result);
    // gsInfo << "\nresult1 =\n" << result << "\n";

    // gsDiffScaledOuterNormalField<real_t> PdiffNormalField(0, boundary::south, pressure, 3*math::sqrt(2));
    // PdiffNormalField.eval_into(pts.transpose(), result);
    // gsInfo << "\nresult2 =\n" << result << "\n";

    // ----------------------------------------------------------------------------

    // int wanted = 4;

    // gsMatrix<int> mat(3,2);
    // mat << 1, 2, 3, 4, 5, 6;

    // bool b = (mat.col(1).array() == wanted).any();
    // std::string str = (b == 1) ? "true" : "false";

    // gsInfo << "Is any element equal to " << wanted << "? - " << str << "\n";

    // ----------------------------------------------------------------------------

    // std::vector<gsMatrix<> > vec(3);
    // gsMatrix<> mat(2, 2);
    // mat << 1, 2, 3, 4;
    // vec[1] = mat;
    
    // gsInfo << "vec[0] = " << vec[0] << "\n";
    // gsInfo << "vec[1] = " << vec[1] << "\n";
    // gsInfo << "vec[2] = " << vec[2] << "\n";

    // ----------------------------------------------------------------------------

    // unsigned flag = NEED_VALUE;

    // gsInfo << "flag = " << flag << "\n";
    // gsInfo << "binary flag = " << std::bitset<10>(flag) << "\n";

    // flag = flag|NEED_MEASURE;

    // gsInfo << "flag = " << flag << "\n";
    // gsInfo << "binary flag = " << std::bitset<10>(flag) << "\n";

    // flag = flag|NEED_MEASURE;

    // gsInfo << "flag = " << flag << "\n";
    // gsInfo << "binary flag = " << std::bitset<10>(flag) << "\n";

    // unsigned flag1 = NEED_DERIV;

    // gsInfo << "flag & flag1 = " << (flag&flag1) << "\n";
    // gsInfo << "binary flag & flag1 = " << std::bitset<10>(flag&flag1) << "\n";

    // ----------------------------------------------------------------------------

    // gsKnotVector<> kv(0, 1, 3, 3);
    // gsTensorBSplineBasis<2> basis(kv, kv);

    // gsVector<real_t> point(2);
    // point << 0, 0; 

    // index_t basisID = 8;
    // gsMatrix<> supp = basis.support(basisID);
    // gsMatrix<index_t> elem = basis.elementSupport(basisID);
    // //size_t elemID = basis.elementIndex(point); // pada

    // gsInfo << "\nkv = " << kv << "\n";
    // gsInfo << "\nsupport " << basisID  << ":\n" << supp << "\n";
    // gsInfo << "\nelementSupport " << basisID  << ":\n" << elem << "\n";
    // //gsInfo << "elementIndex( " << point  << " ):\n" << elemID << "\n";

    // typename gsBasis<>::domainIter domIt = basis.makeDomainIterator(boundary::none);

    // while(domIt->good())
    // {
    //     bool inSupport = true; 

    //     for (index_t d = 0; d < 2; d++)
    //     {
    //         if ( (domIt->lowerCorner()[d] < supp(d,0)) ||  (domIt->upperCorner()[d] > supp(d,1)))
    //         {
    //             inSupport = false;
    //             break;
    //         }
    //     }

    //     gsInfo << "\nElement " << domIt->id() << " is in support of basis fcn " << basisID << ": " << inSupport;

    //     domIt->next();
    // }

    // gsInfo << "\n\n" << supp.middleCols(1,1);

}