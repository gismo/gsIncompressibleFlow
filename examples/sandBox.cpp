#include <gismo.h>
#include <bitset>
#include <unordered_set>
#include <unordered_map>

#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowFieldCreator.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>

using namespace gismo;

int main(int argc, char *argv[])
{ 
    real_t viscosity = 0.1;
    int numRefine = 4;
    int maxIter = 10;
    int nParts = 4;

    int plotPts = 10000;
    bool plotMesh = true;
    std::string outPath = "";

    gsCmdLine cmd(".");
    cmd.addInt("", "maxIt", "maximum Picard iterations", maxIter);
    cmd.addInt("s", "plotSamples", "", plotPts);
    cmd.addInt("r", "refine", "global uniform refine", numRefine);
    cmd.addInt("b", "nParts", "number of blades", nParts);
    cmd.addReal("v", "visc", "viscosity", viscosity);
    cmd.addSwitch("", "mesh", "plot the computational mesh", plotMesh);
    cmd.addString("", "outPath", "path to the output directory", outPath);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsFunctionExpr<> Uin("-y/sqrt(x^2+y^2)", "x/sqrt(x^2+y^2)", 2);
    gsFunctionExpr<> Uwall("0", "0", 2);
    gsFunctionExpr<> f("0", "0", 2);

    real_t phi = -(2. / nParts)*EIGEN_PI;
    real_t cos = math::cos(phi);
    real_t sin = math::sin(phi);
    gsMatrix<real_t> transformMatrix(2, 2);
    transformMatrix(0, 0) = cos;
    transformMatrix(0, 1) = -sin;
    transformMatrix(1, 0) = sin;
    transformMatrix(1, 1) = cos;

    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uin, 0);
    bcInfo.addPeriodic(0, boundary::west, 0, boundary::east, 2);
    bcInfo.setTransformMatrix(transformMatrix);

    // bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uwall, 0);
    // bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);

    gsKnotVector<> kv(0,1,0,3);
    gsMatrix<> C(9,2);
    gsMatrix<> W(9,1);
    real_t w = math::sqrt(2)/2;
    C <<    3, 0,
            3, 3,
            0, 3,
            2, 0,
            2, 2,
            0, 2,
            1, 0,
            1, 1,
            0, 1;
    W << 1, w, 1, 1, w, 1, 1, w, 1;
    gsTensorNurbs<2> mp(kv, kv, C, W);

    gsMultiBasis<> tbasis(mp); 

    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);
    discreteBases.push_back(tbasis);
    discreteBases[0].degreeElevate(1);

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.dirValues = dirichlet::interpolation;
    opt.intStrategy = iFace::glue;

    gsNavStokesPde<real_t> NSpde(mp, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases);
    params.options().setString("assemb.loop", "EbE");

    gsINSSolverSteady<real_t> solver(params);

    // std::vector<gsMatrix<index_t> > elimDof(1);
    // elimDof[0].setZero(1,1);
    // solver.markDofsAsEliminatedZeros(elimDof, 1);

    gsInfo << "numDofs: " << solver.numDofs() << "\n\n";
    gsInfo << "initialization...\n";

    solver.initialize();
    solver.solve(maxIter, 1e-4, 0);

    // gsInfo << "velocity basis size = " << discreteBases[0].size() << "\n";
    // gsInfo << "velocity free size = " << solver.getAssembler()->getMappers()[0].freeSize() << "\n";
    // gsInfo << "pressure basis size = " << discreteBases[1].size() << "\n";
    // gsInfo << "pressure free size = " << solver.getAssembler()->getMappers()[1].freeSize() << "\n";

    // index_t udofs, pdofs;
    // udofs = 4;
    // pdofs = 8;

    // gsInfo << "\n========\nx-velocity:\n\n";
    // gsInfo << "matrix diag =\n" << solver.getAssembler()->matrix().toDense().topLeftCorner(udofs, udofs) << "\n\n";
    // gsInfo << "matrix off-diag =\n" << solver.getAssembler()->matrix().toDense().block(0, udofs, udofs, udofs) << "\n\n";
    // gsInfo << "rhs = " << solver.getAssembler()->rhs().topRows(udofs).transpose() << "\n\n";
    // gsInfo << "solution = " << solver.getSolution().topRows(udofs).transpose() << "\n\n";

    // gsInfo << "========\ny-velocity:\n\n";
    // gsInfo << "matrix diag =\n" << solver.getAssembler()->matrix().toDense().block(udofs, udofs, udofs, udofs) << "\n\n";
    // gsInfo << "matrix off-diag =\n" << solver.getAssembler()->matrix().toDense().block(udofs, 0, udofs, udofs) << "\n\n";
    // gsInfo << "rhs = " << solver.getAssembler()->rhs().middleRows(udofs, udofs).transpose() << "\n\n";
    // gsInfo << "solution = " << solver.getSolution().middleRows(udofs, udofs).transpose() << "\n\n";

    // gsInfo << "========\npressure:\n\n";
    // gsInfo << "matrix B1t =\n" << solver.getAssembler()->matrix().toDense().block(0, 2*udofs, udofs, pdofs) << "\n\n";
    // gsInfo << "matrix B2t =\n" << solver.getAssembler()->matrix().toDense().block(udofs, 2*udofs, udofs, pdofs) << "\n\n";
    // gsInfo << "matrix B1 =\n" << solver.getAssembler()->matrix().toDense().block(2*udofs, 0, pdofs, udofs) << "\n\n";
    // gsInfo << "matrix B2 =\n" << solver.getAssembler()->matrix().toDense().block(2*udofs, udofs, pdofs, udofs) << "\n\n";
    // gsInfo << "diag block =\n" << solver.getAssembler()->matrix().toDense().block(2*udofs, 2*udofs, pdofs, pdofs) << "\n\n";
    // gsInfo << "rhs = " << solver.getAssembler()->rhs().bottomRows(pdofs).transpose() << "\n\n";
    // gsInfo << "solution = " << solver.getSolution().bottomRows(pdofs).transpose() << "\n\n";

    // gsInfo << "\nAssembly time:" << solver.getInitAssemblyTime() + solver.getAssemblyTime() << "\n";
    // gsInfo << "Solve time:" << solver.getSolveTime() << "\n";
    // gsInfo << "Solver setup time:" << solver.getSolverSetupTime() << "\n";

    gsField<> velocity = solver.constructSolution(0);
    gsField<> pressure = solver.constructSolution(1);
    gsWriteParaview<>(velocity, outPath + "velocity", plotPts, plotMesh);
    gsWriteParaview<>(pressure, outPath + "pressure", plotPts, false);

    // -------------------------------
    // -------------------------------

    // int numRef = 2;
    // int deg = 1;

    // gsCmdLine cmd("");
    // cmd.addInt("r", "ref", "", numRef);
    // cmd.addInt("d", "deg", "", deg);
    // try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // gsKnotVector<> kv(0, 1, 0, deg+1);
    // gsBSplineBasis<> basis(kv);

    // for (int i = 0; i < numRef; i++)
    //    basis.uniformRefine();

    // gsMatrix<> pts(1,4);
    // pts << 0.3, 0.4, 0.6, 0.7;

    // gsMatrix<int> actives, activesUnique;
    // gsMatrix<> basisVals, basisVals1;

    // basis.active_into(pts, actives);
    // basis.eval_into(pts, basisVals);

    // activesUnique = createVectorOfUniqueIndices(actives);

    // gsInfo << "actives\n" << actives << "\n\n";
    // gsInfo << "activesUnique\n" << activesUnique << "\n\n";

    // std::unordered_map<int, int> activesUnique_val_to_ID;
    // for (int i = 0; i < activesUnique.size(); ++i)
    //     activesUnique_val_to_ID[activesUnique(i)] = i;

    // basisVals1.setZero(activesUnique.rows(), pts.cols());

    // for (int i = 0; i < actives.rows(); i++)
    //     for (int j = 0; j < actives.cols(); j++)
    //         basisVals1(activesUnique_val_to_ID[actives(i, j)], j) = basisVals(i,j);

    // gsInfo << "basisVals\n" << basisVals << "\n\n";
    // gsInfo << "basisVals1\n" << basisVals1 << "\n\n";


    // -------------------------------
    // -------------------------------
    // eval_into vs. evalSingle_into:

    // int numRef = 1;
    // int deg = 2;
    // int dim = 2;

    // gsCmdLine cmd("");
    // cmd.addInt("r", "ref", "", numRef);
    // cmd.addInt("d", "deg", "", deg);
    // try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // gsKnotVector<> kv(0, 1, 0, deg+1);
    // gsTensorBSplineBasis<2, real_t> basis(kv, kv);

    // for (int i = 0; i < numRef; i++)
    //     basis.uniformRefine();

    // gsVector<index_t> numQuadNodes(dim); 
    // numQuadNodes.setConstant(basis.maxDegree());
    // gsQuadRule<> rule = gsGaussRule<>(numQuadNodes);

    // // MOTOR VERSION

    // typename gsBasis<>::domainIter domIt = basis.makeDomainIterator(boundary::none);
    // gsMatrix<> quNodes;
    // gsVector<> quWeights;
    // gsMatrix<int> actives;
    // gsMatrix<> funData;
    // gsStopwatch clock;

    // clock.restart();
    // while (domIt->good())
    // {
    //     rule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    //     basis.active_into(quNodes.col(0), actives);
    //     basis.eval_into(quNodes, funData);

    //     domIt->next();
    // }
    // real_t t1 = clock.stop();

    // domIt->reset();

    // clock.restart();
    // while (domIt->good())
    // {
    //     rule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    //     basis.active_into(quNodes.col(0), actives);
    //     basis.deriv_into(quNodes, funData);

    //     domIt->next();
    // }
    // real_t t2 = clock.stop();


    // // gsIF VERSION

    // domIt->reset();

    // clock.restart();
    // while (domIt->good())
    // {
    //     rule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    //     basis.active_into(quNodes.col(0), actives);

    //     int numAct = actives.rows();
    //     funData.setZero(numAct, quNodes.cols());

    //     gsMatrix<> tmpData;
    //     for(index_t i = 0; i < numAct; i++)
    //     {
    //         basis.evalSingle_into(actives(i), quNodes, tmpData);
    //         funData.row(i) = tmpData;
    //     }

    //     domIt->next();
    // }
    // real_t t3 = clock.stop();

    // domIt->reset();

    // clock.restart();
    // while (domIt->good())
    // {
    //     rule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    //     basis.active_into(quNodes.col(0), actives);

    //     int numAct = actives.rows();
    //     funData.setZero(dim*numAct, quNodes.cols());

    //     gsMatrix<> tmpData;
    //     for(index_t i = 0; i < numAct; i++)
    //     {
    //         basis.derivSingle_into(actives(i), quNodes, tmpData);
    //         funData.middleRows(dim*i, dim) = tmpData;
    //     }

    //     domIt->next();
    // }
    // real_t t4 = clock.stop();

    // gsInfo << "motor: eval time = " << t1 << "\n";
    // gsInfo << "gsIF: eval time = " << t3 << "\n\n";
    // gsInfo << "eval time ratio = " << (t3/t1) << "\n\n";

    // gsInfo << "motor: deriv time = " << t2 << "\n";
    // gsInfo << "gsIF: deriv time = " << t4 << "\n\n";
    // gsInfo << "deriv time ratio = " << (t4/t2) << "\n\n";


    // -------------------------------
    // -------------------------------
    // triplets test:

    // int dim = 2;
    // int numRefine = 0;
    // int numIncrease = 0;

    // gsCmdLine cmd("Comparison of computational times for classic matrix assembly and assembly from triplets");
    // cmd.addInt("d", "dim", "space dimension", dim);
    // cmd.addInt("r", "refine", "uniform refine", numRefine);
    // cmd.addInt("i", "degIncr", "degree increase", numIncrease);
    // try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // real_t a = 8;
    // real_t b = 2;
    // real_t a_in = 1;

    // gsMultiPatch<> patches;
    
    // switch(dim)
    // {
    //     case 2:
    //         gsInfo << "Domain: 2D step\n";
    //         patches = BSplineStep2D<real_t>(1, a, b, a_in);
    //         break;

    //     case 3:
    //         gsInfo << "Domain: 3D step\n";
    //         patches = BSplineStep3D<real_t>(1, a, b, b, a_in);
    //         break;

    //     default: GISMO_ERROR("Wrong dimension!");
    // }

    // gsInfo << patches;

    // gsBoundaryConditions<> bcInfo;
    // std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall;
    // defineBCs_step(bcInfo, bndIn, bndOut, bndWall, dim);

    // gsMultiBasis<> basis(patches);
    
    // for (int i = 0; i < numRefine; i++)
    //     basis.uniformRefine();

    // basis.degreeIncrease(numIncrease);

    // gsDofMapper mapper;
    // basis.getMapper(dirichlet::elimination, iFace::glue, bcInfo, mapper, 0);

    // int numDofs = mapper.freeSize();
    // gsInfo << "numDofs = " << numDofs << "\n";

    // gsStopwatch clock, clockIn;

    // // ---------------------------
    // // fill matrix as in visitors:

    // gsInfo << "----------\n\n";

    // gsInfo << "\nAssembling matrix as in visitors - ColMajor.\n";

    // real_t reserveT1, assembleT1;
    // real_t insertT1 = 0;

    // gsSparseMatrix<real_t, ColMajor> globalMat1(numDofs, numDofs);

    // int nnzPerCol = 1;
    // for (int i = 0; i < dim; i++)
    //     nnzPerCol *= 2 * basis.maxDegree(i) + 1;

    // clock.restart();
    // globalMat1.reserve(gsVector<index_t>::Constant(numDofs, nnzPerCol));
    // reserveT1 = clock.stop();

    // gsInfo << "\nallocated size = " << globalMat1.data().allocatedSize() << "\n";

    // clock.restart();
    // for(size_t p = 0; p < patches.nPatches(); p++)
    // {
    //     gsVector<index_t> numQuadNodes(dim); 
    //     numQuadNodes.setConstant(basis.piece(p).maxDegree()+1);
    //     gsGaussRule<> quRule(numQuadNodes);
    //     gsMatrix<> quNodes;
    //     gsVector<> quWeights;

    //     typename gsBasis<>::domainIter domIt = basis.piece(p).makeDomainIterator(boundary::none);

    //     while (domIt->good())
    //     {
    //         quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

    //         gsMatrix<index_t> actives;
    //         basis.piece(p).active_into(quNodes.col(0), actives);

    //         int numAct = actives.rows();

    //         gsMatrix<> localMat(numAct, numAct);
    //         localMat.setOnes();

    //         mapper.localToGlobal(actives, p, actives);

    //         for (index_t i = 0; i < numAct; ++i)
    //         {
    //             const index_t ii = actives(i);

    //             if (mapper.is_free_index(ii))
    //             {
    //                 for (index_t j = 0; j < numAct; ++j)
    //                 {
    //                     const index_t jj = actives(j);

    //                     if (mapper.is_free_index(jj))
    //                     {
    //                         clockIn.restart();
    //                         globalMat1.coeffRef(ii, jj) = localMat(i, j);
    //                         insertT1 += clockIn.stop();
    //                     }
    //                 }
    //             }
    //         }

    //         domIt->next();
    //     }
    // }
    // assembleT1 = clock.stop();

    // globalMat1.makeCompressed();

    // gsInfo << "nonzeros = " << globalMat1.nonZeros() << "\n\n";
    // gsInfo << "reserve time = " << reserveT1 << "\n";
    // gsInfo << "insert time = " << insertT1 << "\n\n";

    // gsInfo << "----------\n\n";

    // gsInfo << "Assembling matrix as in visitors - RowMajor.\n";

    // real_t reserveT2, assembleT2;
    // real_t insertT2 = 0;

    // gsSparseMatrix<real_t, RowMajor> globalMat2(numDofs, numDofs);

    // clock.restart();
    // globalMat2.reserve(gsVector<index_t>::Constant(numDofs, nnzPerCol));
    // reserveT2 = clock.stop();

    // gsInfo << "\nallocated size = " << globalMat2.data().allocatedSize() << "\n";

    // clock.restart();
    // for(size_t p = 0; p < patches.nPatches(); p++)
    // {
    //     gsVector<index_t> numQuadNodes(dim); 
    //     numQuadNodes.setConstant(basis.piece(p).maxDegree()+1);
    //     gsGaussRule<> quRule(numQuadNodes);
    //     gsMatrix<> quNodes;
    //     gsVector<> quWeights;

    //     typename gsBasis<>::domainIter domIt = basis.piece(p).makeDomainIterator(boundary::none);

    //     while (domIt->good())
    //     {
    //         quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

    //         gsMatrix<index_t> actives;
    //         basis.piece(p).active_into(quNodes.col(0), actives);

    //         int numAct = actives.rows();

    //         gsMatrix<> localMat(numAct, numAct);
    //         localMat.setOnes();

    //         mapper.localToGlobal(actives, p, actives);

    //         for (index_t i = 0; i < numAct; ++i)
    //         {
    //             const index_t ii = actives(i);

    //             if (mapper.is_free_index(ii))
    //             {
    //                 for (index_t j = 0; j < numAct; ++j)
    //                 {
    //                     const index_t jj = actives(j);

    //                     if (mapper.is_free_index(jj))
    //                     {
    //                         clockIn.restart();
    //                         globalMat2.coeffRef(ii, jj) = localMat(i, j);
    //                         insertT2 += clockIn.stop();
    //                     }
    //                 }
    //             }
    //         }

    //         domIt->next();
    //     }
    // }
    // assembleT2 = clock.stop();

    // globalMat2.makeCompressed();

    // gsInfo << "nonzeros = " << globalMat2.nonZeros() << "\n\n";
    // gsInfo << "reserve time = " << reserveT2 << "\n";
    // gsInfo << "insert time = " << insertT2 << "\n\n";

    // gsInfo << "----------\n\n";

    // // ---------------------------
    // // fill matrix using triplets:

    // gsInfo << "Assembling matrix using triplets.\n";

    // real_t reserveT3, assembleT3, createTripletsT3, createMatT3;

    // gsSparseMatrix<real_t> globalMat3(numDofs, numDofs);

    // std::vector<gsEigen::Triplet<real_t> > triplets;

    // clock.restart();
    // triplets.reserve(numDofs * nnzPerCol);
    // reserveT3 = clock.stop();

    // clock.restart();

    // clockIn.restart();
    // for(size_t p = 0; p < patches.nPatches(); p++)
    // {
    //     gsVector<index_t> numQuadNodes(dim); 
    //     numQuadNodes.setConstant(basis.piece(p).maxDegree()+1);
    //     gsGaussRule<> quRule(numQuadNodes);
    //     gsMatrix<> quNodes;
    //     gsVector<> quWeights;

    //     typename gsBasis<>::domainIter domIt = basis.piece(p).makeDomainIterator(boundary::none);

    //     while (domIt->good())
    //     {
    //         quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

    //         gsMatrix<index_t> actives;
    //         basis.piece(p).active_into(quNodes.col(0), actives);

    //         int numAct = actives.rows();

    //         gsMatrix<> localMat(numAct, numAct);
    //         localMat.setOnes();

    //         mapper.localToGlobal(actives, p, actives);

    //         for (index_t i = 0; i < numAct; ++i)
    //         {
    //             const index_t ii = actives(i);

    //             if (mapper.is_free_index(ii))
    //             {
    //                 for (index_t j = 0; j < numAct; ++j)
    //                 {
    //                     const index_t jj = actives(j);

    //                     if (mapper.is_free_index(jj))
    //                         triplets.push_back(gsEigen::Triplet<real_t>(ii, jj, localMat(i, j)));
    //                 }
    //             }
    //         }

    //         domIt->next();
    //     }
    // }
    // createTripletsT3 = clockIn.stop();

    // clockIn.restart();
    // globalMat3.setFromTriplets(triplets.begin(), triplets.end());
    // createMatT3 = clockIn.stop();

    // assembleT3 = clock.stop();

    // gsInfo << "num triplets = " << triplets.size() << "\n";
    // gsInfo << "nonzeros = " << globalMat3.nonZeros() << "\n\n";
    // gsInfo << "reserve time = " << reserveT3 << "\n";
    // gsInfo << "create triplets time = " << createTripletsT3 << "\n";
    // gsInfo << "setFromTriplets time = " << createMatT3 << "\n\n";

    // gsInfo << "==========\n\n";
    // gsInfo << "classic ColMajor total time = " << assembleT1 << "\n";
    // gsInfo << "classic RowMajor total time = " << assembleT2 << "\n";
    // gsInfo << "triplet total time = " << assembleT3 << "\n\n";
    
    // -------------------------------
    // -------------------------------
    // test of jacobian check:

    // int npts = 10;

    // int deg = 2;
    // int dim = 2;
    // int numRefine = 0;
    // int maxIt = 5;
    // real_t viscosity = 0.1;
    // real_t tol = 1e-5;

    // real_t a = 8;
    // real_t b = 2;
    // real_t a_in = 1;

    // gsMultiPatch<> patches;
    // gsFunctionExpr<> f;
    
    // switch(dim)
    // {
    //     case 2:
    //         patches = BSplineStep2D<real_t>(deg, a, b, a_in);
    //         // patches.patch(0).coef(1,1) += 1.1;
    //         f = gsFunctionExpr<>("0", "0", 2);
    //         break;

    //     case 3:
    //         patches = BSplineStep3D<real_t>(deg, a, b, b, a_in);
    //         f = gsFunctionExpr<>("0", "0", "0", 3);
    //         break;

    //     default: GISMO_ERROR("Wrong dimension!");
    // }

    // patches.uniformRefine();

    // gsMultiBasis<> basis(patches);
    // gsBoundaryConditions<> bcInfo;
    // std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall; // containers of patch sides corresponding to inflow, outflow and wall boundaries
    
    // defineBCs_step(bcInfo, bndIn, bndOut, bndWall, dim); 
    // refineBasis_step(basis, numRefine, 0, 0, 0, 0, dim, a, b, b);

    // std::vector< gsMultiBasis<> >  discreteBases;
    // discreteBases.push_back(basis); // basis for velocity
    // discreteBases.push_back(basis); // basis for pressure
    // discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)

    // gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    // gsFlowSolverParams<real_t> params(NSpde, discreteBases);
    // gsINSSolverSteady<real_t> NSsolver(params);

    // gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";
    // gsInfo << "\ninitialization...\n";
    // NSsolver.initialize();
    // //NSsolver.solve(maxIt, tol);

    // int check = NSsolver.checkGeoJacobian(npts);

    // gsWriteParaview<>(patches, "domain", 30000, true);

    // gsField<> velocity = NSsolver.constructSolution(0);
    // gsField<> pressure = NSsolver.constructSolution(1);

    // gsWriteParaview<>(velocity, "step_velocity", 30000, true);
    // gsWriteParaview<>(pressure, "step_pressure", 30000);

    // -------------------------------
    // test of new boundary condition:

    // real_t pTarget = 1.8;
    // gsDiffScaledOuterNormalField<real_t> Uin(2, boundary::west, pressure, pTarget);
    // gsBoundaryConditions<> bcInfo1;
    // gsFunctionExpr<>Uwall("0", "0", 2);

    // for (size_t i = 0; i < bndWall.size(); i++)
    //     bcInfo1.addCondition(bndWall[i].first, bndWall[i].second, condition_type::dirichlet, Uwall, 0);

    // for (size_t i = 0; i < bndIn.size(); i++)
    //     bcInfo1.addCondition(bndIn[i].first, bndIn[i].second, condition_type::dirichlet, Uin, 0, true);

    // gsNavStokesPde<real_t> NSpde1(patches, bcInfo1, &f, viscosity);
    // gsFlowSolverParams<real_t> params1(NSpde1, discreteBases);
    // gsINSSolverSteady<real_t> NSsolver1(params1);

    // gsInfo << "numDofs: " << NSsolver1.numDofs() << "\n";
    // gsInfo << "\ninitialization...\n";
    // NSsolver1.initialize();
    // NSsolver1.solve(1, tol);

    // gsField<> velocity1 = NSsolver1.constructSolution(0);

    // gsWriteParaview<>(velocity1, "step_velocity_newBC", 30000, true);

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