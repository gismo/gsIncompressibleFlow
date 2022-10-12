#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "uwbReadWriteOpt.h"

using namespace gismo;

int main(int argc, char *argv[]){

    int numRef, addRefV, numRefU, numRefV, numKnotU, numKnotV, maxIt, maxPicardIt, linMaxIt, precmaxIt, precfill;
    real_t visc, omega, maxUin, tol, tolPicard, linTol, prectol, precdropTol, precgamma;
    bool plot, preciter, preclumpingM, preclumpingA;
    std::string precType, mainDir, matrixDir, inletBCDir, gammaOptDir, precondTestDir, precondTestTolDir;
    std::map <std::string, gsVector<std::string> > paramList;
    gsOptionList optBCin, optMatrixGen, optPTest, optPTestTol, optGammaOpt;

    std::string inFilePath = MOTOR_DATA_DIR "/uwb-pilsen/linearTesting.txt";
    paramList = readInputFile(inFilePath);
    get_parameter(paramList,mainDir,"mainDir");
    get_parameter(paramList,matrixDir,"matrixDir");
    get_parameter(paramList,inletBCDir,"inletBCDir");
    get_parameter(paramList,gammaOptDir,"gammaOptDir");
    get_parameter(paramList,precondTestDir,"precondTestDir");
    get_parameter(paramList,precondTestTolDir,"precondTestTolDir");

    gsFileManager::mkdir(mainDir);
    gsFileManager::mkdir(mainDir + matrixDir);
    gsFileManager::mkdir(mainDir + inletBCDir);
    gsFileManager::mkdir(mainDir + precondTestDir);
    gsFileManager::mkdir(mainDir + precondTestTolDir);
    gsFileManager::mkdir(mainDir + gammaOptDir);

    //==================================== BC ===============================================

    get_parameter(paramList,numRef,"numRef"); std::string st_numRef = util::to_string(numRef);
    get_parameter(paramList,addRefV,"addRefV"); std::string st_addRefV = util::to_string(addRefV);
    get_parameter(paramList,numRefU,"numRefU"); std::string st_numRefU = util::to_string(numRefU);
    get_parameter(paramList,numRefV,"numRefV"); std::string st_numRefV = util::to_string(numRefV);
    get_parameter(paramList,visc,"visc"); std::string st_visc = util::to_string(visc); st_visc.replace(1,1,"-");

    std::string addID = "_ref" + st_numRef + st_addRefV + st_numRefU + st_numRefV;

    get_parameter(paramList,numKnotU,"numKnotU");
    get_parameter(paramList,numKnotV,"numKnotV");
    get_parameter(paramList,maxIt,"maxIt");
    get_parameter(paramList,maxPicardIt,"maxPicardIt");
    get_parameter(paramList,linMaxIt,"linMaxIt");
    get_parameter(paramList,maxUin,"maxUin");
    get_parameter(paramList,tol,"tol");
    get_parameter(paramList,tolPicard,"tolPicard");
    get_parameter(paramList,linTol,"linTol");
    get_parameter(paramList,plot,"plot");
    get_parameter(paramList,precType,"precType");
    get_parameter(paramList,precmaxIt,"precmaxIt");
    get_parameter(paramList,precfill,"precfill");
    get_parameter(paramList,prectol,"prectol");
    get_parameter(paramList,precdropTol,"precdropTol");
    get_parameter(paramList,precgamma,"precgamma");
    get_parameter(paramList,preciter,"preciter");
    get_parameter(paramList,preclumpingM,"preclumpingM");
    get_parameter(paramList,preclumpingA,"preclumpingA");

    optBCin.addInt("numRef","Number of uniform refinements",numRef);
    optBCin.addInt("addRefV","Additional uniform refinements in v-direction",addRefV);
    optBCin.addInt("numRefU","Number of local refinements in u-direction (near blade)",numRefU);
    optBCin.addInt("numRefV","Number of local refinements in v-direction (near turbine case)",numRefV);
    optBCin.addReal("visc","viscosity",visc);
    optBCin.addString("addID","String to add to the output file name",addID);
    optBCin.addInt("numKnotU","Number of knot spans to be refined in u-direction",numKnotU);
    optBCin.addInt("numKnotV","Number of knot spans to be refined in v-direction",numKnotV);
    optBCin.addInt("maxIt","Maximum number of iterations",maxIt);
    optBCin.addInt("maxPicardIt","Maximum number of Picard iterations",maxPicardIt);
    optBCin.addInt("linMaxIt","Maximum number of iterations for solving linear systems",linMaxIt);
    optBCin.addReal("maxUin","maximum of the inlet velocity magnitude",maxUin);
    optBCin.addReal("tol","Stopping tolerance",tol);
    optBCin.addReal("tolPicard","Stopping tolerance for Picard",tolPicard);
    optBCin.addReal("linTol","Stopping tolerance for linear solvers",linTol);
    optBCin.addSwitch("plot","Plot results",plot);
    optBCin.addString("outPath","Path to the output file",mainDir + inletBCDir);
    optBCin.addString("precType","Preconditioner name",precType);
    optBCin.addInt("prec.maxIt","Max number of iterations for inner solves",precmaxIt);
    optBCin.addInt("prec.fill","Fill factor for ILUT precond",precfill);
    optBCin.addReal("prec.tol","Stopping tolerance for inner solves",prectol);
    optBCin.addReal("prec.dropTol","Drop tolerance for ILUT precond",precdropTol);
    optBCin.addReal("prec.gamma","Parameter gamma for AL",precgamma);
    optBCin.addSwitch("prec.iter","Subsystems iteratively",preciter);
    optBCin.addSwitch("prec.lumpingM","Use lumped diagonal mass matrices in preconditioners",preclumpingM);
    optBCin.addSwitch("prec.lumpingA","Use lumped diagonal blockA in SIMPLE-type preconditioners",preclumpingA);

    //============================================== MatrixGen ==================================================

    int GnumRef, GnumRefU, numRefW, degElev, GmaxIt, GmaxPicardIt, GlinMaxIt, GnumKnotU, numKnotW, sliceKnot, numVanes, GprecmaxIt, Gprecfill;
    real_t timeStep, Gtol, GtolPicard, GlinTol, Gprecgamma, Gprectol, GprecdropTol;
    std::string inputIC, GprecType;
    bool steady, linDomain, periodic, Gplot, saveMat, computeIC, Gpreciter, GpreclumpingM, GpreclumpingA;

    get_parameter(paramList,GnumRef,"GnumRef");
    get_parameter(paramList,GnumRefU,"GnumRefU");
    get_parameter(paramList,numRefW,"numRefW"); std::string st_numRefW = util::to_string(numRefW);
    get_parameter(paramList,degElev,"degElev"); std::string st_degElev = util::to_string(degElev);
    get_parameter(paramList,omega,"omega"); std::string st_omega = util::to_string(omega);
    get_parameter(paramList,timeStep,"timeStep");
    get_parameter(paramList,steady,"steady");
    get_parameter(paramList,GmaxIt,"GmaxIt");
    get_parameter(paramList,GmaxPicardIt,"GmaxPicardIt");
    get_parameter(paramList,GlinMaxIt,"GlinMaxIt");
    get_parameter(paramList,linDomain,"linDomain");
    get_parameter(paramList,periodic,"periodic");
    get_parameter(paramList,Gplot,"Gplot");
    get_parameter(paramList,saveMat,"saveMat");
    get_parameter(paramList,GnumKnotU,"GnumKnotU");
    get_parameter(paramList,numKnotW,"numKnotW");
    get_parameter(paramList,sliceKnot,"sliceKnot");
    get_parameter(paramList,numVanes,"numVanes");
    get_parameter(paramList,Gtol,"Gtol");
    get_parameter(paramList,GtolPicard,"GtolPicard");
    get_parameter(paramList,GlinTol,"GlinTol");
    get_parameter(paramList,computeIC,"computeIC");
    get_parameter(paramList,inputIC,"inputIC");
    get_parameter(paramList,GprecType,"GprecType");
    get_parameter(paramList,GprecmaxIt,"GprecmaxIt");
    get_parameter(paramList,Gprecfill,"Gprecfill");
    get_parameter(paramList,Gprectol,"Gprectol");
    get_parameter(paramList,GprecdropTol,"GprecdropTol");
    get_parameter(paramList,Gprecgamma,"Gprecgamma");
    get_parameter(paramList,Gpreciter,"Gpreciter");
    get_parameter(paramList,GpreclumpingM,"GpreclumpingM");
    get_parameter(paramList,GpreclumpingA,"GpreclumpingA");

    std::string inletBC = mainDir + inletBCDir + "guideVanesSol_visc_" + st_visc + addID + ".xml";

    std::string mainFiles = "runner_visc_" + st_visc + "_om_" + st_omega + "_ref_" + st_numRef + "_" + st_numRefU + "_" + st_numRefW + "_elev_" + st_degElev;
    std::string matrixFile = "runner_visc_" + st_visc + "_om_" + st_omega + addID + "_elev_" + st_degElev;
    if (steady){
        mainFiles += "_st";
        matrixFile += "_st";
    } else {
        mainFiles += "_unst";
        matrixFile += "_unst";
    }

    if (linDomain){
        mainFiles += "_lin";
        matrixFile += "_lin";
    }

    if (!periodic){
        mainFiles += "_nonper";
        matrixFile += "_nonper";
    }

    std::string exampleDir = mainDir + mainFiles + "/";

    int retsys;

    gsFileManager::mkdir(exampleDir);

    optMatrixGen.addInt("numRef","Number of uniform refinements",GnumRef);
    optMatrixGen.addInt("numRefU","Number of local refinements in u-direction (near turbine case)",GnumRefU);
    optMatrixGen.addInt("numRefW","Number of local refinements in w-direction (near blade)",numRefW);
    optMatrixGen.addInt("degElev","Number of degree elevations",degElev);
    optMatrixGen.addReal("visc","viscosity",visc);
    optMatrixGen.addReal("timeStep","Time step",timeStep);
    optMatrixGen.addReal("omega","Omega",omega);
    optMatrixGen.addSwitch("steady","Steady computation",steady);
    optMatrixGen.addString("inletBC","Path to the xml file with solution from guide vanes",inletBC);
    optMatrixGen.addInt("maxIt","Maximum number of iterations",GmaxIt);
    optMatrixGen.addInt("maxPicardIt","Maximum number of Picard iterations",GmaxPicardIt);
    optMatrixGen.addInt("linMaxIt","Maximum number of iterations for solving linear systems",GlinMaxIt);
    optMatrixGen.addSwitch("linDomain","Use linearized domain",linDomain);
    optMatrixGen.addSwitch("periodic","Periodic domain",periodic);
    optMatrixGen.addSwitch("plot","Plot results",Gplot);
    optMatrixGen.addSwitch("saveMat","Save matrices",saveMat);
    optMatrixGen.addInt("numKnotU","Number of knot spans to be refined in u-direction",GnumKnotU);
    optMatrixGen.addInt("numKnotW","Number of knot spans to be refined in w-direction",numKnotW);
    optMatrixGen.addInt("sliceKnot","Index of the knot on the guide vanes patch where it intersects the runner wheel domain",sliceKnot);
    optMatrixGen.addInt("numVanes","Number of guide vanes",numVanes);
    optMatrixGen.addReal("tol","Stopping tolerance",Gtol);
    optMatrixGen.addReal("tolPicard","Stopping tolerance for Picard",GtolPicard);
    optMatrixGen.addReal("linTol","Stopping tolerance for linear solvers",GlinTol);
    optMatrixGen.addSwitch("computeIC","Compute the initial condition",computeIC);
    optMatrixGen.addString("inputIC","Path to the xml file with initial solution","");
    optMatrixGen.addString("outPath","Path to the output file",exampleDir);
    optMatrixGen.addString("matPath","Path to the saved matrices",mainDir + matrixDir);
    optMatrixGen.addString("precType","Preconditioner name",GprecType);
    optMatrixGen.addInt("prec.maxIt","Max number of iterations for inner solves",GprecmaxIt);
    optMatrixGen.addInt("prec.fill","Fill factor for ILUT precond",Gprecfill);
    optMatrixGen.addReal("prec.tol","Stopping tolerance for inner solves",Gprectol);
    optMatrixGen.addReal("prec.dropTol","Drop tolerance for ILUT precond",GprecdropTol);
    optMatrixGen.addReal("prec.gamma","Parameter gamma for AL",Gprecgamma);
    optMatrixGen.addSwitch("prec.iter","Subsystems iteratively",Gpreciter);
    optMatrixGen.addSwitch("prec.lumpingM","Use lumped diagonal mass matrices in preconditioners",GpreclumpingM);
    optMatrixGen.addSwitch("prec.lumpingA","Use lumped diagonal blockA in SIMPLE-type preconditioners",GpreclumpingA);

    //================================================ PrecondTest ================================

    int dim, PprecmaxIt, Pprecfill, itStep, PmaxIt;
    real_t Ptol, Pprectol, PprecdropTol, Pprecgamma, PprecalphaP;
    std::string solver;
    bool lu, Ppreciter, PpreclumpingM, PpreclumpingA;
    gsVector<std::string> precN;

    std::string matMp = mainFiles + "_matPresMass.xml";
    std::string matMv = mainFiles + "_matVelMass.xml";
    std::string matNS = mainFiles + "_matNS.xml";

    get_parameter(paramList,solver,"solver");
    get_parameter(paramList,dim,"dim");
    get_parameter(paramList,itStep,"itStep");
    get_parameter(paramList,PmaxIt,"PmaxIt");
    get_parameter(paramList,Ptol,"Ptol");
    get_parameter(paramList,lu,"lu");
    get_parameter(paramList,precN,"precN");
    get_parameter(paramList,PprecmaxIt,"PprecmaxIt");
    get_parameter(paramList,Pprecfill,"Pprecfill");
    get_parameter(paramList,Pprectol,"Pprectol");
    get_parameter(paramList,PprecdropTol,"PprecdropTol");
    get_parameter(paramList,Pprecgamma,"Pprecgamma");
    get_parameter(paramList,PprecalphaP,"PprecalphaP");
    get_parameter(paramList,Ppreciter,"Ppreciter");
    get_parameter(paramList,PpreclumpingM,"PpreclumpingM");
    get_parameter(paramList,PpreclumpingA,"PpreclumpingA");

    optPTest.addString("matMp","Pressure mass matrix file name",matMp);
    optPTest.addString("matMv","Velocity mass matrix file name",matMv);
    optPTest.addString("matNS","NS matrix file name",matNS);
    optPTest.addString("matPath","Path to the the matrices (xml files)",mainDir + matrixDir);
    optPTest.addString("solver","Name of linear solver to be used",solver);
    optPTest.addInt("dim","Problem dimension",dim);
    optPTest.addInt("itStep","Number of iterations for printing res, err, time",itStep);
    optPTest.addInt("maxIt","Maximum number of iterations",PmaxIt);
    optPTest.addReal("visc","viscosity",visc);
    optPTest.addReal("tol","Stopping tolerance",Ptol);
    optPTest.addSwitch("lu","Compute direct solution with LU",lu);
    for (int i = 0; i < precN.size(); i++){
        optPTest.addString("precN." + util::to_string(i),"Names of preconditioners to be used",precN(i));
    }
    optPTest.addInt("precN.Size","Number of preconditioners",precN.size());
    optPTest.addInt("prec.maxIt","Max number of iterations for inner solves",PprecmaxIt);
    optPTest.addInt("prec.fill","Fill factor for ILUT precond",Pprecfill);
    optPTest.addReal("prec.tol","Stopping tolerance for inner solves",Pprectol);
    optPTest.addReal("prec.dropTol","Drop tolerance for ILUT precond",PprecdropTol);
    optPTest.addReal("prec.gamma","Parameter gamma for AL",Pprecgamma);
    optPTest.addReal("prec.alphaP","Pressure relaxation parameter for SIMPLE-type",PprecalphaP);
    optPTest.addSwitch("prec.iter","Subsystems iteratively",Ppreciter);
    optPTest.addSwitch("prec.lumpingM","Use lumped diagonal mass matrices in preconditioners",PpreclumpingM);
    optPTest.addSwitch("prec.lumpingA","Use lumped diagonal blockA in SIMPLE-type preconditioners",PpreclumpingA);

    std::string mainFilesTest = mainFiles + "_periodicTesting";
    std::string mainFilesTol = mainFiles;
    if (PpreclumpingA){
        mainFilesTest += "_lumping";
        mainFilesTol += "_lumping";
    }
    std::string ofPathTest = mainDir + precondTestDir + mainFilesTest + ".txt";

    optPTest.addString("ofPath","Path to the output file",ofPathTest);

    //========================================= PrecondTestTol====================================

    real_t precgammaM;

    get_parameter(paramList,precgammaM,"precgammaM");

    optPTestTol.addString("matMp","Pressure mass matrix file name",matMp);
    optPTestTol.addString("matMv","Velocity mass matrix file name",matMv);
    optPTestTol.addString("matNS","NS matrix file name",matNS);
    optPTestTol.addString("matPath","Path to the the matrices (xml files)",mainDir + matrixDir);
    optPTestTol.addString("solver","Name of linear solver to be used",solver);
    optPTestTol.addInt("dim","Problem dimension",dim);
    optPTestTol.addInt("itStep","Number of iterations for printing res, err, time",itStep);
    optPTestTol.addInt("maxIt","Maximum number of iterations",PmaxIt);
    optPTestTol.addReal("visc","viscosity",visc);
    optPTestTol.addReal("tol","Stopping tolerance",Ptol);
    optPTestTol.addSwitch("lu","Compute direct solution with LU",lu);
    for (int i = 0; i < precN.size(); i++){
        optPTestTol.addString("precN." + util::to_string(i),"Names of preconditioners to be used",precN(i));
    }
    optPTestTol.addInt("precN.Size","Number of preconditioners",precN.size());
    optPTestTol.addInt("prec.maxIt","Max number of iterations for inner solves",PprecmaxIt);
    optPTestTol.addInt("prec.fill","Fill factor for ILUT precond",Pprecfill);
    optPTestTol.addReal("prec.tol","Stopping tolerance for inner solves",Pprectol);
    optPTestTol.addReal("prec.dropTol","Drop tolerance for ILUT precond",PprecdropTol);
    optPTestTol.addReal("prec.gamma","Parameter gamma for AL",Pprecgamma);
    optPTestTol.addReal("prec.alphaP","Pressure relaxation parameter for SIMPLE-type",PprecalphaP);
    optPTestTol.addSwitch("prec.iter","Subsystems iteratively",Ppreciter);
    optPTestTol.addSwitch("prec.lumpingM","Use lumped diagonal mass matrices in preconditioners",PpreclumpingM);
    optPTestTol.addSwitch("prec.lumpingA","Use lumped diagonal blockA in SIMPLE-type preconditioners",PpreclumpingA);

    std::string ofPathTol = mainDir + precondTestTolDir + mainFilesTol + ".txt";

    optPTestTol.addString("ofPath","Path to the output file",ofPathTol);
    optPTestTol.addReal("prec.gammaM","Parameter gamma for MAL",precgammaM);

    //======================================= gammaOpt ===============================================

    int OmaxIt, OprecmaxIt, Oprecfill;
    real_t gammaMax, gammaMin, gammaStep, Otol, Oprectol, OprecdropTol;
    bool Opreciter, OpreclumpingM;
    gsVector<std::string> OprecN;

    get_parameter(paramList,OmaxIt,"OmaxIt");
    get_parameter(paramList,gammaMax,"gammaMax");
    get_parameter(paramList,gammaMin,"gammaMin");
    get_parameter(paramList,gammaStep,"gammaStep");
    get_parameter(paramList,Otol,"Otol");
    get_parameter(paramList,OprecN,"OprecN");
    get_parameter(paramList,OprecmaxIt,"OprecmaxIt");
    get_parameter(paramList,Oprecfill,"Oprecfill");
    get_parameter(paramList,Oprectol,"Oprectol");
    get_parameter(paramList,OprecdropTol,"OprecdropTol");
    get_parameter(paramList,Opreciter,"Opreciter");
    get_parameter(paramList,OpreclumpingM,"OpreclumpingM");

    optGammaOpt.addString("matMp","Pressure mass matrix file name",matMp);
    optGammaOpt.addString("matMv","Velocity mass matrix file name",matMv);
    optGammaOpt.addString("matNS","NS matrix file name",matNS);
    optGammaOpt.addString("matPath","Path to the the matrices (xml files)",mainDir + matrixDir);
    optGammaOpt.addString("solver","Name of linear solver to be used",solver);
    optGammaOpt.addInt("dim","Problem dimension",dim);
    optGammaOpt.addInt("maxIt","Maximum number of iterations",OmaxIt);
    optGammaOpt.addReal("gammaMax","Maximum gamma",gammaMax);
    optGammaOpt.addReal("gammaMin","Minimum gamma",gammaMin);
    optGammaOpt.addReal("gammaStep","Gamma increment",gammaStep);
    optGammaOpt.addReal("visc","viscosity",visc);
    optGammaOpt.addReal("tol","Stopping tolerance",Otol);
    for (int i = 0; i < OprecN.size(); i++){
        optGammaOpt.addString("precN." + util::to_string(i),"Names of preconditioners to be used",OprecN(i));
    }
    optGammaOpt.addInt("precN.Size","Number of preconditioners",OprecN.size());
    optGammaOpt.addInt("prec.maxIt","Max number of iterations for inner solves",OprecmaxIt);
    optGammaOpt.addInt("prec.fill","Fill factor for ILUT precond",Oprecfill);
    optGammaOpt.addReal("prec.tol","Stopping tolerance for inner solves",Oprectol);
    optGammaOpt.addReal("prec.dropTol","Drop tolerance for ILUT precond",OprecdropTol);
    optGammaOpt.addSwitch("prec.iter","Subsystems iteratively",Opreciter);
    optGammaOpt.addSwitch("prec.lumpingM","Use lumped diagonal mass matrices in preconditioners",OpreclumpingM);

    std::string ofPathGamma = mainDir + gammaOptDir + mainFiles + ".txt";

    optGammaOpt.addString("ofPathGamma","Path to the output file",ofPathGamma);

    gsWrite(optBCin,mainDir + inletBCDir + "linSolveRunnerBC" + addID + "_visc_" + st_visc + ".xml");
    gsWrite(optMatrixGen,mainDir + matrixDir + "M" + matrixFile + ".xml");
    gsWrite(optPTest,mainDir + precondTestDir + "T" + mainFilesTest + ".xml");
    gsWrite(optPTestTol,mainDir + precondTestTolDir + "TT" + mainFilesTol + ".xml");
    gsWrite(optGammaOpt,mainDir + gammaOptDir + "G" + mainFiles + ".xml");

    std::ifstream inFile(inFilePath);
    std::ofstream cpFile(exampleDir + "linearTesting.txt");
    cpFile << inFile.rdbuf();

    return 0;
}
