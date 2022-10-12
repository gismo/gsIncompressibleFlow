/*
    Gradient-free swarm method for optimization problem. Miminization of objective function.
    Author(s): J. Egermaier
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "uwbSwarm7.h"
#include "uwb2DRF_max7.h"
#include "uwbReadWriteOpt.h"

using namespace gismo;


/*
void write_input_data(std::ofstream & ofile, std::stringstream & filename, std::string inputfile)
{ 
    std::ifstream infile(inputfile);
    std::string content = "";
    int count;
    for(count=0 ; infile.eof()!=true ; count++){ content += infile.get();}
    count--;
    content.erase(content.end()-1);     // erase last character
    infile.close();
    ofile.open(filename.str(), std::ios_base::app);
    ofile << content;                 // output
    ofile.close();
}
*/

void writeToFile(std::ofstream & ofile, std::stringstream & filename, swarm<real_t> S, int i, gsVector<std::string> nameOfOpen)
{
    int sizenOO = nameOfOpen.size();
    ofile.open(filename.str(), std::ios_base::app);
    ofile << "member " << i << ": x = " << S.member[i].x.asRowVector() << "\n"
          << "f = "<< S.member[i].f << ", fp = " << S.member[i].get_bestValue() << ", penaltyF = " << S.member[i].penaltyF << "\n";
    for (int j = 0; j < sizenOO; j++){
        ofile << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << "\n" 
              << nameOfOpen(j) + ".f_eff(clear efficiency) = " << S.member[i].open[j].f_eff << ", " << nameOfOpen(j) + ".f(weighted) = " << S.member[i].open[j].f << "\n"
              << nameOfOpen(j) + ".f_parts = " << S.member[i].open[j].f_parts << "\n"
              << nameOfOpen(j) + ".N_iter = " << S.member[i].open[j].N_iter.asRowVector() << "\n" + nameOfOpen(j) + ".CH = " << S.member[i].open[j].CH.asRowVector() 
              << "\n" + nameOfOpen(j) + ".CH_f = " << S.member[i].open[j].CH_f.asRowVector() << "\n";
    }
    ofile << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << "\n";
    ofile.close();
}     
    
void writeToFileFg(std::ofstream & ffile, std::stringstream & filename, swarm<real_t> S, int i, gsVector<std::string> nameOfOpen)
{
    int sizenOO = nameOfOpen.size();
    ffile.open(filename.str(), std::ios_base::app);
    ffile << S.member[i].f; 
    for (int j = 0; j < sizenOO; j++){
        ffile << ", " << S.member[i].open[j].f_eff << ", " << S.member[i].open[j].f;
    }
    ffile << "\n";
    for (int j = 0; j < sizenOO; j++){
        ffile << S.member[i].open[j].f_parts << "\n";
    } 
    ffile.close();
}

void writeToFileG (std::ofstream & ofile, std::stringstream & filename, swarm<real_t> S)
{
    ofile.open(filename.str(), std::ios_base::app);
    gsVector<real_t> g = S.get_bestPosition();
    int bi = S.get_bestIndex();
    ofile << "fg = " << S.get_bestValue() << ", member " << bi << "\n" << "g = " << g.asRowVector() << "\n";
    ofile.close();
}
    
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

int main(int argc, char *argv[]){

    srand (time(NULL));

//================= initial parameters ===========================

    bool IpOpt, load_solution;
    int nM, nCycles, nSwarmCycles, nP, nDV, nthr, bi;
    gsVector<real_t> lowerBounds, upperBounds, x0, max_tol;
    std::map <std::string, gsVector<std::string>> paramList;
    std::string grid;
    gsVector<std::string> draw_solution;
    real_t chyba;
    real_t hCoef = 1.0; //penalty parameter
    uwbOptFlowSolver<real_t> optimizer = uwbOptFlowSolver<real_t>(); //optimization object

    get_initialValues(x0, lowerBounds, upperBounds, max_tol, nP, nDV);
    optimizer.nP = nP;
    int nDesignV = nP*nDV;

    paramList = readInputFile("inputData.txt");
    get_parameter(paramList,optimizer.allOnameOfOpen,"name_of_open");
    int nOO = optimizer.allOnameOfOpen.size();
    optimizer.numOfOpens = nOO;

    optimizer.initialize_vectors();
    get_parameter(paramList,nM,"size_of_population");
    get_parameter(paramList,nCycles,"number_of_cycles");
    get_parameter(paramList,nSwarmCycles,"number_of_generations");
    get_parameter(paramList,nthr,"number_of_threads");
    get_parameter(paramList,IpOpt,"IpOpt");
    get_parameter(paramList,optimizer.loftEfficiency,"lofted_efficiency");
    get_parameter(paramList,optimizer.weightLift,"weightLift");
    get_parameter(paramList,optimizer.weightVelocity,"weightVelocity");
    get_parameter(paramList,optimizer.weightPressure,"weightPressure");
    get_parameter(paramList,optimizer.weightEfficiency,"weightEfficiency");
    get_parameter(paramList,load_solution,"load_solution");
    get_parameter(paramList,draw_solution,"draw_solution");
    get_parameter(paramList,optimizer.viscosity,"viscosity");
    get_parameter(paramList,optimizer.turbIntensity,"turbulence_intensity");
    get_parameter(paramList,optimizer.timeStep,"time_step");
    get_parameter(paramList,optimizer.minNumIter,"minimum_number_of_iterations");
    get_parameter(paramList,optimizer.numIterSteadyNS,"number_of_iterations_steadyNS");
    get_parameter(paramList,optimizer.numIterKOmegaSteady,"number_of_iterations_kOmegaSteady");
    get_parameter(paramList,optimizer.numIterFirst,"number_of_iterations_start");
    get_parameter(paramList,optimizer.numIterAll,"number_of_iterations_optimization");
    get_parameter(paramList,optimizer.numIterGrad,"number_of_iterations_gradient");  
    get_parameter(paramList,optimizer.velAVG,"velocity_by_average");  
    get_parameter(paramList,optimizer.rr,"radius_of_profiles");                                           check_vector_length(optimizer.rr, nP);
    get_parameter(paramList,optimizer.tolSteady,"tolerance_relative_error_of_steadyNS_solution");
    get_parameter(paramList,optimizer.tolRelNorm,"tolerance_relative_error_of_RANS_solution");
    get_parameter(paramList,optimizer.tolObjRelVal,"tolerance_relative_error_of_objective_function");
    get_parameter(paramList,optimizer.tmEvaluator,"turbulence_model");
    get_parameter(paramList,grid,"grid");

    if (grid == "7000Dofs"){
        optimizer.refinement = 1;
        optimizer.numRefine = 1;
        optimizer.numRefineUniformLocal_v = 0; optimizer.numUniformKnot_v = 1;
        optimizer.numRefineLocal_v = 4;        optimizer.numKnot_v = 1;
        optimizer.numRefineLocal_u0 = 4;       optimizer.numKnot_u0 = 1;
        optimizer.numRefineLocal_u1_s = 4;     optimizer.numKnot_u1_s = 1;
        optimizer.numRefineLocal_u1_e = 0;     optimizer.numKnot_u1_e = 1;
        optimizer.numRefineLocal_u2 = 0;       optimizer.numKnot_u2 = 1;
        optimizer.numRefineLocalFirstKnot_v = 0;
    } else if (grid == "17000Dofs"){
        optimizer.refinement = 2;
        optimizer.numRefine = 1;
        optimizer.numRefineUniformLocal_v = 0; optimizer.numUniformKnot_v = 1;
        optimizer.numRefineLocal_v = 4;        optimizer.numKnot_v = 2;
        optimizer.numRefineLocal_u0 = 4;       optimizer.numKnot_u0 = 1;
        optimizer.numRefineLocal_u1_s = 4;     optimizer.numKnot_u1_s = 1;
        optimizer.numRefineLocal_u1_e = 0;     optimizer.numKnot_u1_e = 1;
        optimizer.numRefineLocal_u2 = 0;       optimizer.numKnot_u2 = 1;
        optimizer.numRefineLocalFirstKnot_v = 5;
    } else {
        gsInfo << "Wrong choice of the grid!" << "\n";
    }
    for (int i = 0; i < nOO; i++){
        get_parameter(paramList,optimizer.allOdiffToOptOpen(i),"angle_difference_to_optimal_open_" + optimizer.allOnameOfOpen(i));
        get_parameter(paramList,optimizer.allOviscosityRatio(i),"viscosity_ratio_" + optimizer.allOnameOfOpen(i));        check_vector_length(optimizer.allOviscosityRatio(i), nP);
        get_parameter(paramList,optimizer.allOflowRate(i),"flow_rate_" + optimizer.allOnameOfOpen(i));
        get_parameter(paramList,optimizer.allOvelocityRelativeX(i),"input_velocity_x_" + optimizer.allOnameOfOpen(i));    check_vector_length(optimizer.allOvelocityRelativeX(i), nP);
        get_parameter(paramList,optimizer.allOvelocityRelativeY(i),"input_velocity_y_" + optimizer.allOnameOfOpen(i));    check_vector_length(optimizer.allOvelocityRelativeY(i), nP);
        get_parameter(paramList,optimizer.allOvOutTarget(i),"output_velocity_tangential_" + optimizer.allOnameOfOpen(i)); check_vector_length(optimizer.allOvOutTarget(i), nP);
        get_parameter(paramList,optimizer.allOuniformityParam(i),"uniformity_parameter_" + optimizer.allOnameOfOpen(i)); 
    }

//=====================================================================================================
gsInfo << "Draw_solution = " << draw_solution(0) << ", velikosti " << draw_solution.size() << "\n";
    if (draw_solution(0) != ""){
        for (int i = 0; i < draw_solution.size(); i++){
            gsFileData<real_t> xRead(draw_solution(i) + "_profile_designParameters.xml");
            gsMatrix<real_t> loadX = *(xRead.getFirst< gsMatrix<real_t> >());
            gsVector<real_t> x_draw = loadX.col(0);
            optimizer.flowSolver_load(1, draw_solution(i));//start with precomputed solution
            swarmMember<real_t> drawMember = swarmMember<real_t>(x_draw, nP, nOO);
            for (int j = 0; j < nOO; j++){ 
                drawMember.open[j].previousSolution = optimizer.allOpreviousSolutionBest(j);
            }
            optimizer.flowSolver_draw(drawMember, draw_solution(i));
        }
    }

else { // =========================== !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ======================================================
    
    std::stringstream fileHybrid; fileHybrid << "hybrid.txt"; std::ofstream ofile;
    std::stringstream fileSummary; fileSummary << "summary.txt"; std::ofstream sfile;
    std::stringstream fileFg; fileFg << "fg.txt"; std::ofstream ffile;
    ofile.open(fileHybrid.str());sfile.open(fileSummary.str());ffile.open(fileFg.str());
    ofile.close();sfile.close();ffile.close();

    write_input_data(ofile, fileHybrid, "inputData.txt"); 
    write_input_data(sfile, fileSummary, "inputData.txt");

//================= start solution - computed here ====================================
    if (!load_solution){
        swarmMember<real_t> startMember = swarmMember<real_t>(x0, nP, nOO);
        real_t nothing = optimizer.flowSolver_swarm(nM, startMember, true); 

        ofile.open(fileHybrid.str(), std::ios_base::app);
        for (int j = 0; j < nOO; j++){
            ofile << "Iniial flow solution =============================" 
                  << optimizer.allOnameOfOpen(j) + ".N_iter = " << startMember.open[j].N_iter.asRowVector() << "\n" 
                  << optimizer.allOnameOfOpen(j) + ".CH = " << startMember.open[j].CH.asRowVector() << "\n"
                  << optimizer.allOnameOfOpen(j) + ".CH_f = " << startMember.open[j].CH_f.asRowVector() << "\n"
                  << "==================================================" << "\n";
        }
        ofile.close();    
        //only for previouSolution
        optimizer.flowSolver_save(x0,"init",true);
    }
//================= start solution - load from file ===================================

    optimizer.flowSolver_load(nM, "init");//start with precomputed solution

    ofile.open(fileHybrid.str(), std::ios_base::app);
    for (int j = 0; j < nOO; j++){
        ofile << optimizer.allOnameOfOpen(j) + ".pTarget = " << optimizer.allOpTarget(j) << "\n";
        ofile << optimizer.allOnameOfOpen(j) + ".initObjectiveParts = " << optimizer.allOinitObjParts(j) << "\n";
    }
    ofile.close();

    gsVector<real_t> g (nDesignV); g.setZero(nDesignV); //best design variables

    //================================= optimization cycles ================================

    for (int cycle = 0; cycle < nCycles; cycle++){
        swarm<real_t> S = swarm<real_t>(); // swarm constructor 
        S.initBounds(nDesignV); // initialize both bounds to +-infinity

        for(int i = 0; i < nDesignV; i++){ // set bounds with respect to x0 values 
            S.set_lowerBound(i,std::max(x0[i]*(1.0 - max_tol(i)*(1.0 - (double)cycle/(double)nCycles)),lowerBounds[i])); //decrease starting interval
            S.set_upperBound(i,std::min(x0[i]*(1.0 + max_tol(i)*(1.0 - (double)cycle/(double)nCycles)),upperBounds[i]));
        }
        S.trialMembers(nM, x0, 0.3, nP, nOO); // random design variables with suitable penalty function, x0 and build swarm;
        //S.bestMember();

        for (int i = 0; i < nM; i++){
            for (int j = 0; j < nOO; j++){ 
                S.member[i].open[j].previousSolution = optimizer.allOpreviousSolutionBest(j);
            }
        }

        //========================= swarm Iterations ======================================

        for(int iter = 0; iter < nSwarmCycles; iter++){
            if (iter > 0){S.next_gen();} // new generation of the swarm, actualize fp

            //omp_set_nested(1);
            #pragma omp parallel for schedule(dynamic) num_threads(nthr)
            for(int i = 0; i < nM; i++){ 
                S.member[i].f = optimizer.flowSolver_swarm(nM, S.member[i], false); 
	        real_t H = S.penaltyConstraints7(S.member[i].x, nP);
                for (int IOP = 0; IOP < nP; IOP++){
                    for (int j = 0; j < nOO; j++){
                        if (S.member[i].open[j].N_iter(IOP) == optimizer.numIterAll){
                            chyba = S.member[i].open[j].CH_f(IOP) + S.member[i].open[j].CH(IOP);
                            S.member[i].open[j].f += chyba;
                            S.member[i].f += chyba;
                        }
                    }
                }
                S.member[i].penaltyF = hCoef*H;
                S.member[i].f += hCoef*H;
                S.member[i].actual_best(); // actualization fp (best value of the member)
            }
            //==================================== write to files =================================================================    
 
           ofile.open(fileHybrid.str(), std::ios_base::app);sfile.open(fileSummary.str(), std::ios_base::app);
            ofile << iter+1 << ". iterace: =========================================================================================================" << "\n";
            sfile << iter+1 << ". iterace: =========================================================================================================" << "\n";
            ofile.close(); sfile.close();
            for(int i = 0; i < nM; i++){  
                writeToFile(ofile, fileHybrid, S, i, optimizer.allOnameOfOpen);
            }
            if(iter == 0){
                writeToFileFg(ffile, fileFg, S, 0, optimizer.allOnameOfOpen);
                writeToFile(sfile, fileSummary, S, 0, optimizer.allOnameOfOpen);
            }
            //==================================== updationg previousSolutionBest ==================================================
 
           real_t fgOld = S.get_bestValue();
            S.bestMember();// actualization of g
            if (fgOld > S.get_bestValue()){
                bi = S.get_bestIndex();
                for (int j = 0; j < nOO; j++){ 
                    optimizer.allOpreviousSolutionBest(j) = S.member[bi].open[j].previousSolution;
                }
                for(int i = 0; i < nM; i++){
                    for (int IOP = 0; IOP < nP; IOP++){
                        for (int j = 0; j < nOO; j++){ 
                            if(S.member[i].open[j].N_iter(IOP) == optimizer.numIterAll && S.member[bi].open[j].N_iter(IOP) < optimizer.numIterAll){
                                S.member[i].open[j].previousSolution.col(IOP) = optimizer.allOpreviousSolutionBest(j).col(IOP);
                            }
                        }
                    }
                }
                writeToFile(sfile, fileSummary, S, bi, optimizer.allOnameOfOpen);
                writeToFileFg(ffile, fileFg, S, bi, optimizer.allOnameOfOpen);
            }
            else {
                ffile.open(fileFg.str(), std::ios_base::app); ffile << S.get_bestValue() << "\n"; ffile.close();
            }
            //======================================= update and write g ====================================================
 
           g = S.get_bestPosition(); // best design parameters of the swarm
            writeToFileG(ofile, fileHybrid, S);                
            writeToFileG(sfile, fileSummary, S);                
        }

        x0 = g; 
        //=========================================== end of cycle and possible of restart ============================================================
 
       if (cycle < (nCycles-1)){
            ofile.open(fileHybrid.str(), std::ios_base::app);sfile.open(fileSummary.str(), std::ios_base::app);
            ofile << cycle+1 << ". restart: =========================================================================================================" << "\n";
            sfile << cycle+1 << ". restart: =========================================================================================================" << "\n";
            ofile.close(); sfile.close();
            std::ostringstream strs_cycle;  strs_cycle << cycle + 1;
            saveDesignVariables(S.member[0].x, nP, "designVariables_middle" + strs_cycle.str() + ".txt"); 
            optimizer.flowSolver_save(x0,"middle" + strs_cycle.str(),false);//saving best member parameters and solution
        }

    } //============================================ end of optimization =======================================================================

    saveDesignVariables(x0, nP, "designVariables_final.txt");
    optimizer.flowSolver_save(x0,"final",false);

} // !draw_solution

    return 0;
}

