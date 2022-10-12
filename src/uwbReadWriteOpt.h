/*
    Author(s): J. Egermaier
*/

#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <gismo.h>

using namespace gismo;

/*
void createSH(std::string dirName, std::string fileName)
{
    std::ofstream ofile;
    ofile.open(fileName);
    ofile << "dirPath=$(dirname $0)" << "\n";
//    ofile << "dirPath='" + dirName + "'" << "\n";
    ofile << "cd ${dirPath}" << "\n";  
    ofile << "mkdir pokusicek" << "\n";
    ofile.close();
}
*/

std::map <std::string, gsVector<std::string>> readInputFile(std::string inputFile)
{
    std::map <std::string, gsVector<std::string>> param_list;

    std::ifstream cFile;
    cFile.open(inputFile.c_str());
    if (cFile.is_open()){
        std::string line;
        while(getline(cFile, line)){
            line.erase(remove_if(line.begin(), line.end(), isspace),line.end());
            if(line[0] == '#' || line.empty())
                continue;

            int delimiterPos = line.find("=");
            std::string name = line.substr(0, delimiterPos);
            line = line.substr(delimiterPos + 1);
                 
            std::vector<std::string> var;
            while ( line.find(",") != std::string::npos){
                int commaPos = line.find(",");
                var.push_back(line.substr(0, commaPos));
                line = line.substr(commaPos + 1);
            }    
            var.push_back(line);                 

            int varSize = var.size();
            gsVector<std::string> variable (varSize);
            for (int i = 0; i < varSize; i++){
                variable(i) = var[i];
            }

            param_list[name] = variable;
        }
        cFile.close();
    }
    else { gsWarn << "Unable to open param file: " << inputFile << '\n'; }
    
    return param_list;
}

//===================================================================================================================

void get_parameter(std::map <std::string, gsVector<std::string>> param_list, int & P, std::string parameter)
{
    P = stoi(param_list[parameter](0));
}

void get_parameter(std::map <std::string, gsVector<std::string>> param_list, gsVector<int> & P, std::string parameter)
{
    P.resize(param_list[parameter].size());
    for (int i = 0; i < param_list[parameter].size(); i++){
        P(i) = stoi(param_list[parameter](i));
    }
}

void get_parameter(std::map <std::string, gsVector<std::string>> param_list, real_t & R, std::string parameter)
{
    R = atof(param_list[parameter](0).c_str());    
}

void get_parameter(std::map <std::string, gsVector<std::string>> param_list, gsVector<real_t> & R, std::string parameter)
{
    R.resize(param_list[parameter].size());
    for (int i = 0; i < param_list[parameter].size(); i++){
        R(i) = atof(param_list[parameter](i).c_str()); 
    }
}

void get_parameter(std::map <std::string, gsVector<std::string>> param_list, std::string & S, std::string parameter)
{
    S = param_list[parameter](0);    
}

void get_parameter(std::map <std::string, gsVector<std::string>> param_list, gsVector<std::string> & S, std::string parameter)
{
    S.resize(param_list[parameter].size());
    for (int i = 0; i < param_list[parameter].size(); i++){
        S(i) = param_list[parameter](i); 
    }
}

void get_parameter(std::map <std::string, gsVector<std::string>> param_list, bool & B, std::string parameter)
{
    if ((param_list[parameter](0) == "true")|(param_list[parameter](0) == "1")){
        B = true;
    } else {
        B = false;
    }    
}

void get_parameter(std::map <std::string, gsVector<std::string>> param_list, gsVector<bool> & B, std::string parameter)
{
    B.resize(param_list[parameter].size());
    for(int i = 0; i < param_list[parameter].size(); i++){
        if ((param_list[parameter](i) == "true")|(param_list[parameter](i) == "1")){
            B(i) = true;
        } else {
            B(i) = false;
        }
    }
}

//=================================================================================================

void check_vector_length(gsVector<real_t> V, int L)
{
    GISMO_ENSURE( V.size() == L, "Parameter vector has wrong length! "); 
}

void check_vector_length(gsVector<int> V, int L)
{
    GISMO_ENSURE( V.size() == L, "Parameter vector has wrong length!"); 
}

void check_vector_length(gsVector<std::string> V, int L)
{
    GISMO_ENSURE( V.size() == L, "Parameter vector has wrong length!"); 
}

void check_vector_length(gsVector<bool> V, int L)
{
    GISMO_ENSURE( V.size() == L, "Parameter vector has wrong length!"); 
}


//=================================================================================================

std::map<std::string, real_t> readSolution (int nn, std::string solutionFile)
{
    std::map <std::string, real_t> param_list;
    gsVector<std::string> name (nn), value_string (nn);
    gsVector<real_t> value (nn);

    std::ifstream cFile;
    cFile.open(solutionFile.c_str());
    if (cFile.is_open()){
        int cur = 0;
        int start, end;
        std::string line;
        std::getline(cFile, line);
        for (int i = 0; i < nn; i++){
            while (isspace(line[cur])){cur++;}
            start = cur;
            while (!isspace(line[cur])){cur++;}
            end = cur; 
            name(i) = line.substr(start, end-start);
        }
        std::getline(cFile, line);
        std::getline(cFile, line);
        cur = 0;
        for (int i = 0; i < nn; i++){
            while (isspace(line[cur])){cur++;}
            start = cur; 
            while (!isspace(line[cur])){cur++;}
            end = cur; 
            value_string[i] = line.substr(start, end-start);
            value(i) = atof(value_string[i].c_str());
        }
        cFile.close();
    }
    else{ gsWarn << "Unable to open param file." << '\n'; }
    for (int i = 0; i < nn; i++){
        param_list[name[i]] = value[i];
    }
    return param_list;  
}            

//=============================================

void readCSV(gsVector<real_t> & prutok, gsVector<real_t> & ucinnost, std::string csvFile)
{
    std::ifstream file;
    std::string line;
    std::vector<real_t> U, Q;

    file.open(csvFile.c_str());
    std::getline(file, line);
    while(file.good()){
        std::getline(file, line, ',');
        std::getline(file, line, ',') ;
        std::getline(file, line, ',') ;
        std::getline(file, line, ',') ;
        Q.push_back(atof(line.c_str()));
        std::getline(file, line, ',') ; 
        std::getline(file, line) ; 
        U.push_back(atof(line.c_str()));
    }

    int delka = U.size();
    ucinnost.resize(delka);
    prutok.resize(delka);
    for (int i = 0; i < delka; i++){
        ucinnost(i) = U[i];
        prutok(i) = Q[i];
    }
    file.close();   
}

//============================================================

void copyAndRunSH(std::string beta, std::string dirName)
{        
    //int retSys;
    std::string inputSH = beta + "_opt.sh";
    std::string outputSH = dirName + "/" + beta + "_opt.sh";
    std::string copySH = "cp " + inputSH + " " + outputSH;
    //retSys = system(copySH.c_str());
 
    std::string chmodFile = "chmod +x " + outputSH; 
    //retSys = system(chmodFile.c_str());
    std::string runSH = "./" + outputSH;
    //retSys = system(runSH.c_str());
}  

//=============================================

void get_initialValues(gsVector<real_t> & x0, gsVector<real_t> & lowerBounds, gsVector<real_t> & upperBounds, gsVector<real_t> & max_tol, int & nP, int & nDV)
{
    gsVector<real_t> camber_x, camber_y, leading_angle, trailing_angle, thickness_x, thickness_y, output_angle, radius, angle, ending_offset, chord_length, rotation_center_x, rotation_center_y;
    gsVector<real_t> LBcamber_x, LBcamber_y, LBleading_angle, LBtrailing_angle, LBthickness_x, LBthickness_y, LBoutput_angle, LBradius,
                     LBangle, LBending_offset, LBchord_length, LBrotation_center_x, LBrotation_center_y;
    gsVector<real_t> UBcamber_x, UBcamber_y, UBleading_angle, UBtrailing_angle, UBthickness_x, UBthickness_y, UBoutput_angle, UBradius, 
                     UBangle, UBending_offset, UBchord_length, UBrotation_center_x, UBrotation_center_y;
    gsVector<real_t> RDcamber_x, RDcamber_y, RDleading_angle, RDtrailing_angle, RDthickness_x, RDthickness_y, RDoutput_angle, RDradius, 
                     RDangle, RDending_offset, RDchord_length, RDrotation_center_x, RDrotation_center_y;

    std::map <std::string, gsVector<std::string>> paramList;

    paramList = readInputFile("designVariables.txt");
    get_parameter(paramList,camber_x,"camber_x");
    nP = camber_x.size();
    get_parameter(paramList,camber_y,"camber_y");                   check_vector_length(camber_y, nP);
    get_parameter(paramList,leading_angle,"leading_angle");         check_vector_length(leading_angle, nP);
    get_parameter(paramList,trailing_angle,"trailing_angle");       check_vector_length(trailing_angle, nP);
    get_parameter(paramList,thickness_x,"thickness_x");             check_vector_length(thickness_x, nP);
    get_parameter(paramList,thickness_y,"thickness_y");             check_vector_length(thickness_y, nP);
    get_parameter(paramList,output_angle,"output_angle");           check_vector_length(output_angle, nP);
    get_parameter(paramList,radius,"radius");                       check_vector_length(radius, nP);
    get_parameter(paramList,angle,"angle");                         check_vector_length(angle, nP);
    get_parameter(paramList,ending_offset,"ending_offset");         check_vector_length(ending_offset, nP);
    get_parameter(paramList,chord_length,"chord_length");           check_vector_length(chord_length, nP);
    get_parameter(paramList,rotation_center_x,"rotation_center_x"); check_vector_length(rotation_center_x, nP);
    get_parameter(paramList,rotation_center_y,"rotation_center_y"); check_vector_length(rotation_center_y, nP);

    paramList = readInputFile("lowerBounds.txt");
    get_parameter(paramList,LBcamber_x,"camber_x");                   check_vector_length(LBcamber_x, nP);
    get_parameter(paramList,LBcamber_y,"camber_y");                   check_vector_length(LBcamber_y, nP);
    get_parameter(paramList,LBleading_angle,"leading_angle");         check_vector_length(LBleading_angle, nP);
    get_parameter(paramList,LBtrailing_angle,"trailing_angle");       check_vector_length(LBtrailing_angle, nP);
    get_parameter(paramList,LBthickness_x,"thickness_x");             check_vector_length(LBthickness_x, nP);
    get_parameter(paramList,LBthickness_y,"thickness_y");             check_vector_length(LBthickness_y, nP);
    get_parameter(paramList,LBoutput_angle,"output_angle");           check_vector_length(LBoutput_angle, nP);
    get_parameter(paramList,LBradius,"radius");                       check_vector_length(LBradius, nP);
    get_parameter(paramList,LBangle,"angle");                         check_vector_length(LBangle, nP);
    get_parameter(paramList,LBending_offset,"ending_offset");         check_vector_length(LBending_offset, nP);
    get_parameter(paramList,LBchord_length,"chord_length");           check_vector_length(LBchord_length, nP);
    get_parameter(paramList,LBrotation_center_x,"rotation_center_x"); check_vector_length(LBrotation_center_x, nP);
    get_parameter(paramList,LBrotation_center_y,"rotation_center_y"); check_vector_length(LBrotation_center_y, nP);

    paramList = readInputFile("upperBounds.txt");
    get_parameter(paramList,UBcamber_x,"camber_x");                   check_vector_length(UBcamber_x, nP);
    get_parameter(paramList,UBcamber_y,"camber_y");                   check_vector_length(UBcamber_y, nP);
    get_parameter(paramList,UBleading_angle,"leading_angle");         check_vector_length(UBleading_angle, nP);
    get_parameter(paramList,UBtrailing_angle,"trailing_angle");       check_vector_length(UBtrailing_angle, nP);
    get_parameter(paramList,UBthickness_x,"thickness_x");             check_vector_length(UBthickness_x, nP);
    get_parameter(paramList,UBthickness_y,"thickness_y");             check_vector_length(UBthickness_y, nP);
    get_parameter(paramList,UBoutput_angle,"output_angle");           check_vector_length(UBoutput_angle, nP);
    get_parameter(paramList,UBradius,"radius");                       check_vector_length(UBradius, nP);
    get_parameter(paramList,UBangle,"angle");                         check_vector_length(UBangle, nP);
    get_parameter(paramList,UBending_offset,"ending_offset");         check_vector_length(UBending_offset, nP);
    get_parameter(paramList,UBchord_length,"chord_length");           check_vector_length(UBchord_length, nP);
    get_parameter(paramList,UBrotation_center_x,"rotation_center_x"); check_vector_length(UBrotation_center_x, nP);
    get_parameter(paramList,UBrotation_center_y,"rotation_center_y"); check_vector_length(UBrotation_center_y, nP);

    paramList = readInputFile("relativeDisplacements.txt");
    get_parameter(paramList,RDcamber_x,"camber_x");                   check_vector_length(RDcamber_x, nP);
    get_parameter(paramList,RDcamber_y,"camber_y");                   check_vector_length(RDcamber_y, nP);
    get_parameter(paramList,RDleading_angle,"leading_angle");         check_vector_length(RDleading_angle, nP);
    get_parameter(paramList,RDtrailing_angle,"trailing_angle");       check_vector_length(RDtrailing_angle, nP);
    get_parameter(paramList,RDthickness_x,"thickness_x");             check_vector_length(RDthickness_x, nP);
    get_parameter(paramList,RDthickness_y,"thickness_y");             check_vector_length(RDthickness_y, nP);
    get_parameter(paramList,RDoutput_angle,"output_angle");           check_vector_length(RDoutput_angle, nP);
    get_parameter(paramList,RDradius,"radius");                       check_vector_length(RDradius, nP);
    get_parameter(paramList,RDangle,"angle");                         check_vector_length(RDangle, nP);
    get_parameter(paramList,RDending_offset,"ending_offset");         check_vector_length(RDending_offset, nP);
    get_parameter(paramList,RDchord_length,"chord_length");           check_vector_length(RDchord_length, nP);
    get_parameter(paramList,RDrotation_center_x,"rotation_center_x"); check_vector_length(RDrotation_center_x, nP);
    get_parameter(paramList,RDrotation_center_y,"rotation_center_y"); check_vector_length(RDrotation_center_y, nP);

    nDV = paramList.size();
    x0.resize(nP*nDV); lowerBounds.resize(nP*nDV); upperBounds.resize(nP*nDV); max_tol.resize(nP*nDV);

    int j = 0;
    x0.middleRows(j,nP) = camber_x; j += nP;
    x0.middleRows(j,nP) = camber_y; j += nP;
    x0.middleRows(j,nP) = leading_angle; j += nP;
    x0.middleRows(j,nP) = trailing_angle; j += nP;
    x0.middleRows(j,nP) = thickness_x; j += nP;
    x0.middleRows(j,nP) = thickness_y; j += nP;
    x0.middleRows(j,nP) = output_angle; j += nP;
    x0.middleRows(j,nP) = radius; j += nP;
    x0.middleRows(j,nP) = angle; j += nP;
    x0.middleRows(j,nP) = ending_offset; j += nP;
    x0.middleRows(j,nP) = chord_length; j += nP;
    x0.middleRows(j,nP) = rotation_center_x; j += nP;
    x0.middleRows(j,nP) = rotation_center_y; j += nP;

    j = 0;
    lowerBounds.middleRows(j,nP) = LBcamber_x; j += nP;
    lowerBounds.middleRows(j,nP) = LBcamber_y; j += nP;
    lowerBounds.middleRows(j,nP) = LBleading_angle; j += nP;
    lowerBounds.middleRows(j,nP) = LBtrailing_angle; j += nP;
    lowerBounds.middleRows(j,nP) = LBthickness_x; j += nP;
    lowerBounds.middleRows(j,nP) = LBthickness_y; j += nP;
    lowerBounds.middleRows(j,nP) = LBoutput_angle; j += nP;
    lowerBounds.middleRows(j,nP) = LBradius; j += nP;
    lowerBounds.middleRows(j,nP) = LBangle; j += nP;
    lowerBounds.middleRows(j,nP) = LBending_offset; j += nP;
    lowerBounds.middleRows(j,nP) = LBchord_length; j += nP;
    lowerBounds.middleRows(j,nP) = LBrotation_center_x; j += nP;
    lowerBounds.middleRows(j,nP) = LBrotation_center_y; j += nP;

    j = 0;
    upperBounds.middleRows(j,nP) = UBcamber_x; j += nP;
    upperBounds.middleRows(j,nP) = UBcamber_y; j += nP;
    upperBounds.middleRows(j,nP) = UBleading_angle; j += nP;
    upperBounds.middleRows(j,nP) = UBtrailing_angle; j += nP;
    upperBounds.middleRows(j,nP) = UBthickness_x; j += nP;
    upperBounds.middleRows(j,nP) = UBthickness_y; j += nP;
    upperBounds.middleRows(j,nP) = UBoutput_angle; j += nP;
    upperBounds.middleRows(j,nP) = UBradius; j += nP;
    upperBounds.middleRows(j,nP) = UBangle; j += nP;
    upperBounds.middleRows(j,nP) = UBending_offset; j += nP;
    upperBounds.middleRows(j,nP) = UBchord_length; j += nP;
    upperBounds.middleRows(j,nP) = UBrotation_center_x; j += nP;
    upperBounds.middleRows(j,nP) = UBrotation_center_y; j += nP;

    j = 0;
    max_tol.middleRows(j,nP) = RDcamber_x; j += nP;
    max_tol.middleRows(j,nP) = RDcamber_y; j += nP;
    max_tol.middleRows(j,nP) = RDleading_angle; j += nP;
    max_tol.middleRows(j,nP) = RDtrailing_angle; j += nP;
    max_tol.middleRows(j,nP) = RDthickness_x; j += nP;
    max_tol.middleRows(j,nP) = RDthickness_y; j += nP;
    max_tol.middleRows(j,nP) = RDoutput_angle; j += nP;
    max_tol.middleRows(j,nP) = RDradius; j += nP;
    max_tol.middleRows(j,nP) = RDangle; j += nP;
    max_tol.middleRows(j,nP) = RDending_offset; j += nP;
    max_tol.middleRows(j,nP) = RDchord_length; j += nP;
    max_tol.middleRows(j,nP) = RDrotation_center_x; j += nP;
    max_tol.middleRows(j,nP) = RDrotation_center_y; j += nP;
}

//=================================================

void saveDesignVariables(gsVector<real_t> x, int nP, std::string file)
{
    std::ofstream ofile; ofile.open(file);
    int j = 0;
    ofile << "camber_x = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "camber_y = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "leading_angle = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "trailing_angle = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "thickness_x = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "thickness_y = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "output_angle = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "radius =  "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "angle = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "ending_offset = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "chord_length = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "rotation_center_x = "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n"; j += nP;
    ofile << "rotation_center_y =  "; for (int i = 0; i < nP; i++){ if (i>0){ ofile << ", ";} ofile << x.middleRows(j,nP)(i);} ofile << "\n";
    ofile.close();
}

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

void write_sparse_xml_to_dat(std::string outputfile, std::string inputfile)
{ 
    std::ofstream ofile; ofile.open(outputfile);

    std::ifstream infile(inputfile);
    std::string content = "";
    std::string newChar, oldline, line;
    int count, endPos, varInt;
    //bool write = false;
    gsVector<std::string> var(2);

    while (getline(infile, line)){
        if (line.size() > 3){
        if ((line[0]==' ')||(line[0]=='<')){
            endPos = line.find(">");
            line.erase(0, endPos + 1);
            newChar = line[0];
            if ((newChar == "0")||(newChar == "1")||(newChar == "2")||(newChar == "3")||(newChar == "4")||(newChar == "5")||(newChar == "6")||(newChar == "7")||(newChar == "8")||(newChar == "9")){
                for (int i = 0; i < 2; i++){
                    int commaPos = line.find(" ");
                    var(i) = line.substr(0, commaPos);
                    line = line.substr(commaPos + 1);
                    varInt = stoi(var(i));
                    varInt ++;
                    std::ostringstream oss; oss << varInt; var(i) = oss.str();
                }                   
                oldline = line;
                line = var(0) + " " + var(1) + " " + oldline;  
                content += line + "\n";
            }
        } else {
            for (int i = 0; i < 2; i++){
                int commaPos = line.find(" ");
                var(i) = line.substr(0, commaPos);
                line = line.substr(commaPos + 1);
                varInt = stoi(var(i));
                varInt ++;
                std::ostringstream oss; oss << varInt; var(i) = oss.str();
            }                   
            oldline = line;
            line = var(0) + " " + var(1) + " " + oldline;
            content += line + "\n";
        }
    }
    }
    count--;
    content.erase(content.end()-1);     // erase last character
    infile.close();
    ofile << content;                 // output
    ofile.close();
}

void write_sparse_direct_to_dat(std::string outputfile, gsSparseMatrix<real_t> M)
{
    gsFileData<real_t> fd;
    fd << M;
    std::string inputfile = "tempM.xml";
    fd.save(inputfile);

    std::ofstream ofile; ofile.open(outputfile);

    std::ifstream infile(inputfile);
    std::string content = "";
    std::string newChar, oldline, line;
    int count, endPos, varInt;
    //bool write = false;
    gsVector<std::string> var(2);

    while (getline(infile, line)){
        if (line.size() > 3){
        if ((line[0]==' ')||(line[0]=='<')){
            endPos = line.find(">");
            line.erase(0, endPos + 1);
            newChar = line[0];
            if ((newChar == "0")||(newChar == "1")||(newChar == "2")||(newChar == "3")||(newChar == "4")||(newChar == "5")||(newChar == "6")||(newChar == "7")||(newChar == "8")||(newChar == "9")){
                for (int i = 0; i < 2; i++){
                    int commaPos = line.find(" ");
                    var(i) = line.substr(0, commaPos);
                    line = line.substr(commaPos + 1);
                    varInt = stoi(var(i));
                    varInt ++;
                    std::ostringstream oss; oss << varInt; var(i) = oss.str();
                }
                oldline = line;
                line = var(0) + " " + var(1) + " " + oldline;
                content += line + "\n";
            }
        } else {
            for (int i = 0; i < 2; i++){
                int commaPos = line.find(" ");
                var(i) = line.substr(0, commaPos);
                line = line.substr(commaPos + 1);
                varInt = stoi(var(i));
                varInt ++;
                std::ostringstream oss; oss << varInt; var(i) = oss.str();
            }
            oldline = line;
            line = var(0) + " " + var(1) + " " + oldline;
            content += line + "\n";
        }
    }
    }
    count--;
    content.erase(content.end()-1);     // erase last character
    infile.close();
    ofile << content;                 // output
    ofile.close();

    std::string remInput = "rm " + inputfile;
    int retsys = system(remInput.c_str());
}

gsVector<std::string> combine_names(gsVector<std::string> start, gsVector<std::string> ref, gsVector<std::string> body, gsVector<std::string> end)
{
    int count = 0;
    int outSize = start.size()*ref.size()*body.size()*end.size();
    gsVector<std::string> output (outSize);
    for (int i = 0; i < start.size(); i++){
        for (int j = 0; j < ref.size(); j++){
            for (int k = 0; k < body.size(); k++){
                for (int l = 0; l < end.size(); l++){
                    output(count) = start(i) + ref(j) + body(k) + end(l);
                    count++;
                }
            }
        }
    }
    return output;
}

void saveToXMLandDat(std::vector<gsSparseMatrix<real_t>> A, std::string directory, std::string name)
{
    for (unsigned i = 0; i < A.size(); i++){
        gsFileData<real_t> fd;
        fd << A[i];
        std::string filename = directory + name + "(" + util::to_string(i) + ")";
        fd.save(filename);
        write_sparse_xml_to_dat(filename + ".dat", filename + ".xml");
    }
}

void saveToXMLandDat(std::vector<gsSparseMatrix<real_t, RowMajor>> A, std::string directory, std::string name)
{
    for (unsigned i = 0; i < A.size(); i++){
        gsFileData<real_t> fd;
        gsSparseMatrix<real_t> sM = A[i];
        fd << sM;
        std::string filename = directory + name + "(" + util::to_string(i) + ")";
        fd.save(filename);
        write_sparse_xml_to_dat(filename + ".dat", filename + ".xml");
    }
}

void saveMatrixCSRtoTXT(gsSparseMatrix<real_t> A, std::string filename){
    A.makeCompressed();
    int nnz = A.nonZeros();
    int ncols = A.cols();
    gsMatrix<real_t> V (nnz,1);
    gsMatrix<int> II (nnz,1);
    gsMatrix<int> OI (ncols+1,1);
    real_t * v = A.valuePtr();
    int * oi = A.outerIndexPtr();
    int * ii = A.innerIndexPtr();

    std::string fileV = "";
    std::string fileII = "";
    std::string fileOI = "";
    for (int i = 0; i < nnz; i++){
        fileV += util::to_string(v[i]) + "\n";
        fileII += util::to_string(ii[i]) + "\n";
    }
    for (int i = 0; i < ncols+1; i++){
        fileOI += util::to_string(oi[i]) + "\n";
    }
    gsInfo << v[0] << ", " << util::to_string(v[0]) << "\n";

    std::ofstream ofile;
    ofile.open(filename + "_v.txt");
    ofile << fileV;
    ofile.close();
    ofile.open(filename + "_ii.txt");
    ofile << fileII;
    ofile.close();
    ofile.open(filename + "_oi.txt");
    ofile << fileOI;
    ofile.close();
}


