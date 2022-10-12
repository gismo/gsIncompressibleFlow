/*
    Author(s): J. Egermaier
*/

#ifdef _MSC_VER // to be removed
#define _USE_MATH_DEFINES
#endif

#ifndef UWBGEOMETRYSTEP2D_H
#define UWBGEOMETRYSTEP2D_H

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <gismo.h>

#include <math.h>
#include <string.h>
#include <sstream>

using namespace gismo;

template<class T> void geometryStep(std::map<std::string, gsVector<std::string>> parList, std::string outputDIR)
{
    int addRefPart, dim;
    gsVector<int> numRefine, numRefineU, numRefineWalls, numRefineCorner, deg;
    T a, b, c, a_in;

    std::string GSF, geomType;
    get_parameter(parList,geomType,"geometry_type");
    get_parameter(parList,GSF,"geometry_settings_file");
    std::map<std::string, gsVector<std::string>> paramList = readInputFile(GSF);

    get_parameter(paramList,numRefine,"numRefine");
    get_parameter(paramList,numRefineU,"numRefineU");
    get_parameter(paramList,numRefineWalls,"numRefineWalls");
    get_parameter(paramList,numRefineCorner,"numRefineCorner");
    get_parameter(paramList,deg,"deg");
    get_parameter(paramList,a,"a");
    get_parameter(paramList,a_in,"a_in");
    get_parameter(paramList,b,"b");
    get_parameter(paramList,c,"c");
    get_parameter(paramList,addRefPart,"addRefPart");

    if (geomType == "step3D"){
        dim = 3;
    } else {
        dim = 2;
    }

    for (int i = 0; i < numRefine.size(); i++){
        for (int j = 0; j < numRefineCorner.size(); j++){
            for (int k = 0; k < numRefineWalls.size(); k++){
                for (int l = 0; l < deg.size(); l++){
                    gsMultiPatch<T> patches;
                    switch(dim)
                    {
                    case 2: patches = BSplineStep2D<T>(deg(l), a, b, a_in);
                        c = 0.;
                        break;
                    case 3: patches = BSplineStep3D<T>(deg(l), a, b, c, a_in);
                        break;
                    default: GISMO_ERROR("Wrong dimension!");
                        break;
                    }
                    gsMultiBasis<T> tbasis(patches);

                    refineBasis_step<T>(tbasis, numRefine(i), numRefineWalls(k), numRefineCorner(j), numRefineU(0), addRefPart, dim, a, b, c);

                    std::string fileName = geomType + "_ref_" + util::to_string(numRefine(i)) + "_" +
                                            util::to_string(numRefineCorner(j)) + "_" + util::to_string(numRefineWalls(k)) + "_deg_" + util::to_string(deg(l));
                    gsFileData<> fdpatches;
                    fdpatches << patches;
                    fdpatches.save(outputDIR + fileName + "_patches.xml");

                    gsFileData<> fdtbasis;
                    fdtbasis << tbasis;
                    fdtbasis.save(outputDIR + fileName + "_tbasis.xml");
                }
            }
        }
    }
 }

#endif
