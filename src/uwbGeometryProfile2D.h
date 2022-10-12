/*
    Author(s): J. Egermaier
*/

#ifdef _MSC_VER // to be removed
#define _USE_MATH_DEFINES
#endif

#ifndef UWBGEOMETRYPROFILE2D_H
#define UWBGEOMETRYPROFILE2D_H

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <gismo.h>

//problem settings
#include "uwbRANSExamplesSetting.h"
#include "uwbGeometryCreators.h"
#include "uwbDataFileReader.h"

#include <math.h>
#include <string.h>
#include <sstream>

using namespace gismo;

void geometryProfile2D(std::map<std::string, gsVector<std::string>> paramList, std::string outputDIR)
{
    std::string geomSettingsFile;

    get_parameter(paramList,geomSettingsFile,"geometry_settings_file");

    gsVector<int> inGeoInt(7);
    gsMatrix<int> refineUniformSettings(0,0);
    gsMatrix<int> refineLocalSettings(0,0);
    gsVector<real_t> geomParams(0);
    gsVector<bool> inGeoBool(2);
    std::vector<real_t> kvfit_knots;

    readInitialGeometry(geomSettingsFile, inGeoInt, refineUniformSettings, refineLocalSettings, geomParams, inGeoBool, kvfit_knots);

    int geomChoice = inGeoInt(0);
    int index_of_profile = inGeoInt(1);
    int uniformRefine = inGeoInt(3);

    bool coarse = inGeoBool(0);
    bool makeLinearBool = inGeoBool(1);

    int makeLinearNumOfRefine = inGeoInt(2);

    gsVector<int> kvfit_print(kvfit_knots.size());
    for (int i = 0; i < kvfit_knots.size(); i++)
    {
        kvfit_print(i) = kvfit_knots[i];
    }

    int setting_of_blade_rotation = 0; //0-optimal, 1-maximal, 2-minimal
    int setting_of_guide_vanes = 0; //0-compute with the formula, 1- 56° GV for optimal and 60° GV for maximal, 2-68° GV for optimal a 72° GV for maximal

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.intStrategy = iFace::glue;
    opt.dirValues = dirichlet::interpolation;
    //opt.dirValues = dirichlet::l2Projection; //default

    std::string profile = "_blade" + std::to_string(index_of_profile);//strs_profile.str();

    const unsigned num_blades = 4;
    const unsigned num_bladeprofiles = 7;

    real_t viscosity;
    get_parameter(paramList,viscosity,"visc");
    int plot_pts = 1700;

    uwbRANSProfileExample<real_t> problemSettings(viscosity, profile, index_of_profile, num_bladeprofiles, num_blades, geomChoice, kvfit_knots, coarse, geomParams, plot_pts);
    problemSettings.setParameters(setting_of_blade_rotation, setting_of_guide_vanes);

    gsMultiPatch<real_t> patches_start = problemSettings.makeMultiPatch();
//    gsFileData<> fdpatches_start;
//    fdpatches_start << patches_start;
//    fdpatches_start.save("filePatches_created.xml");
//    gsWriteParaview( patches_start, "patches_start", 50000, true);

    // ========================================= Define basis and refine =========================================
    gsMultiBasis<> tbasis_old(patches_start); // basis for RANS equations

    for(int i = 0; i < refineUniformSettings.rows(); i++)
        problemSettings.refineBasisUniformZones(tbasis_old, refineUniformSettings(i,0), refineUniformSettings(i,1), refineUniformSettings(i,2), refineUniformSettings(i,3), refineUniformSettings(i,4));

    for(int i = 0; i < uniformRefine; i++)
        tbasis_old.uniformRefine();

    for(int i = 0; i < refineLocalSettings.rows(); i++)
        problemSettings.refineBasisLocalZones(tbasis_old, refineLocalSettings(i,0), refineLocalSettings(i,1), refineLocalSettings(i,2), refineLocalSettings(i,3), refineLocalSettings(i,4));

    //-------------------------------------------------------------------------------------------------

    gsMultiPatch<> patches;

    int US = 0;
    int LS = 0;
    std::string isLin, nRefUn;
    if (makeLinearBool){
        patches = makeLinearBasis(patches_start, tbasis_old, makeLinearNumOfRefine);
        isLin = "_deg_1";
        nRefUn = util::to_string(uniformRefine-1);//(makeLinearNumOfRefine-2);
    } else {
        patches = patches_start;
        isLin = "_deg_3";
        nRefUn = util::to_string(uniformRefine-1);
    }

    gsMultiBasis<> tbasis(patches);
    if (!makeLinearBool){
        tbasis = tbasis_old;
    }

    if (refineUniformSettings.rows() > 0){US = refineUniformSettings(0,4);}
    if (refineLocalSettings.rows() > 0){LS = refineLocalSettings(0,4);}

    std::string fileName = "profile2D_profile" + util::to_string(index_of_profile) + "_ref_" + nRefUn + "_" +
                            util::to_string(US-2) + "_" + util::to_string(LS) + isLin;
    //gsWriteParaview(patches, "patches_beforeBasisRefine", 50000, true);
    gsFileData<> fdpatches;
    fdpatches << patches;
    fdpatches.save(outputDIR + fileName + "_patches.xml");

    gsFileData<> fdtbasis;
    fdtbasis << tbasis;
    fdtbasis.save(outputDIR + fileName + "_tbasis.xml");
    //==================================================================================================

 }

#endif
