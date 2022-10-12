#ifndef UWBDATAFILEREADER_H
#define UWBDATAFILEREADER_H

void readInitialSetting(std::string fileName,  gsVector<int>& inInt,  gsVector<real_t>& inRealT, std::string & inString,  gsVector<bool> & inBool)
{
    std::ifstream inFile;
    inFile.open(fileName);

    if (!inFile) {
        gsInfo << "Unable to open file \n";
        exit(1); // terminate with error
    }

    std::string str;
    for(int i = 0; i < 2*inInt.rows(); i++)
    {
        str.clear();
        std::getline(inFile,str);
       // gsInfo << str << "\n";
        if (str[0] == '#')
            continue;
        else
        {
            std::istringstream number(str);
            number >> inInt (i/2);
        }
    }

    for(int i = 0; i < 2*inRealT.rows(); i++)
    {
        str.clear();
        std::getline(inFile,str);
        if (str[0] == '#')
            continue;
        else
        {
            std::istringstream number(str);
            number >> inRealT(i/2);
        }
    }

    str.clear();
    std::getline(inFile,str);
    std::getline(inFile,str);
    inString = str;

    for(int i = 0; i < 2*inBool.rows(); i++)
    {
        str.clear();
        std::getline(inFile,str);
        if (str[0] == '#')
            continue;
        else
        {
            std::istringstream number(str);
            number >> inBool(i/2);
        }
    }
    inFile.close();
    return;
}

std::vector<double> getKnotVector(std::string const& Line)
{
  std::istringstream iss(Line);

  return std::vector<double>{
    std::istream_iterator<double>(iss),
    std::istream_iterator<double>()
  };
}

void readInitialGeometry(std::string fileName, gsVector<int> &inGeoInt, gsMatrix<int> &refineUniformSettings, gsMatrix<int> &refineLocalSettings, gsVector<real_t> &geomParams, gsVector<bool> &inGeoBool,  std::vector<real_t> &kvfit_knots)
{


    std::ifstream inFileGeom;
    inFileGeom.open(fileName);

    if (!inFileGeom) {
        gsInfo << "Unable to open file";
        exit(1); // terminate with error
    }


     std::string str;

    for(int i = 0; i < 2*2; i++)
    {
        str.clear();
        std::getline(inFileGeom,str);

       if (str[0] == '#')
           continue;
       else
        {  std::istringstream number(str);
            number >> inGeoInt(i/2);}
    }


    str.clear();
    std::getline(inFileGeom,str);
    str.clear();
    std::getline(inFileGeom,str);
    kvfit_knots = getKnotVector(str);




    for(int i = 0; i < 2*2; i++)
    {
        str.clear();
        std::getline(inFileGeom,str);


       if (str[0] == '#')
           continue;
       else
        {  std::istringstream number(str);
            number >> inGeoBool(i/2);}
    }



    for(int i = 2*2; i < 2*5; i++)
    {
        str.clear();
        std::getline(inFileGeom,str);

       if (str[0] == '#')
           continue;
       else
        {  std::istringstream number(str);
            number >> inGeoInt(i/2);}
    }

    gsInfo << inGeoInt;

    refineUniformSettings.resize(inGeoInt(4),5);
    if (inGeoInt(4) == 0){
        gsInfo << "There is not =uniform= refinement of patch. \n";
    }
    else{
        for (int i = 0; i < refineUniformSettings.rows(); i++)
          for (int j = 0; j < 2*(refineUniformSettings.cols()); j++)
             {
              str.clear();
              std::getline(inFileGeom,str);
             if (str[0] == '#')
                 continue;
             else
              {  std::istringstream number(str);
                  number >> refineUniformSettings(i,j/2);}
              }
    }

  //  gsInfo << refineUniformSettings << "\n";
    str.clear();
    std::getline(inFileGeom,str);
    str.clear();
    std::getline(inFileGeom,str);
    std::istringstream number4(str);
    number4 >> inGeoInt(5);



    refineLocalSettings.resize(inGeoInt(5),5);

      if (inGeoInt(5) == 0){
          gsInfo << "There is not =local recursive= refinement of patch.";
      }
      else{
          for (int i = 0; i < refineLocalSettings.rows(); i++)
            for (int j = 0; j < 2*(refineLocalSettings.cols()); j++)
            {
             str.clear();
             std::getline(inFileGeom,str);
             //gsInfo << str << "\n";

            if (str[0] == '#')
                continue;
            else
             {  std::istringstream number(str);
                 number >> refineLocalSettings(i,j/2);}
         }
      }



      str.clear();
      std::getline(inFileGeom,str);
      str.clear();
      std::getline(inFileGeom,str);
      std::istringstream number6(str);
      number6 >> inGeoInt(6);

      geomParams.resize(inGeoInt(6));
      if (inGeoInt(6) == 0){
          gsInfo << "There are not geomParams.";
      }
      else{
      for (int j = 0; j < 2*geomParams.rows(); j++)
      {
       str.clear();
       std::getline(inFileGeom,str);
      if (str[0] == '#')
          continue;
      else
       {  std::istringstream number(str);
           number >> geomParams(j/2);}
       }
      }



 inFileGeom.close();
    return;
}


#endif // UWBDATAFILEREADER_H
