#ifndef UWBKAPLANTURBINEGUIDEVANE
#define UWBKAPLANTURBINEGUIDEVANE

using namespace gismo;

template <class T>
class KaplanTurbineGuideVane
{
public:
    // Constructors and destructors
    KaplanTurbineGuideVane(T thisCamberX, T thisCamberY, T thisLeadingAngle, T thisTrailingAngle,
                             T thisThicknessX, T thisThicknessY, T thisEndingOffset,
                             T thisOutputAngle, T thisChordLength, T thisAngle,
                             T thisRotationCenterX, T thisRotationCenterY, gsMatrix<T> thisCircularMeshSetting, gsMatrix<T> thisConeSetting);
    KaplanTurbineGuideVane(T thisCamberX, T thisCamberY, T thisLeadingAngle, T thisTrailingAngle,
                             T thisThicknessX, T thisThicknessY, T thisEndingOffset,
                             T thisOutputAngle, T thisRadius, T thisChordLength, T thisAngle,
                             T thisRotationCenterX, T thisRotationCenterY, gsMatrix<T> thisCircularMeshSetting, gsMatrix<T> thisConeSetting);
    ~KaplanTurbineGuideVane() {

            delete m_BladeProfiles;


    }

    // Setters
    void setSuctionSideSurface(gsTensorBSpline<2,T> bspl) { mySuctionSideSurface = bspl; }
    gsTensorBSpline<2,T> getSuctionSideSurface() { return mySuctionSideSurface; }
    void setPressureSideSurface(gsTensorBSpline<2,T> bspl) { myPressureSideSurface = bspl; }
    gsTensorBSpline<2,T> getPressureSideSurface() { return myPressureSideSurface; }
    void setSuctionSideSurfaceAfterTrimm(gsTensorBSpline<2,T> bspl) { mySuctionSideSurfaceAfterTrimm = bspl; }
    gsTensorBSpline<2,T> getSuctionSideSurfaceAfterTrimm() { return mySuctionSideSurfaceAfterTrimm; }
    void setPressureSideSurfaceAfterTrimm(gsTensorBSpline<2,T> bspl) { myPressureSideSurfaceAfterTrimm = bspl; }
    gsTensorBSpline<2,T> getPressureSideSurfaceAfterTrimm() { return myPressureSideSurfaceAfterTrimm; }

    // Other functions
    int compute(bool plot);
    int trimming(bool plot);

private:

    T myCamberX;
    T myCamberY;
    T myLeadingAngle;
    T myTrailingAngle;
    T myThicknessX;
    T myThicknessY;
    T myEndingOffset;
    T myOutputAngle;
    T myRadius;
    T myChordLength;
    T myAngle;
    T myRotationCenterX;
    T myRotationCenterY;
    gsMatrix<T> myCircularMeshSetting;
    gsMatrix<T> myConeSetting;

    gsTensorBSpline<2,T> myPressureSideSurface;
    gsTensorBSpline<2,T> mySuctionSideSurface;
    gsTensorBSpline<2,T> myPressureSideSurfaceAfterTrimm;
    gsTensorBSpline<2,T> mySuctionSideSurfaceAfterTrimm;


    BladeProfile<T> * m_BladeProfiles;

};

    template <class T>
  KaplanTurbineGuideVane<T>::KaplanTurbineGuideVane(T thisCamberX, T thisCamberY, T thisLeadingAngle, T thisTrailingAngle,
                          T thisThicknessX, T thisThicknessY, T thisEndingOffset,
                          T thisOutputAngle, T thisChordLength, T thisAngle,
                          T thisRotationCenterX, T thisRotationCenterY, gsMatrix<T> thisCircularMeshSetting, gsMatrix<T> thisConeSetting){
     myCamberX = thisCamberX;
     myCamberY = thisCamberY;
     myLeadingAngle = thisLeadingAngle;
     myTrailingAngle = thisTrailingAngle;
     myThicknessX = thisThicknessX;
     myThicknessY = thisThicknessY;
     myEndingOffset = thisEndingOffset;
     myOutputAngle = thisOutputAngle;
     myChordLength = thisChordLength;
     myAngle = thisAngle;
     myRotationCenterX = thisRotationCenterX;
     myRotationCenterY = thisRotationCenterY;
     myCircularMeshSetting = thisCircularMeshSetting;
     myConeSetting = thisConeSetting;
 }

  template <class T>
   KaplanTurbineGuideVane<T>::KaplanTurbineGuideVane(T thisCamberX, T thisCamberY, T thisLeadingAngle, T thisTrailingAngle,
                                                     T thisThicknessX, T thisThicknessY, T thisEndingOffset,
                                                     T thisOutputAngle, T thisRadius, T thisChordLength, T thisAngle,
                                                     T thisRotationCenterX, T thisRotationCenterY, gsMatrix<T> thisCircularMeshSetting, gsMatrix<T> thisConeSetting){
      myCamberX = thisCamberX;
      myCamberY = thisCamberY;
      myLeadingAngle = thisLeadingAngle;
      myTrailingAngle = thisTrailingAngle;
      myThicknessX = thisThicknessX;
      myThicknessY = thisThicknessY;
      myEndingOffset = thisEndingOffset;
      myOutputAngle = thisOutputAngle;
      myRadius = thisRadius;
      myChordLength = thisChordLength;
      myAngle = thisAngle;
      myRotationCenterX = thisRotationCenterX;
      myRotationCenterY = thisRotationCenterY;
      myCircularMeshSetting = thisCircularMeshSetting;
      myConeSetting = thisConeSetting;
  }

   template <class T>
   int KaplanTurbineGuideVane<T>::compute(bool plot) {

       int num_samples = 30;
       int num_samples3d = 50;

       gsBSpline<T> suction_side_curve;
       gsBSpline<T> pressure_side_curve;
       gsBSpline<T> suction_side_curve3d;
       gsBSpline<T> pressure_side_curve3d;
       gsKnotVector<T> kvfit(0, 1, 4, 4);
         gsKnotVector<T> kv(0, 1, 0, 2);
         gsTensorBSplineBasis<2, T> basis(kvfit,kv);
       BladeProfile<T> * pBladeProfile;

       gsInfo << "Construction of planar blade profile.";
       if (myRadius == 0) {
           pBladeProfile = new BladeProfile<T>(myCamberX, myCamberY, myLeadingAngle, myTrailingAngle, myThicknessX,
                                               myThicknessY, myEndingOffset, myOutputAngle, myChordLength,
                                               myAngle, myRotationCenterX, myRotationCenterY, myCircularMeshSetting, myConeSetting);
       }
       else {
           pBladeProfile = new BladeProfile<T>(myCamberX, myCamberY, myLeadingAngle, myTrailingAngle, myThicknessX,
                                               myThicknessY, myEndingOffset, myOutputAngle, myRadius, myChordLength,
                                               myAngle, myRotationCenterX, myRotationCenterY, myCircularMeshSetting, myConeSetting);
       }


       pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 2);
         pBladeProfile->compute3DGuideVane(suction_side_curve3d, pressure_side_curve3d, kvfit, num_samples3d);


        gsMatrix<T> parameter_points(1, num_samples3d+1);
        gsMatrix<T> suction_side_curve_points3d(num_samples3d+1, 3);
        suction_side_curve_points3d.setZero(num_samples3d+1, 3);
        gsMatrix<T> pressure_side_curve_points3d(num_samples3d+1, 3);
        pressure_side_curve_points3d.setZero(num_samples3d+1, 3);
        gsMatrix<T> suction_side_curve_points3d_top(num_samples3d+1, 3);
     suction_side_curve_points3d_top.setZero(num_samples3d+1, 3);
        gsMatrix<T> pressure_side_curve_points3d_top(num_samples3d+1, 3);
        pressure_side_curve_points3d_top.setZero(num_samples3d+1, 3);
        gsMatrix<T> suction_side_curve_points3d_bottom(num_samples3d+1, 3);
         suction_side_curve_points3d_bottom.setZero(num_samples3d+1, 3);
        gsMatrix<T> pressure_side_curve_points3d_bottom(num_samples3d+1, 3);
        pressure_side_curve_points3d_bottom.setZero(num_samples3d+1, 3);

        for (index_t i = 0; i < num_samples3d; i++) {
            //parameter_points(0,i) = (1-cos(pi/2*i/num_samples3d));
            parameter_points(0,i) = pow((real_t) i/num_samples3d,2);
        }

        parameter_points(0,num_samples3d) = 1;
        gsInfo << parameter_points << "\n";
        suction_side_curve_points3d = suction_side_curve3d.eval(parameter_points);
        suction_side_curve_points3d = suction_side_curve_points3d.transpose();
        pressure_side_curve_points3d =  pressure_side_curve3d.eval(parameter_points);
        pressure_side_curve_points3d = pressure_side_curve_points3d.transpose();


        // top profile on the cone
        real_t Sz = myCircularMeshSetting(5,0);
        real_t z0 = myConeSetting(1,2);


        const double PI = 3.14159265;
        real_t beta = ((PI/2) - myCircularMeshSetting(4,0));
        real_t c = math::tan(beta);
        real_t Ax,Ay,Az,frac;

        for (index_t i = 0; i < num_samples3d+1; i++) {


            Ax = suction_side_curve_points3d(i,0);
            Ay = suction_side_curve_points3d(i,1);
            Az = suction_side_curve_points3d(i,2);
             frac =(c*(Sz-z0)*(Sz-z0))/(-c*(Az-Sz)*(Sz-z0)+math::sqrt((Ax*Ax+Ay*Ay)*(Sz-z0)*(Sz-z0)));
            suction_side_curve_points3d_top(i,0) = Ax*frac;
            suction_side_curve_points3d_top(i,1) =  Ay*frac;
           suction_side_curve_points3d_top(i,2) =  Sz+(Az-Sz)*frac;

            Ax = pressure_side_curve_points3d(i,0);
            Ay = pressure_side_curve_points3d(i,1);
            Az = pressure_side_curve_points3d(i,2);
            frac =(c*(Sz-z0)*(Sz-z0))/(-c*(Az-Sz)*(Sz-z0)+math::sqrt((Ax*Ax+Ay*Ay)*(Sz-z0)*(Sz-z0)));
            pressure_side_curve_points3d_top(i,0) = Ax*frac;
            pressure_side_curve_points3d_top(i,1) =  Ay*frac;
            pressure_side_curve_points3d_top(i,2) =  Sz+(Az-Sz)*frac;

        }




        z0 = myConeSetting(2,2);



        for (index_t i = 0; i < num_samples3d+1; i++) {


            Ax = suction_side_curve_points3d(i,0);
            Ay = suction_side_curve_points3d(i,1);
            Az = suction_side_curve_points3d(i,2);
             frac =(c*(Sz-z0)*(Sz-z0))/(-c*(Az-Sz)*(Sz-z0)+math::sqrt((Ax*Ax+Ay*Ay)*(Sz-z0)*(Sz-z0)));
            suction_side_curve_points3d_bottom(i,0) = Ax*frac;
            suction_side_curve_points3d_bottom(i,1) =  Ay*frac;
           suction_side_curve_points3d_bottom(i,2) =  Sz+(Az-Sz)*frac;

            Ax = pressure_side_curve_points3d(i,0);
            Ay = pressure_side_curve_points3d(i,1);
            Az = pressure_side_curve_points3d(i,2);
            frac =(c*(Sz-z0)*(Sz-z0))/(-c*(Az-Sz)*(Sz-z0)+math::sqrt((Ax*Ax+Ay*Ay)*(Sz-z0)*(Sz-z0)));
            pressure_side_curve_points3d_bottom(i,0) = Ax*frac;
            pressure_side_curve_points3d_bottom(i,1) =  Ay*frac;
            pressure_side_curve_points3d_bottom(i,2) =  Sz+(Az-Sz)*frac;

        }



        gsBSpline<T> pressureSide3d_top;
        gsBSpline<T> suctionSide3d_top;
        gsBSpline<T> pressureSide3d_bottom;
        gsBSpline<T> suctionSide3d_bottom;
        gsMatrix<T> parameter_points_suction(1, num_samples3d+1);
        gsMatrix<T> parameter_points_pressure(1, num_samples3d+1);

        parameter_points_suction = centripetalParameterization(suction_side_curve_points3d_top);
        suctionSide3d_top = curveFittingWithBoundary(suction_side_curve_points3d_top, parameter_points_suction, kvfit);
        parameter_points_suction= centripetalParameterization(suction_side_curve_points3d_bottom);
        suctionSide3d_bottom = curveFittingWithBoundary(suction_side_curve_points3d_bottom, parameter_points_suction, kvfit);

        parameter_points_pressure = centripetalParameterization(pressure_side_curve_points3d_top);
        pressureSide3d_top = curveFittingWithBoundary(pressure_side_curve_points3d_top, parameter_points_pressure, kvfit);
        parameter_points_pressure= centripetalParameterization(pressure_side_curve_points3d_bottom);
        pressureSide3d_bottom = curveFittingWithBoundary(pressure_side_curve_points3d_bottom, parameter_points_pressure, kvfit);



        gsInfo<<"top and bottom curves computed \n";

        gsInfo<< "suction top \n" <<suctionSide3d_bottom.coefs() << "\n";
        gsInfo<< "pressure top \n" <<pressureSide3d_bottom.coefs() << "\n";
        gsInfo<< "suction bottom \n" <<suctionSide3d_top.coefs() << "\n";
        gsInfo<< "pressure bottom \n" <<pressureSide3d_top.coefs() << "\n";


       gsMatrix<T> coefs_suction((kvfit.size()-4) * 2, 3);
       coefs_suction.block(0,0,kvfit.size()-4,3) = suctionSide3d_bottom.coefs();
       coefs_suction.block(kvfit.size()-4,0,kvfit.size()-4,3) = suctionSide3d_top.coefs();

       real_t incline =myCircularMeshSetting(4,0);
       gsMatrix<T> coefs_suction_rot((kvfit.size()-4) * 2, 3);
          for (unsigned i = 0; i < (kvfit.size()-4) * 2; i++) {


           Ax = coefs_suction(i,0);
           Ay = coefs_suction(i,1);
           Az = coefs_suction(i,2);

           coefs_suction_rot(i,0)=Ax*math::sin( incline)*math::sin( incline) + math::cos(incline)*(-Ay*math::sin(myAngle)+(Az-Sz)*math::sin(incline)) + math::cos(incline)*math::cos(myAngle)*(Ax*math::cos(incline)+(-Az+Sz)*math::sin(incline));
           coefs_suction_rot(i,1)=Ay*math::cos(myAngle)+math::sin(myAngle)*(Ax*math::cos(incline)+(-Az+Sz)*math::sin(incline));
            coefs_suction_rot(i,2)=Sz+(Az-Sz)*math::cos( incline)*math::cos( incline) + Ay*math::sin(incline)*math::sin(myAngle) + Az*math::cos(myAngle)*math::sin(incline)*math::sin(incline)-Sz*math::cos(myAngle)*math::sin(incline)*math::sin(incline)+Ax*math::sin(myAngle/2)*math::sin(myAngle/2)*math::sin(2*incline);

       }


      // coefs_suction.block(0,0,kvfit.size()-4,3) = poefs.block(0,0,kvfit.size()-4,3);
      // coefs_suction.block(kvfit.size()-4,0,kvfit.size()-4,3) = poefs.block(kvfit.size()-4,0,kvfit.size()-4,3);
       gsTensorBSpline<2, T>  suction_side_surface(basis, coefs_suction_rot);
       suction_side_surface.degreeElevate(2,1);

       gsMatrix<T> coefs_pressure((kvfit.size()-4) * 2, 3);
       coefs_pressure.block(0,0,kvfit.size()-4,3) = pressureSide3d_bottom.coefs();
       coefs_pressure.block(kvfit.size()-4,0,kvfit.size()-4,3) = pressureSide3d_top.coefs();
   //    coefs_pressure.block(0,0,kvfit.size()-4,3) = poefs.block(2*kvfit.size()-8,0,kvfit.size()-4,3);
     //  coefs_pressure.block(kvfit.size()-4,0,kvfit.size()-4,3) = poefs.block(3*kvfit.size()-12,0,kvfit.size()-4,3);

       gsMatrix<T> coefs_pressure_rot((kvfit.size()-4) * 2, 3);
          for (unsigned i = 0; i < (kvfit.size()-4) * 2; i++) {


           Ax = coefs_pressure(i,0);
           Ay = coefs_pressure(i,1);
           Az = coefs_pressure(i,2);

           coefs_pressure_rot(i,0)=Ax*math::sin( incline)*math::sin( incline) + math::cos(incline)*(-Ay*math::sin(myAngle)+(Az-Sz)*math::sin(incline)) + math::cos(incline)*math::cos(myAngle)*(Ax*math::cos(incline)+(-Az+Sz)*math::sin(incline));
           coefs_pressure_rot(i,1)=Ay*math::cos(myAngle)+math::sin(myAngle)*(Ax*math::cos(incline)+(-Az+Sz)*math::sin(incline));
            coefs_pressure_rot(i,2)=Sz+(Az-Sz)*math::cos( incline)*math::cos( incline) + Ay*math::sin(incline)*math::sin(myAngle) + Az*math::cos(myAngle)*math::sin(incline)*math::sin(incline)-Sz*math::cos(myAngle)*math::sin(incline)*math::sin(incline)+Ax*math::sin(myAngle/2)*math::sin(myAngle/2)*math::sin(2*incline);

       }
       gsTensorBSpline<2, T>  pressure_side_surface(basis, coefs_pressure_rot);
       pressure_side_surface.degreeElevate(2,1);

       gsInfo << "surfaces computed \n";

       gsInfo << "suction_side_surface \n" << suction_side_surface << "\n";
       gsInfo << "pressure_side_surface \n" << pressure_side_surface << "\n";

       gsInfo << "suction_side_surface mesh \n" << suction_side_surface.coefs() << "\n";
       gsInfo << "pressure_side_surface mesh \n" << pressure_side_surface.coefs() << "\n";


       mySuctionSideSurface = suction_side_surface;
       myPressureSideSurface = pressure_side_surface;


       if (plot) {
           gsWriteParaview( suction_side_surface, "suctionSideSurfaceGV", 5000);
           gsWriteParaview( pressure_side_surface, "pressureSideSurfaceGV", 5000);

           gsMesh<> suction_side_mesh;
           gsMesh<> pressure_side_mesh;
           suction_side_surface.controlNet(suction_side_mesh);
           pressure_side_surface.controlNet(pressure_side_mesh);
           gsWriteParaview( suction_side_mesh, "suctionSideMeshGV");
           gsWriteParaview( pressure_side_mesh, "pressureSideMeshGV");
       }

       return 0;

   }



   template <class T>
   int KaplanTurbineGuideVane<T>::trimming(bool plot) {


       int num_samples3d = 50;


       gsKnotVector<T> kvfit(0, 1, 4, 4);
         gsKnotVector<T> kv(0, 1, 0, 2);
         gsTensorBSplineBasis<2, T> basis(kvfit,kv);




        gsMatrix<T> parameter_points_bottom(2, num_samples3d+1);
         gsMatrix<T> parameter_points_top(2, num_samples3d+1);
        gsMatrix<T> suction_side_curve_points_top(num_samples3d+1, 3);
        suction_side_curve_points_top.setZero(num_samples3d+1, 3);
        gsMatrix<T> suction_side_curve_points_bottom(num_samples3d+1, 3);
        suction_side_curve_points_bottom.setZero(num_samples3d+1, 3);
        gsMatrix<T> pressure_side_curve_points_top(num_samples3d+1, 3);
        pressure_side_curve_points_top.setZero(num_samples3d+1, 3);
        gsMatrix<T> pressure_side_curve_points_bottom(num_samples3d+1, 3);
        pressure_side_curve_points_bottom.setZero(num_samples3d+1, 3);
        gsMatrix<T> suction_side_curve_points3d_top(num_samples3d+1, 3);
     suction_side_curve_points3d_top.setZero(num_samples3d+1, 3);
        gsMatrix<T> pressure_side_curve_points3d_top(num_samples3d+1, 3);
        pressure_side_curve_points3d_top.setZero(num_samples3d+1, 3);
        gsMatrix<T> suction_side_curve_points3d_bottom(num_samples3d+1, 3);
         suction_side_curve_points3d_bottom.setZero(num_samples3d+1, 3);
        gsMatrix<T> pressure_side_curve_points3d_bottom(num_samples3d+1, 3);
        pressure_side_curve_points3d_bottom.setZero(num_samples3d+1, 3);

        for (index_t i = 0; i < num_samples3d; i++) {
            //parameter_points(0,i) = (1-cos(pi/2*i/num_samples3d));
            parameter_points_bottom(0,i) = pow((real_t) i/num_samples3d,2);
            parameter_points_bottom(1,i) = 0;
        }

        parameter_points_bottom(0,num_samples3d) = 1;
        parameter_points_bottom(1,num_samples3d) = 0;

        for (index_t i = 0; i < num_samples3d; i++) {
            //parameter_points(0,i) = (1-cos(pi/2*i/num_samples3d));
            parameter_points_top(0,i) = pow((real_t) i/num_samples3d,2);
            parameter_points_top(1,i) = 1;
        }

        parameter_points_top(0,num_samples3d) = 1;
        parameter_points_top(1,num_samples3d) = 0;


        gsInfo<<"---------------------------------\n";
        pressure_side_curve_points_bottom = myPressureSideSurface.eval(parameter_points_bottom);
        pressure_side_curve_points_top = myPressureSideSurface.eval(parameter_points_top);
        suction_side_curve_points_bottom = mySuctionSideSurface.eval(parameter_points_bottom);
        suction_side_curve_points_top = mySuctionSideSurface.eval(parameter_points_top);
  gsInfo<<"Po eval---------------------------------\n";

        suction_side_curve_points3d_bottom = suction_side_curve_points_bottom.transpose();
        pressure_side_curve_points3d_bottom =  pressure_side_curve_points_bottom.transpose();
        suction_side_curve_points3d_top = suction_side_curve_points_top.transpose();
        pressure_side_curve_points3d_top =  pressure_side_curve_points_top.transpose();


        // top profile on the cone
        real_t Sz = myCircularMeshSetting(5,0);




        real_t Ax,Ay,Az,frac;
        real_t R = myCircularMeshSetting(6,0);

        for (index_t i = 0; i < num_samples3d+1; i++) {


            Ax = suction_side_curve_points3d_top(i,0);
            Ay = suction_side_curve_points3d_top(i,1);
            Az = suction_side_curve_points3d_top(i,2);
             frac =math::sqrt(Ax*Ax+Ay*Ay+(Az-Sz)*(Az-Sz));
            suction_side_curve_points3d_top(i,0) = Ax*R/frac;
            suction_side_curve_points3d_top(i,1) =  Ay*R/frac;
           suction_side_curve_points3d_top(i,2) =  Sz+R*(Az-Sz)/frac;

            Ax = pressure_side_curve_points3d_top(i,0);
            Ay = pressure_side_curve_points3d_top(i,1);
            Az = pressure_side_curve_points3d_top(i,2);
             frac =math::sqrt(Ax*Ax+Ay*Ay+(Az-Sz)*(Az-Sz));
            pressure_side_curve_points3d_top(i,0) = Ax*R/frac;
            pressure_side_curve_points3d_top(i,1) =  Ay*R/frac;
            pressure_side_curve_points3d_top(i,2) =  Sz+(Az-Sz)*R/frac;

        }






          R = myCircularMeshSetting(7,0);


                  for (index_t i = 0; i < num_samples3d+1; i++) {


                      Ax = suction_side_curve_points3d_bottom(i,0);
                      Ay = suction_side_curve_points3d_bottom(i,1);
                      Az = suction_side_curve_points3d_bottom(i,2);
                       frac =math::sqrt(Ax*Ax+Ay*Ay+(Az-Sz)*(Az-Sz));
                      suction_side_curve_points3d_bottom(i,0) = Ax*R/frac;
                      suction_side_curve_points3d_bottom(i,1) =  Ay*R/frac;
                     suction_side_curve_points3d_bottom(i,2) =  Sz+(Az-Sz)*R/frac;

                      Ax = pressure_side_curve_points3d_bottom(i,0);
                      Ay = pressure_side_curve_points3d_bottom(i,1);
                      Az = pressure_side_curve_points3d_bottom(i,2);
                     frac =math::sqrt(Ax*Ax+Ay*Ay+(Az-Sz)*(Az-Sz));
                      pressure_side_curve_points3d_bottom(i,0) = Ax*R/frac;
                      pressure_side_curve_points3d_bottom(i,1) =  Ay*R/frac;
                      pressure_side_curve_points3d_bottom(i,2) =  Sz+(Az-Sz)*R/frac;

                  }



        gsBSpline<T> pressureSide3d_top;
        gsBSpline<T> suctionSide3d_top;
        gsBSpline<T> pressureSide3d_bottom;
        gsBSpline<T> suctionSide3d_bottom;
        gsMatrix<T> parameter_points_suction(1, num_samples3d+1);
        gsMatrix<T> parameter_points_pressure(1, num_samples3d+1);

        parameter_points_suction = centripetalParameterization(suction_side_curve_points3d_top);
        suctionSide3d_top = curveFittingWithBoundary(suction_side_curve_points3d_top, parameter_points_suction, kvfit);
        parameter_points_suction= centripetalParameterization(suction_side_curve_points3d_bottom);
        suctionSide3d_bottom = curveFittingWithBoundary(suction_side_curve_points3d_bottom, parameter_points_suction, kvfit);

        parameter_points_pressure = centripetalParameterization(pressure_side_curve_points3d_top);
        pressureSide3d_top = curveFittingWithBoundary(pressure_side_curve_points3d_top, parameter_points_pressure, kvfit);
        parameter_points_pressure= centripetalParameterization(pressure_side_curve_points3d_bottom);
        pressureSide3d_bottom = curveFittingWithBoundary(pressure_side_curve_points3d_bottom, parameter_points_pressure, kvfit);



        gsInfo<<"top and bottom curves computed \n";

        gsInfo<< "suction top \n" <<suctionSide3d_bottom.coefs() << "\n";
        gsInfo<< "pressure top \n" <<pressureSide3d_bottom.coefs() << "\n";
        gsInfo<< "suction bottom \n" <<suctionSide3d_top.coefs() << "\n";
        gsInfo<< "pressure bottom \n" <<pressureSide3d_top.coefs() << "\n";


       gsMatrix<T> coefs_suction((kvfit.size()-4) * 2, 3);
       coefs_suction.block(0,0,kvfit.size()-4,3) = suctionSide3d_bottom.coefs();
       coefs_suction.block(kvfit.size()-4,0,kvfit.size()-4,3) = suctionSide3d_top.coefs();




      // coefs_suction.block(0,0,kvfit.size()-4,3) = poefs.block(0,0,kvfit.size()-4,3);
      // coefs_suction.block(kvfit.size()-4,0,kvfit.size()-4,3) = poefs.block(kvfit.size()-4,0,kvfit.size()-4,3);
       gsTensorBSpline<2, T>  suction_side_surface(basis, coefs_suction);
       suction_side_surface.degreeElevate(2,1);

       gsMatrix<T> coefs_pressure((kvfit.size()-4) * 2, 3);
       coefs_pressure.block(0,0,kvfit.size()-4,3) = pressureSide3d_bottom.coefs();
       coefs_pressure.block(kvfit.size()-4,0,kvfit.size()-4,3) = pressureSide3d_top.coefs();
   //    coefs_pressure.block(0,0,kvfit.size()-4,3) = poefs.block(2*kvfit.size()-8,0,kvfit.size()-4,3);
     //  coefs_pressure.block(kvfit.size()-4,0,kvfit.size()-4,3) = poefs.block(3*kvfit.size()-12,0,kvfit.size()-4,3);


       gsTensorBSpline<2, T>  pressure_side_surface(basis, coefs_pressure);
       pressure_side_surface.degreeElevate(2,1);

       gsInfo << "surfaces computed \n";

       gsInfo << "suction_side_surface \n" << suction_side_surface << "\n";
       gsInfo << "pressure_side_surface \n" << pressure_side_surface << "\n";

       gsInfo << "suction_side_surface mesh \n" << suction_side_surface.coefs() << "\n";
       gsInfo << "pressure_side_surface mesh \n" << pressure_side_surface.coefs() << "\n";


       mySuctionSideSurfaceAfterTrimm = suction_side_surface;
       myPressureSideSurfaceAfterTrimm = pressure_side_surface;



       if (plot) {
           gsWriteParaview( suction_side_surface, "suctionSideSurfaceAfterTrimmGV", 5000);
           gsWriteParaview( pressure_side_surface, "pressureSideSurfaceAfterTrimmGV", 5000);

           gsMesh<> suction_side_mesh;
           gsMesh<> pressure_side_mesh;
           suction_side_surface.controlNet(suction_side_mesh);
           pressure_side_surface.controlNet(pressure_side_mesh);
           gsWriteParaview( suction_side_mesh, "suctionSideMeshAfterTrimmGV");
           gsWriteParaview( pressure_side_mesh, "pressureSideMeshAfterTrimmGV");
       }


       return 0;

   }

#endif // UWBKAPLANTURBINEGUIDEVANE
