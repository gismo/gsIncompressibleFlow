#ifndef UWBKAPLANTURBINERUNNERBLADE_H
#define UWBKAPLANTURBINERUNNERBLADE_H

using namespace gismo;

template <class T>
class KaplanTurbineRunnerBlade
{
public:
    // Constructors and destructors
    KaplanTurbineRunnerBlade() {}
    KaplanTurbineRunnerBlade(int num_bladeprofiles, gsVector<T> thisCamberX, gsVector<T> thisCamberY, gsVector<T> thisLeadingAngle, gsVector<T> thisTrailingAngle,
                             gsVector<T> thisThicknessX, gsVector<T> thisThicknessY, gsVector<T> thisEndingOffset,
                             gsVector<T> thisOutputAngle, gsVector<T> thisChordLength, gsVector<T> thisAngle,
                             gsVector<T> thisRotationCenterX, gsVector<T> thisRotationCenterY, gsVector<T> thisCylinderRadius);
    KaplanTurbineRunnerBlade(int num_bladeprofiles, gsVector<T> thisCamberX, gsVector<T> thisCamberY, gsVector<T> thisLeadingAngle, gsVector<T> thisTrailingAngle,
                             gsVector<T> thisThicknessX, gsVector<T> thisThicknessY, gsVector<T> thisEndingOffset,
                             gsVector<T> thisOutputAngle, gsVector<T> thisRadius, gsVector<T> thisChordLength, gsVector<T> thisAngle,
                             gsVector<T> thisRotationCenterX, gsVector<T> thisRotationCenterY, gsVector<T> thisCylinderRadius);
    KaplanTurbineRunnerBlade(int num_bladeprofiles, gsVector<T> thisCamberX, gsVector<T> thisCamberY, gsVector<T> thisLeadingAngle, gsVector<T> thisTrailingAngle,
                             gsVector<T> thisThicknessX, gsVector<T> thisThicknessY, gsVector<T> thisEndingOffset,
                             gsVector<T> thisOutputAngle, gsVector<T> thisRadius, gsVector<T> thisChordLength, gsVector<T> thisAngle,
                             gsVector<T> thisRotationCenterX, gsVector<T> thisRotationCenterY, gsVector<T> thisCylinderRadius, bool thisExtrapolate);
    ~KaplanTurbineRunnerBlade() {
        for (index_t i = 0; i < m_num_bladeprofiles; i++) {
            delete m_BladeProfiles[i];
        }
    }

    // Access function
    void setSuctionSideSurface(gsTensorBSpline<2,T> bspl) { mySuctionSideSurface = bspl; }
    gsTensorBSpline<2,T> getSuctionSideSurface() { return mySuctionSideSurface; }
    void setPressureSideSurface(gsTensorBSpline<2,T> bspl) { myPressureSideSurface = bspl; }
    gsTensorBSpline<2,T> getPressureSideSurface() { return myPressureSideSurface; }
    void setSuctionSideSurfaceAfterTrimming(gsTensorBSpline<2,T> bspl) { mySuctionSideSurfaceAfterTrimming = bspl; }
    gsTensorBSpline<2,T> getSuctionSideSurfaceAfterTrimming() { return mySuctionSideSurfaceAfterTrimming; }
    void setPressureSideSurfaceAfterTrimming(gsTensorBSpline<2,T> bspl) { myPressureSideSurfaceAfterTrimming = bspl; }
    gsTensorBSpline<2,T> getPressureSideSurfaceAfterTrimming() { return myPressureSideSurfaceAfterTrimming; }
    void setNumBladeProfiles(int num) { m_num_bladeprofiles = num; }
    int getNumBladeProfiles() { return m_num_bladeprofiles; }
    BladeProfile<T> * getBladeProfile(int i) { return m_BladeProfiles[i]; }
    gsMatrix<T> getSectionParameters() { return mySectionParameters; }
    gsMatrix<T> getSuctionIntersectionPointinParameterDomainonInnerHR() { return mySuctionIntersectionPointinParameterDomainonInnerHR; }
    gsMatrix<T> getPressureIntersectionPointinParameterDomainonInnerHR() { return myPressureIntersectionPointinParameterDomainonInnerHR; }
    gsVector<T> getCylinderRadius() { return myCylinderRadius; }
    bool getExtrapolate() { return myExtrapolate; }

    // Other functions
    int compute(bool plot, bool print_info = false );
    int trimming(gsTensorBSpline<2,T> trim_surface, bool plot, bool print_info = false);
    int exportCurveFile(std::string const & filename, bool plot, bool print_info);

private:

    int m_num_bladeprofiles;

    gsVector<T> myCamberX;
    gsVector<T> myCamberY;
    gsVector<T> myLeadingAngle;
    gsVector<T> myTrailingAngle;
    gsVector<T> myThicknessX;
    gsVector<T> myThicknessY;
    gsVector<T> myEndingOffset;
    gsVector<T> myOutputAngle;
    gsVector<T> myRadius;
    gsVector<T> myChordLength;
    gsVector<T> myAngle;
    gsVector<T> myRotationCenterX;
    gsVector<T> myRotationCenterY;
    gsVector<T> myCylinderRadius;
    bool myExtrapolate;

    gsMatrix<T> mySectionParameters;

    gsMatrix<T> mySuctionIntersectionPointinParameterDomainonInnerHR;
    gsMatrix<T> myPressureIntersectionPointinParameterDomainonInnerHR;
    gsMatrix<T> mySuctionIntersectionPointonInnerHR;
    gsMatrix<T> myPressureIntersectionPointonInnerHR;

    gsTensorBSpline<2,T> myPressureSideSurface;
    gsTensorBSpline<2,T> mySuctionSideSurface;

    gsTensorBSpline<2,T> mySuctionSideSurfaceAfterTrimming;
    gsTensorBSpline<2,T> myPressureSideSurfaceAfterTrimming;

    BladeProfile<T> * m_BladeProfiles[20];

};

template <class T>
KaplanTurbineRunnerBlade<T>::KaplanTurbineRunnerBlade(int num_bladeprofiles, gsVector<T> thisCamberX, gsVector<T> thisCamberY, gsVector<T> thisLeadingAngle, gsVector<T> thisTrailingAngle,
                                                   gsVector<T> thisThicknessX, gsVector<T> thisThicknessY, gsVector<T> thisEndingOffset,
                                                   gsVector<T> thisOutputAngle, gsVector<T> thisChordLength, gsVector<T> thisAngle,
                                                   gsVector<T> thisRotationCenterX, gsVector<T> thisRotationCenterY, gsVector<T> thisCylinderRadius) {

    myExtrapolate = true;
    T added_value;
    gsVector<T> ma(num_bladeprofiles+1);

    ma(0) = thisCylinderRadius(0) - 2*(thisCylinderRadius(1)-thisCylinderRadius(0))/3;
    //ma(0) = added_cylinder_radius;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisCylinderRadius(i); }
    myCylinderRadius = ma;
    extrapolateData(thisCylinderRadius, thisCamberX, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisCamberX(i); }
    myCamberX = ma;
    extrapolateData(thisCylinderRadius, thisCamberY, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisCamberY(i); }
    myCamberY = ma;
    extrapolateData(thisCylinderRadius, thisLeadingAngle, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisLeadingAngle(i); }
    myLeadingAngle = ma;
    extrapolateData(thisCylinderRadius, thisTrailingAngle, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisTrailingAngle(i); }
    myTrailingAngle = ma;
    extrapolateData(thisCylinderRadius, thisThicknessX, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisThicknessX(i); }
    myThicknessX = ma;
    extrapolateData(thisCylinderRadius, thisThicknessY, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisThicknessY(i); }
    myThicknessY = ma;
    extrapolateData(thisCylinderRadius, thisEndingOffset, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisEndingOffset(i); }
    myEndingOffset = ma;
    extrapolateData(thisCylinderRadius, thisOutputAngle, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisOutputAngle(i); }
    myOutputAngle = ma;
    extrapolateData(thisCylinderRadius, thisChordLength, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisChordLength(i); }
    myChordLength = ma;
    extrapolateData(thisCylinderRadius, thisAngle, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisAngle(i); }
    myAngle = ma;
    extrapolateData(thisCylinderRadius, thisRotationCenterX, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisRotationCenterX(i); }
    myRotationCenterX = ma;
    extrapolateData(thisCylinderRadius, thisRotationCenterY, myCylinderRadius(0), added_value, 1e-5, 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisRotationCenterY(i); }
    myRotationCenterY = ma;

    gsInfo << "Checking values for runner blade:\n";
    gsInfo << "=================================\n";
    gsInfo << "Cylinder radii: " << myCylinderRadius.asRowVector() << "\n";
    gsInfo << "Camber x: " << myCamberX.asRowVector() << "\n";
    gsInfo << "Camber y: " << myCamberY.asRowVector() << "\n";
    gsInfo << "Leading angle: " << myLeadingAngle.asRowVector() << "\n";
    gsInfo << "Trailing angle: " << myTrailingAngle.asRowVector() << "\n";
    gsInfo << "Thickness x: " << myThicknessX.asRowVector() << "\n";
    gsInfo << "Thickness y: " << myThicknessY.asRowVector() << "\n";
    gsInfo << "Ending offset: " << myEndingOffset.asRowVector() << "\n";
    gsInfo << "Output angle: " << myOutputAngle.asRowVector() << "\n";
    gsInfo << "Chord length: " << myChordLength.asRowVector() << "\n";
    gsInfo << "Angle: " << myAngle.asRowVector() << "\n";
    gsInfo << "Rotation center x: " << myRotationCenterX.asRowVector() << "\n";
    gsInfo << "Rotation center y: " << myRotationCenterY.asRowVector() << "\n";

    m_num_bladeprofiles = num_bladeprofiles + 1;

}

template <class T>
KaplanTurbineRunnerBlade<T>::KaplanTurbineRunnerBlade(int num_bladeprofiles, gsVector<T> thisCamberX, gsVector<T> thisCamberY, gsVector<T> thisLeadingAngle, gsVector<T> thisTrailingAngle,
                                                   gsVector<T> thisThicknessX, gsVector<T> thisThicknessY, gsVector<T> thisEndingOffset,
                                                   gsVector<T> thisOutputAngle, gsVector<T> thisRadius, gsVector<T> thisChordLength, gsVector<T> thisAngle,
                                                   gsVector<T> thisRotationCenterX, gsVector<T> thisRotationCenterY, gsVector<T> thisCylinderRadius) {
    myExtrapolate = true;
    T added_value;
    gsVector<T> ma(num_bladeprofiles+1);

    ma(0) = thisCylinderRadius(0) - 2*(thisCylinderRadius(1)-thisCylinderRadius(0))/3;
    //ma(0) = added_cylinder_radius;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisCylinderRadius(i); }
    myCylinderRadius = ma;
    extrapolateData(thisCylinderRadius, thisCamberX, myCylinderRadius(0), added_value, static_cast<real_t>(1e-2), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisCamberX(i); }
    myCamberX = ma;
    extrapolateData(thisCylinderRadius, thisCamberY, myCylinderRadius(0), added_value, static_cast<real_t>(1e-4), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisCamberY(i); }
    myCamberY = ma;
    extrapolateData(thisCylinderRadius, thisLeadingAngle, myCylinderRadius(0), added_value, static_cast<real_t>(1e-2), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisLeadingAngle(i); }
    myLeadingAngle = ma;
    extrapolateData(thisCylinderRadius, thisTrailingAngle, myCylinderRadius(0), added_value, static_cast<real_t>(1e-2), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisTrailingAngle(i); }
    myTrailingAngle = ma;
    extrapolateData(thisCylinderRadius, thisThicknessX, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisThicknessX(i); }
    myThicknessX = ma;
    extrapolateData(thisCylinderRadius, thisThicknessY, myCylinderRadius(0), added_value, static_cast<real_t>(1e-4), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisThicknessY(i); }
    myThicknessY = ma;
    extrapolateData(thisCylinderRadius, thisEndingOffset, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisEndingOffset(i); }
    myEndingOffset = ma;
    extrapolateData(thisCylinderRadius, thisOutputAngle, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisOutputAngle(i); }
    myOutputAngle = ma;
    extrapolateData(thisCylinderRadius, thisChordLength, myCylinderRadius(0), added_value, static_cast<real_t>(1e-4), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisChordLength(i); }
    myChordLength = ma;
    extrapolateData(thisCylinderRadius, thisAngle, myCylinderRadius(0), added_value, static_cast<real_t>(1e-2), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisAngle(i); }
    myAngle = ma;
    extrapolateData(thisCylinderRadius, thisRotationCenterX, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisRotationCenterX(i); }
    myRotationCenterX = ma;
    extrapolateData(thisCylinderRadius, thisRotationCenterY, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisRotationCenterY(i); }
    myRotationCenterY = ma;
    extrapolateData(thisCylinderRadius, thisRadius, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
    ma(0) = added_value;
    for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisRadius(i); }
    myRadius = ma;

    gsInfo << "Checking values for runner blade:\n";
    gsInfo << "=================================\n";
    gsInfo << "Cylinder radii: " << myCylinderRadius.asRowVector() << "\n";
    gsInfo << "Camber x: " << myCamberX.asRowVector() << "\n";
    gsInfo << "Camber y: " << myCamberY.asRowVector() << "\n";
    gsInfo << "Leading angle: " << myLeadingAngle.asRowVector() << "\n";
    gsInfo << "Trailing angle: " << myTrailingAngle.asRowVector() << "\n";
    gsInfo << "Thickness x: " << myThicknessX.asRowVector() << "\n";
    gsInfo << "Thickness y: " << myThicknessY.asRowVector() << "\n";
    gsInfo << "Ending offset: " << myEndingOffset.asRowVector() << "\n";
    gsInfo << "Output angle: " << myOutputAngle.asRowVector() << "\n";
    gsInfo << "Osculating radius: " << myRadius.asRowVector() << "\n";
    gsInfo << "Chord length: " << myChordLength.asRowVector() << "\n";
    gsInfo << "Angle: " << myAngle.asRowVector() << "\n";
    gsInfo << "Rotation center x: " << myRotationCenterX.asRowVector() << "\n";
    gsInfo << "Rotation center y: " << myRotationCenterY.asRowVector() << "\n";

    m_num_bladeprofiles = num_bladeprofiles + 1;

}

template <class T>
KaplanTurbineRunnerBlade<T>::KaplanTurbineRunnerBlade(int num_bladeprofiles, gsVector<T> thisCamberX, gsVector<T> thisCamberY, gsVector<T> thisLeadingAngle, gsVector<T> thisTrailingAngle,
                                                   gsVector<T> thisThicknessX, gsVector<T> thisThicknessY, gsVector<T> thisEndingOffset,
                                                   gsVector<T> thisOutputAngle, gsVector<T> thisRadius, gsVector<T> thisChordLength, gsVector<T> thisAngle,
                                                   gsVector<T> thisRotationCenterX, gsVector<T> thisRotationCenterY, gsVector<T> thisCylinderRadius, bool thisExtrapolate) {

        myCylinderRadius = thisCylinderRadius;

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

        myRadius = thisRadius;

        myExtrapolate = thisExtrapolate;


     gsInfo << "Checking values for runner blade:\n";
     gsInfo << "=================================\n";
     gsInfo << "Cylinder radii: " << myCylinderRadius.asRowVector() << "\n";
     gsInfo << "Camber x: " << myCamberX.asRowVector() << "\n";
     gsInfo << "Camber y: " << myCamberY.asRowVector() << "\n";
     gsInfo << "Leading angle: " << myLeadingAngle.asRowVector() << "\n";
     gsInfo << "Trailing angle: " << myTrailingAngle.asRowVector() << "\n";
     gsInfo << "Thickness x: " << myThicknessX.asRowVector() << "\n";
     gsInfo << "Thickness y: " << myThicknessY.asRowVector() << "\n";
     gsInfo << "Ending offset: " << myEndingOffset.asRowVector() << "\n";
     gsInfo << "Output angle: " << myOutputAngle.asRowVector() << "\n";
     gsInfo << "Chord length: " << myChordLength.asRowVector() << "\n";
     gsInfo << "Angle: " << myAngle.asRowVector() << "\n";
     gsInfo << "Rotation center x: " << myRotationCenterX.asRowVector() << "\n";
     gsInfo << "Rotation center y: " << myRotationCenterY.asRowVector() << "\n";

     m_num_bladeprofiles = num_bladeprofiles;

     if(myExtrapolate)
     {
         T added_value;
         gsVector<T> ma(num_bladeprofiles+1);

         ma(0) = thisCylinderRadius(0) - 2*(thisCylinderRadius(1)-thisCylinderRadius(0))/3;
         //ma(0) = added_cylinder_radius;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisCylinderRadius(i); }
         myCylinderRadius = ma;
         extrapolateData(thisCylinderRadius, thisCamberX, myCylinderRadius(0), added_value, static_cast<real_t>(1e-2), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisCamberX(i); }
         myCamberX = ma;
         extrapolateData(thisCylinderRadius, thisCamberY, myCylinderRadius(0), added_value, static_cast<real_t>(1e-4), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisCamberY(i); }
         myCamberY = ma;
         extrapolateData(thisCylinderRadius, thisLeadingAngle, myCylinderRadius(0), added_value, static_cast<real_t>(1e-2), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisLeadingAngle(i); }
         myLeadingAngle = ma;
         extrapolateData(thisCylinderRadius, thisTrailingAngle, myCylinderRadius(0), added_value, static_cast<real_t>(1e-2), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisTrailingAngle(i); }
         myTrailingAngle = ma;
         extrapolateData(thisCylinderRadius, thisThicknessX, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisThicknessX(i); }
         myThicknessX = ma;
         extrapolateData(thisCylinderRadius, thisThicknessY, myCylinderRadius(0), added_value, static_cast<real_t>(1e-4), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisThicknessY(i); }
         myThicknessY = ma;
         extrapolateData(thisCylinderRadius, thisEndingOffset, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisEndingOffset(i); }
         myEndingOffset = ma;
         extrapolateData(thisCylinderRadius, thisOutputAngle, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisOutputAngle(i); }
         myOutputAngle = ma;
         extrapolateData(thisCylinderRadius, thisChordLength, myCylinderRadius(0), added_value, static_cast<real_t>(1e-4), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisChordLength(i); }
         myChordLength = ma;
         extrapolateData(thisCylinderRadius, thisAngle, myCylinderRadius(0), added_value, static_cast<real_t>(1e-2), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisAngle(i); }
         myAngle = ma;
         extrapolateData(thisCylinderRadius, thisRotationCenterX, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisRotationCenterX(i); }
         myRotationCenterX = ma;
         extrapolateData(thisCylinderRadius, thisRotationCenterY, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisRotationCenterY(i); }
         myRotationCenterY = ma;
         extrapolateData(thisCylinderRadius, thisRadius, myCylinderRadius(0), added_value, static_cast<real_t>(1e-3), 20);
         ma(0) = added_value;
         for (index_t i = 0; i < num_bladeprofiles; i++) { ma(i+1) = thisRadius(i); }
         myRadius = ma;

         gsInfo << "Checking values for runner blade:\n";
         gsInfo << "=================================\n";
         gsInfo << "Cylinder radii: " << myCylinderRadius.asRowVector() << "\n";
         gsInfo << "Camber x: " << myCamberX.asRowVector() << "\n";
         gsInfo << "Camber y: " << myCamberY.asRowVector() << "\n";
         gsInfo << "Leading angle: " << myLeadingAngle.asRowVector() << "\n";
         gsInfo << "Trailing angle: " << myTrailingAngle.asRowVector() << "\n";
         gsInfo << "Thickness x: " << myThicknessX.asRowVector() << "\n";
         gsInfo << "Thickness y: " << myThicknessY.asRowVector() << "\n";
         gsInfo << "Ending offset: " << myEndingOffset.asRowVector() << "\n";
         gsInfo << "Output angle: " << myOutputAngle.asRowVector() << "\n";
         gsInfo << "Osculating radius: " << myRadius.asRowVector() << "\n";
         gsInfo << "Chord length: " << myChordLength.asRowVector() << "\n";
         gsInfo << "Angle: " << myAngle.asRowVector() << "\n";
         gsInfo << "Rotation center x: " << myRotationCenterX.asRowVector() << "\n";
         gsInfo << "Rotation center y: " << myRotationCenterY.asRowVector() << "\n";

         m_num_bladeprofiles = num_bladeprofiles + 1;
     }



}

template <class T>
int KaplanTurbineRunnerBlade<T>::compute(bool plot, bool print_info) {

    int num_samples = 30;
    int num_samples3d = 50;

    // 2D and 3D blade profiles computation
    gsInfo << "\n3D blade construction\n";
    gsInfo << "=====================\n";
    gsInfo << "Construction of 3D blade profiles from given shape parameters ...\n";
    //gsBSpline<T> suction_side_curve[this->m_num_bladeprofiles];
    //gsBSpline<T> pressure_side_curve[this->m_num_bladeprofiles];
    //gsBSpline<T> suction_side_curve3d[this->m_num_bladeprofiles];
    //gsBSpline<T> pressure_side_curve3d[this->m_num_bladeprofiles];
    gsBSpline<T> suction_curve;
    gsBSpline<T> pressure_curve;
    gsBSpline<T> suction_curve3d;
    gsBSpline<T> pressure_curve3d;
    std::vector<gsBSpline<T> > suction_side_curve;
    std::vector<gsBSpline<T> > pressure_side_curve;
    std::vector<gsBSpline<T> > suction_side_curve3d;
    std::vector<gsBSpline<T> > pressure_side_curve3d;
    gsKnotVector<T> kvfit(0, 1, 4, 4);
    //BladeProfile<real_t> * pBladeProfiles[num_bladeprofiles] = {};
    BladeProfile<T> * pBladeProfile;

    //gsInfo << myRadius.size() << "\n";
    for (int i = 0; i < m_num_bladeprofiles; i++) {
        //gsInfo << i << "\n";
        if (myRadius.size() == 0) {
            pBladeProfile = new BladeProfile<T>(myCamberX(i), myCamberY(i), myLeadingAngle(i), myTrailingAngle(i), myThicknessX(i),
                                                myThicknessY(i), myEndingOffset(i), myOutputAngle(i), myChordLength(i),
                                                myAngle(i), myRotationCenterX(i), myRotationCenterY(i), myCylinderRadius(i));
        }
        else {
            pBladeProfile = new BladeProfile<T>(myCamberX(i), myCamberY(i), myLeadingAngle(i), myTrailingAngle(i), myThicknessX(i),
                                                myThicknessY(i), myEndingOffset(i), myOutputAngle(i), myRadius(i), myChordLength(i),
                                                myAngle(i), myRotationCenterX(i), myRotationCenterY(i), myCylinderRadius(i));
        }
        //pBladeProfile->compute2D(suction_side_curve[i], pressure_side_curve[i], kvfit, num_samples, 3, 2);
        //pBladeProfile->compute3D(suction_side_curve3d[i], pressure_side_curve3d[i], kvfit, num_samples3d);
        pBladeProfile->compute2D(suction_curve, pressure_curve, kvfit, num_samples, 3, 2);
        suction_side_curve.push_back(suction_curve);
        pressure_side_curve.push_back(pressure_curve);
        pBladeProfile->compute3D(suction_curve3d, pressure_curve3d, kvfit, num_samples3d);
        suction_side_curve3d.push_back(suction_curve3d);
        pressure_side_curve3d.push_back(pressure_curve3d);
        m_BladeProfiles[i] = pBladeProfile;
    }
    gsInfo << "Done.\n\n";

    // Computation of lofting surface from 3D balde profiles for finding suction side and pressure surface
    // describing runner blade
    gsInfo << "Construction of suction and pressure side surface of runner blade from 3D blade profiles ...\n";
    gsKnotVector<T> kvloft(0, 1, m_num_bladeprofiles-4, 4);
    gsTensorBSpline<2,T> suction_side_surface;
    gsTensorBSpline<2,T> pressure_side_surface;
    gsMatrix<T> section_parameters(1, m_num_bladeprofiles);
    for (index_t i = 0; i < m_num_bladeprofiles; i++) {
        section_parameters(i) = (myCylinderRadius[i]-myCylinderRadius[0])/(myCylinderRadius[m_num_bladeprofiles-1]-myCylinderRadius[0]);
    }
    mySectionParameters = section_parameters;
    gsInfo << "Lofting parameters for 3D blade profiles: " << section_parameters << "\n";

    //gsBSpline<T> defcurves_suctionside[m_num_bladeprofiles];
    //gsBSpline<T> defcurves_pressureside[m_num_bladeprofiles];
    std::vector<gsBSpline<T> > defcurves_suctionside;
    std::vector<gsBSpline<T> > defcurves_pressureside;
    for (index_t i = 0; i < m_num_bladeprofiles; i++) {
        //defcurves_suctionside[i] = m_BladeProfiles[i]->getSuctionSide3d();
        //defcurves_pressureside[i] = m_BladeProfiles[i]->getPressureSide3d();
        defcurves_suctionside.push_back(m_BladeProfiles[i]->getSuctionSide3d());
        defcurves_pressureside.push_back(m_BladeProfiles[i]->getPressureSide3d());
    }
    computeLoftSurface(defcurves_suctionside, kvfit, m_num_bladeprofiles, kvloft, section_parameters, suction_side_surface);
    computeLoftSurface(defcurves_pressureside, kvfit, m_num_bladeprofiles, kvloft, section_parameters, pressure_side_surface);
    mySuctionSideSurface = suction_side_surface;
    myPressureSideSurface = pressure_side_surface;
    mySuctionSideSurfaceAfterTrimming = suction_side_surface;
    myPressureSideSurfaceAfterTrimming = pressure_side_surface;
    gsInfo << "Done.\n\n";

    if (plot) {
        std::vector<gsGeometry<>*> curves;
        curves.clear();
        std::vector<gsGeometry<>*> curves3d;
        curves3d.clear();
        for (int i = 0; i < m_num_bladeprofiles; i++) {
            curves.push_back(&suction_side_curve[i]);
            curves3d.push_back(&suction_side_curve3d[i]);
            curves.push_back(&pressure_side_curve[i]);
            curves3d.push_back(&pressure_side_curve3d[i]);
        }
        gsWriteParaview( curves, "2Dprofile", 100);
        gsWriteParaview( curves3d, "3Dprofile", 100);

        gsWriteParaview( suction_side_surface, "suctionSideSurface", 5000);
        gsWriteParaview( pressure_side_surface, "pressureSideSurface", 5000);

        gsMesh<> suction_side_mesh;
        gsMesh<> pressure_side_mesh;
        suction_side_surface.controlNet(suction_side_mesh);
        pressure_side_surface.controlNet(pressure_side_mesh);
        gsWriteParaview( suction_side_mesh, "suctionSideMesh");
        gsWriteParaview( pressure_side_mesh, "pressureSideMesh");
    }

    return 0;

}

template <class T>
int KaplanTurbineRunnerBlade<T>::trimming(gsTensorBSpline<2,T> trim_surface, bool plot, bool print_info) {

    int num_samples3d = 50;
    T par;
    int num_iter_s = 0;
    int max_num_iter_s = 0;
    int num_iter_p = 0;
    int max_num_iter_p = 0;
    gsVector<T> initial_solution(3);
    gsVector<T> solution(3);
    gsKnotVector<T> kvloft = mySuctionSideSurfaceAfterTrimming.knots(0);
    gsKnotVector<T> kvfit = mySuctionSideSurfaceAfterTrimming.knots(1);
    gsMatrix<T> section_parameters(1, m_num_bladeprofiles);
    for (index_t i = 0; i < m_num_bladeprofiles; i++) {
        section_parameters(i) = (myCylinderRadius[i]-myCylinderRadius[0])/(myCylinderRadius[m_num_bladeprofiles-1]-myCylinderRadius[0]);
    }

    gsInfo << "Trimming suction and pressure side surface by given surface\n";
    gsInfo << "===========================================================\n";

    // Computation of intersection of trim_surface with suction and pressure sides of the runner blade
    gsInfo << "Computing intersection points of suction side surface with trim surface ...\n";
    gsMatrix<T> suction_intersection_parameter_points_on_HR(2, num_samples3d);
    gsMatrix<T> suction_intersection_parameter_points_on_blade(2, num_samples3d);
    gsMatrix<T> suction_intersection_points_on_HR(3, num_samples3d);
    gsMatrix<T> suction_intersection_points_on_blade(3, num_samples3d);
    for (index_t i = 0; i < num_samples3d; i++) {
        par = (T) i/(num_samples3d-1);
        initial_solution << 0.6, 0.25, 0.1;
        curveSurfaceIntersectionviaNR(mySuctionSideSurfaceAfterTrimming, par, trim_surface, initial_solution, solution, static_cast<real_t>(1e-5), 100, num_iter_s, false);
        suction_intersection_parameter_points_on_HR(0,i) = solution(0);
        suction_intersection_parameter_points_on_HR(1,i) = solution(1);
        suction_intersection_parameter_points_on_blade(0,i) = solution(2);
        suction_intersection_parameter_points_on_blade(1,i) = par;
        if (max_num_iter_s < num_iter_s) {
            max_num_iter_s = num_iter_s;
        }
    }
    if (print_info) {
        gsInfo << suction_intersection_parameter_points_on_HR << "\n";
        gsInfo << "Max number of NR iteration for suction side: " << max_num_iter_s << "\n\n";
    }
    mySuctionIntersectionPointinParameterDomainonInnerHR = suction_intersection_parameter_points_on_HR;
    mySuctionSideSurfaceAfterTrimming.eval_into(suction_intersection_parameter_points_on_blade, suction_intersection_points_on_blade);
    trim_surface.eval_into(suction_intersection_parameter_points_on_HR, suction_intersection_points_on_HR);
    mySuctionIntersectionPointonInnerHR = suction_intersection_points_on_HR;
    gsInfo << "Done.\n";

    gsInfo << "Computing intersection points of pressure side surface with trim surface ...\n";
    gsMatrix<T> pressure_intersection_parameter_points_on_HR(2, num_samples3d);
    gsMatrix<T> pressure_intersection_parameter_points_on_blade(2, num_samples3d);
    gsMatrix<T> pressure_intersection_points_on_HR(3, num_samples3d);
    gsMatrix<T> pressure_intersection_points_on_blade(3, num_samples3d);
    for (index_t i = 0; i < num_samples3d; i++) {
        par = (T) i/(num_samples3d-1);
        initial_solution << 0.65, 0.25, 0.1;
        curveSurfaceIntersectionviaNR(myPressureSideSurfaceAfterTrimming, par, trim_surface, initial_solution, solution, static_cast<real_t>(1e-5), 100, num_iter_p, false);
        pressure_intersection_parameter_points_on_HR(0,i) = solution(0);
        pressure_intersection_parameter_points_on_HR(1,i) = solution(1);
        pressure_intersection_parameter_points_on_blade(0,i) = solution(2);
        pressure_intersection_parameter_points_on_blade(1,i) = par;
        if (max_num_iter_p < num_iter_p) {
            max_num_iter_p = num_iter_p;
        }
    }
    if (print_info) {
        gsInfo << pressure_intersection_parameter_points_on_HR << "\n";
        gsInfo << "Max number of NR iteration for pressure side: " << max_num_iter_p << "\n\n";
    }
    myPressureIntersectionPointinParameterDomainonInnerHR = pressure_intersection_parameter_points_on_HR;
    myPressureSideSurfaceAfterTrimming.eval_into(pressure_intersection_parameter_points_on_blade, pressure_intersection_points_on_blade);
    trim_surface.eval_into(pressure_intersection_parameter_points_on_HR, pressure_intersection_points_on_HR);
    myPressureIntersectionPointonInnerHR = pressure_intersection_points_on_HR;
    gsInfo << "Done.\n";

    if (plot) {
        gsWriteParaviewPoints<real_t>( suction_intersection_points_on_blade, "suction_intersection_points_on_blade");
        gsWriteParaviewPoints<real_t>( suction_intersection_points_on_HR, "suction_intersection_points_on_HR");
        gsWriteParaviewPoints<real_t>( pressure_intersection_points_on_blade, "pressure_intersection_points_on_blade");
        gsWriteParaviewPoints<real_t>( pressure_intersection_points_on_HR, "pressure_intersection_points_on_HR");
    }

    // Approximation of intersection points of trim surface with suction and pressure sides of the runner blade
    // by B-spline curve
    gsMatrix<T> suction_intersection_points_on_blade2(num_samples3d, 3);
    gsMatrix<T> pressure_intersection_points_on_blade2(num_samples3d, 3);
    gsMatrix<T> suction_intersection_leading_tangent(3,1);
    gsMatrix<T> pressure_intersection_leading_tangent(3,1);
    gsMatrix<T> suction_intersection_tangents = mySuctionSideSurfaceAfterTrimming.deriv(suction_intersection_parameter_points_on_blade);
    //gsInfo << "Rozmery suction_intersection_tangents: " << suction_intersection_tangents.rows() << ", " << suction_intersection_tangents.cols() << "\n";
    suction_intersection_leading_tangent(0) = suction_intersection_tangents(1,0);
    suction_intersection_leading_tangent(1) = suction_intersection_tangents(3,0);
    suction_intersection_leading_tangent(2) = suction_intersection_tangents(5,0);
    pressure_intersection_leading_tangent = -suction_intersection_leading_tangent;
    //gsInfo << suction_intersection_leading_tangent << "\n";
    //gsInfo << pressure_intersection_leading_tangent << "\n";

    gsMatrix<T> parameter_points_suction(1, num_samples3d);
    gsMatrix<T> parameter_points_pressure(1, num_samples3d);
    gsBSpline<T> suctionSideIntersection3d;
    gsBSpline<T> pressureSideIntersection3d;
    pressure_intersection_points_on_blade.col(0) = suction_intersection_points_on_blade.col(0);
    pressure_intersection_points_on_blade.col(num_samples3d-1) = suction_intersection_points_on_blade.col(num_samples3d-1);
    suction_intersection_points_on_blade2 = suction_intersection_points_on_blade.transpose();
    pressure_intersection_points_on_blade2 = pressure_intersection_points_on_blade.transpose();

    //gsInfo << "Rozmery suction_intersection_points_on_blade: " << suction_intersection_points_on_blade.rows() << ", " << suction_intersection_points_on_blade.cols() << "\n";
    //gsInfo << "Rozmery suction_intersection_points_on_blade2: " << suction_intersection_points_on_blade2.rows() << ", " << suction_intersection_points_on_blade2.cols() << "\n";
    //gsInfo << "Rozmery pressure_intersection_points_on_blade: " << pressure_intersection_points_on_blade.rows() << ", " << pressure_intersection_points_on_blade.cols() << "\n";
    //gsInfo << "Rozmery pressure_intersection_points_on_blade2: " << pressure_intersection_points_on_blade2.rows() << ", " << pressure_intersection_points_on_blade2.cols() << "\n";


    parameter_points_suction = centripetalParameterization(suction_intersection_points_on_blade2);
    suctionSideIntersection3d = curveFittingWithBoundaryAndInputTangent(suction_intersection_points_on_blade2, parameter_points_suction, kvfit, suction_intersection_leading_tangent);
    gsMatrix<T> control_points_suction_intersection = suctionSideIntersection3d.coefs();
    parameter_points_pressure = centripetalParameterization(pressure_intersection_points_on_blade2);
    pressureSideIntersection3d = curveFittingWithBoundaryAndInputTangent(pressure_intersection_points_on_blade2, parameter_points_pressure, kvfit, pressure_intersection_leading_tangent);
    gsMatrix<T> control_points_pressure_intersection = pressureSideIntersection3d.coefs();

    // Making intersection curve C1 at leading point
    T l1 = sqrt(pow(control_points_suction_intersection(1,0)-control_points_suction_intersection(0,0),2) + pow(control_points_suction_intersection(1,1)-control_points_suction_intersection(0,1),2));
    T l2 = sqrt(pow(control_points_pressure_intersection(1,0)-control_points_pressure_intersection(0,0),2) + pow(control_points_pressure_intersection(1,1)-control_points_pressure_intersection(0,1),2));
    //gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction_intersection(0,0) << "\t\t\t" << control_points_suction_intersection(0,1) << "\n" <<
    //          control_points_suction_intersection(1,0) << "\t\t" << control_points_suction_intersection(1,1) << "\n" <<
    //          control_points_pressure_intersection(0,0) << "\t\t\t" << control_points_pressure_intersection(0,1) << "\n" <<
    //          control_points_pressure_intersection(1,0) << "\t\t" << control_points_pressure_intersection(1,1) << "\n";
    if (l2 > l1) {
        control_points_suction_intersection(1,0) = control_points_suction_intersection(0,0) + l2*(control_points_suction_intersection(1,0) - control_points_suction_intersection(0,0))/l1;
        control_points_suction_intersection(1,1) = control_points_suction_intersection(0,1) + l2*(control_points_suction_intersection(1,1) - control_points_suction_intersection(0,1))/l1;
    //    gsInfo << "Upravene body:\n" << control_points_suction_intersection << "\n";
        suctionSideIntersection3d = gsBSpline<T>( kvfit, control_points_suction_intersection);
    }
    else {
        control_points_pressure_intersection(1,0) = control_points_pressure_intersection(0,0) + l1/l2*(control_points_pressure_intersection(1,0) - control_points_pressure_intersection(0,0));
        control_points_pressure_intersection(1,1) = control_points_pressure_intersection(0,1) + l1/l2*(control_points_pressure_intersection(1,1) - control_points_pressure_intersection(0,1));
    //    gsInfo << "Upravene body:\n" << control_points_pressure_intersection << "\n";
        pressureSideIntersection3d = gsBSpline<T>( kvfit, control_points_pressure_intersection);
    }

    l1 = sqrt(pow(control_points_suction_intersection(1,0)-control_points_suction_intersection(0,0),2) + pow(control_points_suction_intersection(1,1)-control_points_suction_intersection(0,1),2));
    l2 = sqrt(pow(control_points_pressure_intersection(1,0)-control_points_pressure_intersection(0,0),2) + pow(control_points_pressure_intersection(1,1)-control_points_pressure_intersection(0,1),2));
    //gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction_intersection(0,0) << "\t\t\t" << control_points_suction_intersection(0,1) << "\n" <<
    //         control_points_suction_intersection(1,0) << "\t\t" << control_points_suction_intersection(1,1) << "\n" <<
    //         control_points_pressure_intersection(0,0) << "\t\t\t" << control_points_pressure_intersection(0,1) << "\n" <<
    //         control_points_pressure_intersection(1,0) << "\t\t" << control_points_pressure_intersection(1,1) << "\n";

    //gsInfo << "Vysledne ridici body po uprave - suction:\n" << control_points_suction_intersection << "\n";
    //gsInfo << "Vysledne ridici body po uprave - pressure:\n" << control_points_pressure_intersection << "\n";

    //gsInfo << kvloft << "\n";
    //gsInfo << section_parameters << "\n";
    //gsInfo << suction_intersection_parameter_points_on_blade.row(0) << "\n";
    //gsInfo << pressure_intersection_parameter_points_on_blade.row(0) << "\n";

    int remove_profile = 0;
    bool change = false;
    for (index_t i = 0; i < section_parameters.cols(); i++) {
        change = false;
        for (index_t j = 0; j < suction_intersection_parameter_points_on_blade.cols(); j++) {
            if (suction_intersection_parameter_points_on_blade(0, j) > section_parameters(0, i)) {
                remove_profile = i;
                change = true;
            }
        }
        if (!change) {
            break;
        }
    }
    gsInfo << "First " << remove_profile+1 << " profiles will be removed because of trimming." << "\n";

    gsMatrix<T> section_parameters2(1, m_num_bladeprofiles-remove_profile);
    if (remove_profile > 0) {
        m_BladeProfiles[0]->setSuctionSide3d(suctionSideIntersection3d);
        m_BladeProfiles[0]->setPressureSide3d(pressureSideIntersection3d);
        section_parameters2(0) = 0;
        //gsInfo << "Pred forem\n";
        for (index_t i = 0; i < m_num_bladeprofiles-remove_profile-1; i++) {
            m_BladeProfiles[i+1] = m_BladeProfiles[remove_profile+i+1];
            //suction_side_curve3d[i+1] = suction_side_curve3d[remove_profile+i+1];
            //pressure_side_curve3d[i+1] = pressure_side_curve3d[remove_profile+i+1];
            section_parameters2(i) = (section_parameters(remove_profile+i)-section_parameters(remove_profile))/(section_parameters(m_num_bladeprofiles-1)-section_parameters(remove_profile));
        }
        m_num_bladeprofiles = m_num_bladeprofiles - remove_profile;
        section_parameters2(m_num_bladeprofiles-1) = 1.0;
        //gsInfo << m_num_bladeprofiles << "\n";
    }
    else {
        //suction_side_curve3d[0] = suctionSideIntersection3d;
        //pressure_side_curve3d[0] = pressureSideIntersection3d;
        m_BladeProfiles[0]->setSuctionSide3d(suctionSideIntersection3d);
        m_BladeProfiles[0]->setPressureSide3d(pressureSideIntersection3d);
        section_parameters2 = section_parameters;
    }
    mySectionParameters = section_parameters2;
    //gsInfo << section_parameters2 << "\n";
    gsInfo << "Done.\n";

    // Lofting 2 pro ziskani plochy lopatky po oriznuti vnitrni plochou turbiny
    gsInfo << "Construction of trimmed suction and pressure side surface of runner blade from 3D blade profiles ...\n";
    gsInfo << "Lofting parameters for 3D blade profiles: " << section_parameters2 << "\n";
    gsKnotVector<T> kvloft2(0, 1, m_num_bladeprofiles-4, 4);
    gsTensorBSpline<2,T> suction_side_surface_after_trimming;
    gsTensorBSpline<2,T> pressure_side_surface_after_trimming;
    //gsBSpline<T> defcurves_suctionside_after_trimming[m_num_bladeprofiles];
    //gsBSpline<T> defcurves_pressureside_after_trimming[m_num_bladeprofiles];
    std::vector<gsBSpline<T> > defcurves_suctionside_after_trimming;
    std::vector<gsBSpline<T> > defcurves_pressureside_after_trimming;
    //gsInfo << "a\n";
    for (index_t i = 0; i < m_num_bladeprofiles; i++) {
        //defcurves_suctionside_after_trimming[i] = m_BladeProfiles[i]->getSuctionSide3d();
        //defcurves_pressureside_after_trimming[i] = m_BladeProfiles[i]->getPressureSide3d();
        defcurves_suctionside_after_trimming.push_back(m_BladeProfiles[i]->getSuctionSide3d());
        defcurves_pressureside_after_trimming.push_back(m_BladeProfiles[i]->getPressureSide3d());
    }
    computeLoftSurface(defcurves_suctionside_after_trimming, kvfit, m_num_bladeprofiles, kvloft2, section_parameters2, suction_side_surface_after_trimming);
    computeLoftSurface(defcurves_pressureside_after_trimming, kvfit, m_num_bladeprofiles, kvloft2, section_parameters2, pressure_side_surface_after_trimming);
    mySuctionSideSurfaceAfterTrimming = suction_side_surface_after_trimming;
    myPressureSideSurfaceAfterTrimming = pressure_side_surface_after_trimming;
    gsInfo << "Done.\n\n";

    //gsBSpline<T> suction_side_curve3d[m_num_bladeprofiles];
    //gsBSpline<T> pressure_side_curve3d[m_num_bladeprofiles];
    std::vector<gsBSpline<T> > suction_side_curve3d;
    std::vector<gsBSpline<T> > pressure_side_curve3d;
    if (plot) {
        std::vector<gsGeometry<>*> curves3d;
        curves3d.clear();
        for (int i = 0; i < m_num_bladeprofiles; i++) {
            //suction_side_curve3d[i] = m_BladeProfiles[i]->getSuctionSide3d();
            suction_side_curve3d.push_back(m_BladeProfiles[i]->getSuctionSide3d());
            //pressure_side_curve3d[i] = m_BladeProfiles[i]->getPressureSide3d();
            pressure_side_curve3d.push_back(m_BladeProfiles[i]->getPressureSide3d());
        }
        for (int i = 0; i < m_num_bladeprofiles; i++) {
            curves3d.push_back(&suction_side_curve3d[i]);
            curves3d.push_back(&pressure_side_curve3d[i]);
        }
        gsWriteParaview( curves3d, "3DprofileAfterTrimming", 100);

        gsWriteParaview( suction_side_surface_after_trimming, "suctionSideSurfaceAfterTrimming", 5000);
        gsWriteParaview( pressure_side_surface_after_trimming, "pressureSideSurfaceAfterTrimming", 5000);

        gsMesh<> suction_side_mesh_after_trimming;
        gsMesh<> pressure_side_mesh_after_trimming;
        suction_side_surface_after_trimming.controlNet(suction_side_mesh_after_trimming);
        pressure_side_surface_after_trimming.controlNet(pressure_side_mesh_after_trimming);
        gsWriteParaview( suction_side_mesh_after_trimming, "suctionSideMeshAfterTrimming");
        gsWriteParaview( pressure_side_mesh_after_trimming, "pressureSideMeshAfterTrimming");
    }

    return 0;
}


template <class T>
int KaplanTurbineRunnerBlade<T>::exportCurveFile(std::string const & filename, bool plot, bool print_info) {

    int num_samples3d = 50;
    T par;
    int num_iter_s = 0;
    int max_num_iter_s = 0;
    int num_iter_p = 0;
    int max_num_iter_p = 0;
    T initial_solution;
    T solution;
    T radius = 0.499;

    gsInfo << "Trimming suction and pressure side surface by outer sphere\n";
    gsInfo << "==========================================================\n";

    // Computation of intersection of outer sphere with suction and pressure sides of the runner blade
    gsInfo << "Computing intersection points of suction side surface with outer sphere ...\n";
    gsMatrix<T> suction_intersection_parameter_points_on_blade(2, num_samples3d);
    gsMatrix<T> suction_intersection_points_on_blade(3, num_samples3d);
    for (index_t i = 0; i < num_samples3d; i++) {
        par = (T) i/(num_samples3d-1);
        initial_solution = 0.9;
        curveSphereIntersectionviaNR(mySuctionSideSurface, par, initial_solution, radius, solution, static_cast<real_t>(1e-5), 100, num_iter_s, true);
        suction_intersection_parameter_points_on_blade(0,i) = solution;
        suction_intersection_parameter_points_on_blade(1,i) = par;
        if (max_num_iter_s < num_iter_s) {
            max_num_iter_s = num_iter_s;
        }
    }
    if (print_info) {
        gsInfo << suction_intersection_parameter_points_on_blade << "\n";
        gsInfo << "Max number of NR iteration for suction side: " << max_num_iter_s << "\n\n";
    }
    mySuctionSideSurface.eval_into(suction_intersection_parameter_points_on_blade, suction_intersection_points_on_blade);
    gsInfo << "Done.\n";

    gsInfo << "Computing intersection points of pressure side surface with outer sphere ...\n";
    gsMatrix<T> pressure_intersection_parameter_points_on_blade(2, num_samples3d);
    gsMatrix<T> pressure_intersection_points_on_blade(3, num_samples3d);
    for (index_t i = 0; i < num_samples3d; i++) {
        par = (T) i/(num_samples3d-1);
        initial_solution = 0.9;
        curveSphereIntersectionviaNR(myPressureSideSurface, par, initial_solution, radius, solution, static_cast<real_t>(1e-5), 100, num_iter_p, false);
        pressure_intersection_parameter_points_on_blade(0,i) = solution;
        pressure_intersection_parameter_points_on_blade(1,i) = par;
        if (max_num_iter_p < num_iter_p) {
            max_num_iter_p = num_iter_p;
        }
    }
    if (print_info) {
        gsInfo << pressure_intersection_parameter_points_on_blade << "\n";
        gsInfo << "Max number of NR iteration for pressure side: " << max_num_iter_p << "\n\n";
    }
    myPressureSideSurface.eval_into(pressure_intersection_parameter_points_on_blade, pressure_intersection_points_on_blade);
    gsInfo << "Done.\n";

    if (plot) {
        gsWriteParaviewPoints<real_t>( suction_intersection_points_on_blade, "suction_intersection_points_with_sphere_on_blade");
        gsWriteParaviewPoints<real_t>( pressure_intersection_points_on_blade, "pressure_intersection_points_with_sphere_on_blade");
    }

    // Export blade profile curves of trimmed blade surface to file
    if (true) {
        gsMatrix<T> pars_for_export(1, num_samples3d);
        for (index_t i = 0; i < num_samples3d; i++) {
            pars_for_export(i) = (T) i/(num_samples3d-1);
        }

        std::ofstream exportFile;
        exportFile.open (filename);
        if (exportFile.is_open()) {
            exportFile << "##\n";
            for (index_t i = 0; i < m_num_bladeprofiles; i++) {
                exportFile << "#Rez c. " << i+1 << "\n";
                if (i == 0) {
                    // writing first profile obtained by trimming by inner surface
                    for (index_t j = num_samples3d-1; j >= 0; j--) {
                        exportFile << "\t" << std::setprecision(6) << std::setw(2) << mySuctionIntersectionPointonInnerHR.col(j).row(0) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << mySuctionIntersectionPointonInnerHR.col(j).row(1) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << mySuctionIntersectionPointonInnerHR.col(j).row(2) * 1000
                                   << "\n";
                    }
                    for (index_t j = 0; j < num_samples3d; j++) {
                        exportFile << "\t" << std::setprecision(6) << std::setw(2) << myPressureIntersectionPointonInnerHR.col(j).row(0) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << myPressureIntersectionPointonInnerHR.col(j).row(1) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << myPressureIntersectionPointonInnerHR.col(j).row(2) * 1000
                                   << "\n";
                    }
                }
                else if (i == m_num_bladeprofiles-1) {
                    // writing last profile obtained by trimming by outer sphere
                    for (index_t j = num_samples3d-1; j >= 0; j--) {
                        exportFile << "\t" << std::setprecision(6) << std::setw(2) << suction_intersection_points_on_blade.col(j).row(0) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << suction_intersection_points_on_blade.col(j).row(1) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << suction_intersection_points_on_blade.col(j).row(2) * 1000
                                   << "\n";
                    }
                    for (index_t j = 0; j < num_samples3d; j++) {
                        exportFile << "\t" << std::setprecision(6) << std::setw(2) << pressure_intersection_points_on_blade.col(j).row(0) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << pressure_intersection_points_on_blade.col(j).row(1) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << pressure_intersection_points_on_blade.col(j).row(2) * 1000
                                   << "\n";
                    }
                }
                else {
                    // writing other profiles by sampling profile curves
                    gsBSpline<T> curve_suction;
                    gsBSpline<T> curve_pressure;
                    curve_suction = m_BladeProfiles[i]->getSuctionSide3d();
                    curve_pressure = m_BladeProfiles[i]->getPressureSide3d();

                    gsMatrix<T> point_suction(3, num_samples3d);
                    gsMatrix<T> point_pressure(3, num_samples3d);
                    curve_suction.eval_into(pars_for_export, point_suction);
                    curve_pressure.eval_into(pars_for_export, point_pressure);

                    for (index_t j = num_samples3d-1; j >= 0; j--) {
                        exportFile << "\t" << std::setprecision(6) << std::setw(2) << point_suction.col(j).row(0) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << point_suction.col(j).row(1) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << point_suction.col(j).row(2) * 1000
                                   << "\n";
                    }
                    for (index_t j = 0; j < num_samples3d; j++) {
                        exportFile << "\t" << std::setprecision(6) << std::setw(2) << point_pressure.col(j).row(0) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << point_pressure.col(j).row(1) * 1000
                                   << "\t" << std::setprecision(6) << std::setw(2) << point_pressure.col(j).row(2) * 1000
                                   << "\n";
                    }
                }
            }
            exportFile.close();
        }
        else {
            gsWarn << "Unable to open export file." << '\n';
        }
    }
    //

    return 0;
}


#endif // UWBKAPLANTURBINERUNNERBLADE_H
