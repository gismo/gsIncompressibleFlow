#ifndef UWBBLADEPROFILE_H
#define UWBBLADEPROFILE_H

#include "uwbTurbineUtils.h"

real_t pi = 3.141529;

using namespace gismo;

template <class T>
class BladeProfile
{
public:
    //BladeProfile2D() : myCamberX(0.5), myCamberY(0.1), myLeadingAngle(0.4), myTrailingAngle(0.25), myThicknessX(0.5), myThicknessY(0.1), myEndingOffset(0), myOutputAngle(0.05), myRadius(-1) {}
    BladeProfile() {}
    BladeProfile(T thisCamberX, T thisCamberY, T thisLeadingAngle, T thisTrailingAngle, T thisThicknessX, T thisThicknessY, T thisEndingOffset, T thisOutputAngle,
                 T thisChordLength, T thisAngle, T thisRotationCenterX, T thisRotationCenterY, T thisCylinderRadius) :
        myCamberX(thisCamberX), myCamberY(thisCamberY), myLeadingAngle(thisLeadingAngle), myTrailingAngle(thisTrailingAngle),
        myThicknessX(thisThicknessX), myThicknessY(thisThicknessY), myEndingOffset(thisEndingOffset), myOutputAngle(thisOutputAngle), myRadius(-1),
        myChordLength(thisChordLength), myAngle(thisAngle), myRotationCenterX(thisRotationCenterX), myRotationCenterY(thisRotationCenterY),
        myCylinderRadius(thisCylinderRadius) {}
    BladeProfile(T thisCamberX, T thisCamberY, T thisLeadingAngle, T thisTrailingAngle, T thisThicknessX, T thisThicknessY, T thisEndingOffset, T thisOutputAngle,
                   T thisRadius, T thisChordLength, T thisAngle, T thisRotationCenterX, T thisRotationCenterY, T thisCylinderRadius) :
        myCamberX(thisCamberX), myCamberY(thisCamberY), myLeadingAngle(thisLeadingAngle), myTrailingAngle(thisTrailingAngle),
        myThicknessX(thisThicknessX), myThicknessY(thisThicknessY), myEndingOffset(thisEndingOffset), myOutputAngle(thisOutputAngle), myRadius(thisRadius),
        myChordLength(thisChordLength), myAngle(thisAngle), myRotationCenterX(thisRotationCenterX), myRotationCenterY(thisRotationCenterY),
        myCylinderRadius(thisCylinderRadius) {}
    BladeProfile(T thisCamberX, T thisCamberY, T thisLeadingAngle, T thisTrailingAngle, T thisThicknessX, T thisThicknessY, T thisEndingOffset, T thisOutputAngle,
                 T thisChordLength, T thisAngle, T thisRotationCenterX, T thisRotationCenterY, gsMatrix<T> thisCircularMeshSetting, gsMatrix<T> thisConeSetting) :
        myCamberX(thisCamberX), myCamberY(thisCamberY), myLeadingAngle(thisLeadingAngle), myTrailingAngle(thisTrailingAngle),
        myThicknessX(thisThicknessX), myThicknessY(thisThicknessY), myEndingOffset(thisEndingOffset), myOutputAngle(thisOutputAngle), myRadius(-1),
        myChordLength(thisChordLength), myAngle(thisAngle), myRotationCenterX(thisRotationCenterX), myRotationCenterY(thisRotationCenterY),
        myCircularMeshSetting(thisCircularMeshSetting), myConeSetting(thisConeSetting) {}
    BladeProfile(T thisCamberX, T thisCamberY, T thisLeadingAngle, T thisTrailingAngle, T thisThicknessX, T thisThicknessY, T thisEndingOffset, T thisOutputAngle,
                   T thisRadius, T thisChordLength, T thisAngle, T thisRotationCenterX, T thisRotationCenterY, gsMatrix<T> thisCircularMeshSetting, gsMatrix<T> thisConeSetting) :
        myCamberX(thisCamberX), myCamberY(thisCamberY), myLeadingAngle(thisLeadingAngle), myTrailingAngle(thisTrailingAngle),
        myThicknessX(thisThicknessX), myThicknessY(thisThicknessY), myEndingOffset(thisEndingOffset), myOutputAngle(thisOutputAngle), myRadius(thisRadius),
        myChordLength(thisChordLength), myAngle(thisAngle), myRotationCenterX(thisRotationCenterX), myRotationCenterY(thisRotationCenterY),
          myCircularMeshSetting(thisCircularMeshSetting), myConeSetting(thisConeSetting) {}
    ~BladeProfile() {}

    // Access function
    T getCamberX() const { return myCamberX; }
    T getCamberY() const { return myCamberY; }
    T getLeadingAngle() const { return myLeadingAngle; }
    T getTrailingAngle() const { return myTrailingAngle; }
    void setCamberX(T thisCamberX) { myCamberX = thisCamberX; }
    void setCamberY(T thisCamberY) { myCamberY = thisCamberY; }
    void setLeadingAngle(T thisLeadingAngle) { myLeadingAngle = thisLeadingAngle; }
    void setTrailingAngle(T thisTrailingAngle) { myTrailingAngle = thisTrailingAngle; }
    T getThicknessX() const { return myThicknessX; }
    T getThicknessY() const { return myThicknessY; }
    T getEndingOffset() const { return myEndingOffset; }
    T getOutputAngle() const { return myOutputAngle; }
    T getRadius() const { return myRadius; }
    void setThicknessX(T thisThicknessX) { myThicknessX = thisThicknessX; }
    void setThicknessY(T thisThicknessY) { myThicknessY = thisThicknessY; }
    void setEndingOffset(T thisEndingOffset) { myEndingOffset = thisEndingOffset; }
    void setOutputAngle(T thisOutputAngle) { myOutputAngle = thisOutputAngle; }
    void setRadius(T thisRadius) { myRadius = thisRadius; }
    gsKnotVector<T> getKnotVector() const { return myKnotVector; }
    void setKnotVector(gsKnotVector<T> kv) {myKnotVector = kv;}
    gsBSpline<T> getSuctionSide() const { return myBladeProfile2DSuctionSide; }
    gsBSpline<T> getPressureSide() const { return myBladeProfile2DPressureSide; }
    void setSuctionSide(gsBSpline<T> bspl) { myBladeProfile2DSuctionSide = bspl; }
    void setPressureSide(gsBSpline<T> bspl) { myBladeProfile2DPressureSide = bspl; }
    T getChordLength() const { return myChordLength; }
    T getAngle() const { return myAngle; }
    T getRotationCenterX() const { return myRotationCenterX; }
    T getRotationCenterY() const { return myRotationCenterY; }
    void setCylinderRadius(T rr) {myCylinderRadius = rr;}
    T getCylinderRadius() const { return myCylinderRadius; }
    void setCircularMeshSetting(gsMatrix<T> mat) {myCircularMeshSetting = mat;}
    gsMatrix<T> getCircularMeshSetting() const { return myCircularMeshSetting; }
    void setConSetting(gsMatrix<T> mat) {myConeSetting = mat;}
    gsMatrix<T> getConeSetting() const { return myConeSetting; }
    gsBSpline<T> getSuctionSide3d() const { return myBladeProfile3DSuctionSide; }
    gsBSpline<T> getPressureSide3d() const { return myBladeProfile3DPressureSide; }
    void setSuctionSide3d(gsBSpline<T> bspl) { myBladeProfile3DSuctionSide = bspl; }
    void setPressureSide3d(gsBSpline<T> bspl) { myBladeProfile3DPressureSide = bspl; }
    gsBSpline<T> getSuctionSideOffset() const { return myBladeProfile2DSuctionSideOffset; }
    gsBSpline<T> getPressureSideOffset() const { return myBladeProfile2DPressureSideOffset; }
    void setSuctionSideOffset(gsBSpline<T> bspl) { myBladeProfile2DSuctionSideOffset = bspl; }
    void setPressureSideOffset(gsBSpline<T> bspl) { myBladeProfile2DPressureSideOffset = bspl; }

    // Other functions
    int computeCamberLine(gsBSpline<T> & curveclc, int deg, int cnt);
    int computeThicknessFunction(gsBSpline<T> & ThicknessFunction, int deg, int cnt);
    int compute2D(gsBSpline<T> & suctionSide, gsBSpline<T> & pressureSide, int num_samples, int deg, int cnt);
    int compute2D(gsBSpline<T> & suctionSide, gsBSpline<T> & pressureSide, gsKnotVector<T> kvfit, int num_samples, int deg, int cnt);
    int compute3D(gsBSpline<T> & suctionSide3d, gsBSpline<T> & pressureSide3d, gsKnotVector<T> kvfit, int num_samples3d);
    int compute3DGuideVane(gsBSpline<T> & suctionSide3d, gsBSpline<T> & pressureSide3d, gsKnotVector<T> kvfit, int num_samples3d);
    int computeOffset(T const offsetdist, gsBSpline<T> & suctionSideOffset, gsBSpline<T> & pressureSideOffset, int num_samples);
    int computeOffset(T const offsetdist, gsBSpline<T> & suctionSideOffset, gsBSpline<T> & pressureSideOffset, gsKnotVector<T> kvfit, int num_samples);

    // additional math functions
    T csc(T x) { return 1/sin(x);}
    T cot(T x) { return 1/tan(x);}


private:
    T myCamberX;
    T myCamberY;
    T myLeadingAngle;
    T myTrailingAngle;
    gsBSpline<T> myCamberLine;
    T myThicknessX;
    T myThicknessY;
    T myEndingOffset;
    T myOutputAngle;
    T myRadius;
    gsBSpline<T> myThicknessFunction;
    gsKnotVector<T> myKnotVector;
    gsBSpline<T> myBladeProfile2DSuctionSide;
    gsBSpline<T> myBladeProfile2DPressureSide;
    gsBSpline<T> myBladeProfile2DSuctionSideOffset;
    gsBSpline<T> myBladeProfile2DPressureSideOffset;
    T myChordLength;
    T myAngle;
    T myRotationCenterX;
    T myRotationCenterY;
    T myCylinderRadius;
    gsMatrix<T> myCircularMeshSetting;
    gsMatrix<T> myConeSetting;
    gsBSpline<T> myBladeProfile3DSuctionSide;
    gsBSpline<T> myBladeProfile3DPressureSide;
};

template <class T>
int BladeProfile<T>::computeCamberLine(gsBSpline<T> & CLine, int deg, int cnt)
{
    if (atan(myCamberY/myCamberX) > myLeadingAngle) {
        gsWarn << "Leading edge angle too small with respect to maximal camber and its position. Possible problems with shape of camber line!\n";
    }

    if (atan(myCamberY/(1-myCamberX)) > myTrailingAngle) {
        gsWarn << "Trailing edge angle too small with respect to maximal camber and its position. Possible problems with shape of camber line!\n";
    }
    switch (deg)
    {
        case 2:
        {
            if (myLeadingAngle == 0) {
                GISMO_ERROR("Zero leading edge angle entered - G1 quadratic camber line cannot be constructed.\n");
                return 1;
            }
            if (myTrailingAngle == 0) {
                GISMO_ERROR("Zero trailing edge angle entered - G1 quadratic camber line cannot be constructed.\n");
                return 1;
            }
            gsKnotVector<T> kvclp(0, 1, 1, 3);//start,end,interior knots, start/end multiplicites of knots1
            kvclp.insert(0.5, 1);
            gsMatrix<T> coefsclp(5, 2);
            coefsclp <<    0, 0,
                        myCamberY * 1/tan(myLeadingAngle), myCamberY,
                        myCamberX, myCamberY,
                        1 - (myCamberY * 1/tan(myTrailingAngle)), myCamberY,
                        1, 0;

            CLine = gsBSpline<T>( kvclp, give(coefsclp));
            myCamberLine = CLine;

            return 0;

        }
        break;

        case 3:
        {

            switch (cnt)
            {
            case 1:
            {
                // C1 continuous cubic camber line

                T P2x;
                T a;
                T b;

                gsKnotVector<T> kvclp(0, 1, 1, 4);
                kvclp.insert(0.5, 2);

                P2x = (3 - 14*myCamberX + (-1 + 2*myCamberX)*cos(2*myTrailingAngle) - myCamberY*(sin(2*myLeadingAngle) + sin(2*myTrailingAngle)))/
                        (-14 + cos(2*myLeadingAngle) + cos(2*myTrailingAngle));
                a = (myCamberY*(-15 + cos(2*myTrailingAngle))*sin(myLeadingAngle) + cos(myLeadingAngle)*(3 - 14*myCamberX + (-1 + 2*myCamberX)*cos(2*myTrailingAngle) - myCamberY*sin(2*myTrailingAngle)))/
                        (2.*(-14 + cos(2*myLeadingAngle) + cos(2*myTrailingAngle)));
                b = -(cos(myTrailingAngle)*(myCamberY*cos(myLeadingAngle)*sin(myLeadingAngle) + 2*pow(sin(myLeadingAngle),2)*(3 - 4*myCamberX + 4*myCamberY*tan(myTrailingAngle)) +
                        pow(cos(myLeadingAngle),2)*(5 - 6*myCamberX + 7*myCamberY*tan(myTrailingAngle))))/(-14 + cos(2*myLeadingAngle) + cos(2*myTrailingAngle));

                gsInfo << " P2x = " << P2x << "\n a = " << a << "\n b = " << b << "\n";

                if (a < 0) {
                    GISMO_ERROR("Leading edge tangent vector multiplier negative. Possible problems with shape of camber line!\n");
                    return 1;
                }

                if (b < 0) {
                    GISMO_ERROR("Trailing edge tangent vector multiplier negative. Possible problems with shape of camber line!\n");
                    return 1;
                }

                gsMatrix<T> coefsclp(7, 2);
                coefsclp << 0, 0,
                            a*cos(myLeadingAngle), a*sin(myLeadingAngle),
                            P2x, myCamberY,
                            myCamberX, myCamberY,
                            2*myCamberX-P2x, myCamberY,
                            1 - b*cos(myTrailingAngle), b*sin(myTrailingAngle),
                            1, 0;

                CLine = gsBSpline<T>( kvclp, give(coefsclp));
                myCamberLine = CLine;
            }
            break;

            case 2:
            {
                // C2 continuous cubic camber line

                /*T P2x;
                T a;

                gsKnotVector<T> kvclp(0, 1, 1, 4);
                kvclp.insert(0.5, 2);

                P2x = ((-1 + 6*myCamberX)*(-3 + cos(2*myLeadingAngle)) - 2*myCamberY*sin(2*myLeadingAngle))/(-20 + 8*cos(2*myLeadingAngle));
                a = (1/cos(myLeadingAngle)*(-1 + 6*myCamberX + 14*myCamberY*tan(myLeadingAngle)))/(4.*(3 + 7*pow(tan(myLeadingAngle),2)));

                if (a < 0) {
                    GISMO_ERROR("Leading edge tangent vector multiplier negative. Possible problems with shape of camber line!");
                    return 1;

                gsMatrix<T> coefsclp(7, 2);
                coefsclp << 0, 0,
                            a*cos(myLeadingAngle), a*sin(myLeadingAngle),
                            P2x, myCamberY,
                            myCamberX, myCamberY,
                            2*myCamberX-P2x, myCamberY,
                            4*myCamberX - 4*P2x + a*cos(myLeadingAngle), a*sin(myLeadingAngle),
                            1, 0;
                }*/

                T P2x;
                T a;
                T b;

                gsKnotVector<T> kvclp(0, 1, 1, 4);
                kvclp.insert(0.5, 2);

                P2x = (33 - 138*myCamberX + 2*(-7 + 34*myCamberX)*cos(2*myLeadingAngle) + (-26 + 92*myCamberX)*cos(2*myTrailingAngle) + (7 - 22*myCamberX)*cos(2*(myLeadingAngle + myTrailingAngle)) -
                       48*myCamberY*sin(myLeadingAngle)*sin(myTrailingAngle)*sin(myLeadingAngle + myTrailingAngle))/(-138 + 80*cos(2*myLeadingAngle) + 80*cos(2*myTrailingAngle) - 22*cos(2*(myLeadingAngle + myTrailingAngle)));
                a = (1 - 4*myCamberX + 4*P2x)*csc(myLeadingAngle + myTrailingAngle)*sin(myTrailingAngle);
                b = (1 - 4*myCamberX + 4*P2x)*csc(myLeadingAngle + myTrailingAngle)*sin(myLeadingAngle);

                gsMatrix<T> coefsclp(7, 2);
                coefsclp << 0, 0,
                            a*cos(myLeadingAngle), a*sin(myLeadingAngle),
                            P2x, myCamberY,
                            myCamberX, myCamberY,
                            2*myCamberX-P2x, myCamberY,
                            1 - b*cos(myTrailingAngle), b*sin(myTrailingAngle),
                            1, 0;

                CLine = gsBSpline<T>( kvclp, give(coefsclp));
                myCamberLine = CLine;
            }
            break;
            }

            return 0;
        }
        break;

        default:
        {
            GISMO_ERROR("Unsupported degree of camber line entered!\n");
            gsKnotVector<T> kvclp(0, 1, 0, 1);
            gsMatrix<T> coefsclp(2, 2);
            coefsclp << 0, 0,
                        0, 0;

            CLine = gsBSpline<T>( kvclp, give(coefsclp));
            myCamberLine = CLine;

            return 1;
        }
        break;

    }
}

template <class T>
int BladeProfile<T>::computeThicknessFunction(gsBSpline<T> & THFunction, int deg, int cnt)
{
    if (atan((myThicknessY-myEndingOffset)/(1-myThicknessX)) > myOutputAngle) {
        gsWarn << "Output angle of thickness function too small with respect to maximal thickness and its position. Possible problems with shape of thickness function!\n";
    }

    switch (deg)
    {
        case 2:
        {
            gsKnotVector<T> kvclp(0, 1, 1, 3);//start,end,interior knots, start/end multiplicites of knots1
            kvclp.insert(0.5, 1);
            gsMatrix<T> coefsclp(5, 2);
            coefsclp <<    0, 0,
                        0, myThicknessY,
                        myThicknessX, myThicknessY,
                        1 - myThicknessY/tan(myOutputAngle), myThicknessY,
                        1, myEndingOffset;


            THFunction = gsBSpline<T>( kvclp, give(coefsclp));
            myThicknessFunction = THFunction;

            return 0;

        }
        break;

        case 3:
        {
            switch (cnt)
            {
            case 1:
            {
                // C1 continuous cubic camber line
                // Derived for special fixed choice of P2x

                T P2x;
                T b;
                T P4x;

                gsKnotVector<T> kvclp(0, 1, 1, 4);
                kvclp.insert(0.5, 2);

                if (myThicknessX < 0.5) {
                    P2x = myThicknessX/2;
                    b = ((2 - 3*myThicknessX)*cos(myOutputAngle) + 2*(-myEndingOffset + myThicknessY)*sin(myOutputAngle))/4;
                    P4x = 3*myThicknessX/2;
                }
                else {
                    P2x = (-1 + 3*myThicknessX)/2;
                    b = (cos(myOutputAngle) - myThicknessX*cos(myOutputAngle) - 2*myEndingOffset*sin(myOutputAngle) + 2*myThicknessY*sin(myOutputAngle))/4;
                    P4x = (1 + myThicknessX)/2;
                }

                gsInfo << "P1y = " << sqrt(myRadius)*sqrt(myThicknessX)/sqrt(2) << ", b = " << b << ", P4x = " << P4x << ", P2x = " << P2x << "\n";

                if (b < 0) {
                    GISMO_ERROR("Trailing edge tangent vector multiplier negative. Possible problems with shape of camber line!\n");
                    return 1;
                }

                gsMatrix<T> coefsclp(7, 2);
                coefsclp << 0, 0,
                            0, sqrt(myRadius)*sqrt(myThicknessX)/sqrt(2),
                            P2x, myThicknessY,
                            myThicknessX, myThicknessY,
                            P4x, myThicknessY,
                            1 - b*cos(myOutputAngle), myEndingOffset + b*sin(myOutputAngle),
                            1, myEndingOffset;

                THFunction = gsBSpline<T>( kvclp, give(coefsclp));
                myThicknessFunction = THFunction;

            }
            break;

            case 2:
            {
                // C2 continuous cubic camber line

                /*
                gsKnotVector<T> kvclp(0, 1, 1, 4);
                kvclp.insert(0.5, 2);

                gsMatrix<T> coefsclp(7, 2);
                if (myThicknessX > 0.5) {
                    coefsclp << 0, 0,
                                0, sqrt(myRadius)*sqrt(myThicknessX+(myThicknessX-1)/2),
                                (3*myThicknessX-1)/2, myThicknessY,
                                myThicknessX, myThicknessY,
                                (1+myThicknessX)/2, myThicknessY,
                                2-2*myThicknessX, sqrt(myRadius)*sqrt(myThicknessX+(myThicknessX-1)/2),
                                1, myEndingOffset;
                }
                else {
                    coefsclp << 0, 0,
                                0, sqrt(myRadius)*sqrt(myThicknessX)/sqrt(2),
                                myThicknessX/2, myThicknessY,
                                myThicknessX, myThicknessY,
                                3*myThicknessX/2, myThicknessY,
                                2*myThicknessX, sqrt(myRadius)*sqrt(myThicknessX)/sqrt(2),
                                1, myEndingOffset;
                }
                */

                gsKnotVector<T> kvclp(0, 1, 1, 4);
                kvclp.insert(0.5, 2);

                T P1ya, P1yb, P1y;
                T P2x;
                T b;
                T d1;
                T d2;

                //P1y = sqrt(2/3)*sqrt(myRadius)*sqrt(P2x);
                d1 = myRadius*(-24 + 96*myThicknessX - 24*myEndingOffset*cot(myOutputAngle) + myRadius*pow(cot(myOutputAngle),2))*pow(sin(myOutputAngle),2);
                d2 = myRadius*(24 - 96*myThicknessX + 24*myEndingOffset*cot(myOutputAngle) + myRadius*pow(cot(myOutputAngle),2))*pow(sin(myOutputAngle),2);
                if (d1 >=0) {
                    P1ya = ((sqrt(d1) + myRadius*cos(myOutputAngle))*csc(myOutputAngle))/12.;
                    P1yb = -((sqrt(d1) - myRadius*cos(myOutputAngle))*csc(myOutputAngle))/12.;
                    if (P1ya > 0)
                        P1y = P1ya;
                    else
                        P1y = P1yb;
                }
                else if (d2 >= 0) {
                    P1ya = (-(myRadius*cot(myOutputAngle)) + sqrt(d2)*csc(myOutputAngle))/12.;
                    P1yb = -((sqrt(d2) + myRadius*cos(myOutputAngle))*csc(myOutputAngle))/12.;
                    if (P1ya > 0)
                        P1y = P1ya;
                    else
                        P1y = P1yb;
                }
                else {
                    GISMO_ERROR("Something went wrong during computation of thickness function. Exiting ...\n");
                    return 1;
                }
                P2x = -0.25 + myThicknessX + ((-myEndingOffset + P1y)*cot(myOutputAngle))/4.;
                b = (-myEndingOffset + P1y)*csc(myOutputAngle);

                gsMatrix<T> coefsclp(7, 2);
                coefsclp << 0, 0,
                            0, P1y,
                            P2x, myThicknessY,
                            myThicknessX, myThicknessY,
                            2*myThicknessX-P2x, myThicknessY,
                            1-b*cos(myOutputAngle), myEndingOffset+b*sin(myOutputAngle),
                            1, myEndingOffset;

                THFunction = gsBSpline<T>( kvclp, give(coefsclp));
                myThicknessFunction = THFunction;
            }
            break;
            }

            return 0;
        }
        break;

        default:
        {
            gsInfo << "Unsupported degree of thickness function entered!\n";
            gsKnotVector<T> kvclp(0, 1, 0, 1);
            gsMatrix<T> coefsclp(2, 2);
            coefsclp << 0, 0,
                        0, 0;

            THFunction = gsBSpline<T>( kvclp, give(coefsclp));
            myThicknessFunction = THFunction;
            return 1;
        }
        break;

    }
}

template <class T>
int BladeProfile<T>::compute2D(gsBSpline<T> & suctionSide, gsBSpline<T> & pressureSide, int num_samples, int deg, int cnt)
{
    int num_samples_fortf = 1000;

    gsKnotVector<T> kvfit = myKnotVector;
    gsBSpline<T> CLcurve;
    gsBSpline<T> THcurve;
    gsMatrix<T> parameter_points(1, num_samples+1);
    gsMatrix<T> parameter_points_forth(1, num_samples_fortf+1);
    gsMatrix<T> parameter_points_suction(1, num_samples+1);
    gsMatrix<T> parameter_points_pressure(1, num_samples+1);
    //gsMatrix<T> parameter_points_thickness(1, num_samples+1);
    gsMatrix<T> control_points_suction(kvfit.size()-kvfit.degree()-1, 2);
    gsMatrix<T> control_points_pressure(kvfit.size()-kvfit.degree()-1, 2);
    gsMatrix<T> control_points_pom(2*(kvfit.size()-kvfit.degree()-1), 2);
    gsMatrix<T> camber_line_points(2, num_samples+1);
    gsMatrix<T> camber_line_tangents(2, num_samples+1);
    gsMatrix<T> camber_line_normals(2, num_samples+1);
    gsMatrix<T> thickness_function_points_initial(2, num_samples_fortf+1);
    gsMatrix<T> thickness_function_points(2, num_samples+1);
    gsMatrix<T> suction_side_points(2, num_samples+1);
    gsMatrix<T> pressure_side_points(2, num_samples+1);
    gsMatrix<T> camber_line_normal_leading(2,1);
    gsMatrix<T> camber_line_normal_leading2(2,1);
    T alpha;
    int i, j;

    computeCamberLine(CLcurve, deg, cnt);
    //gsInfo << CLcurve << "\n";
    computeThicknessFunction(THcurve, deg, cnt);
    //gsInfo << THcurve << "\n";

    //gsInfo << "Blade: thicknessx = " << getThicknessX() << ", thicknessy = " << getThicknessY() << ", output angle = " << getOutputAngle() << ", ending offset = " << getEndingOffset()
    //       << ", radius = " << getRadius() << "\n";

    // Sampling point in the parameter domain of camber line
    for (i = 0; i < num_samples; i++) {
        parameter_points(0,i) = (1-cos(pi/2*i/num_samples));
    }
    parameter_points(0,num_samples) = 1;

    camber_line_points = CLcurve.eval(parameter_points);

    // Finding corresponding point on thickness function for camber line points
    for (i = 0; i < num_samples_fortf; i++) {
        parameter_points_forth(0,i) = (1-cos(pi/2*i/num_samples_fortf));
    }
    parameter_points_forth(0,num_samples_fortf) = 1;
    thickness_function_points_initial = THcurve.eval(parameter_points);
    for (i=0; i<=num_samples; i++) {
        for (j=0; j<num_samples_fortf; j++) {
            if ((camber_line_points(0,i) >= thickness_function_points_initial(0,j)) && (camber_line_points(0,i) <= thickness_function_points_initial(0,j+1))) {
                break;
            }
        }
        alpha = (camber_line_points(0,i) - thickness_function_points_initial(0,j))/(thickness_function_points_initial(0,j+1) - thickness_function_points_initial(0,j));
        //parameter_points_thickness(0, i) = (1 - alpha)*parameter_points(0, j) + alpha*parameter_points(0, j+1);
        thickness_function_points(0, i) = (1 - alpha)*thickness_function_points_initial(0, j) + alpha*thickness_function_points_initial(0, j+1);
        thickness_function_points(1, i) = (1 - alpha)*thickness_function_points_initial(1, j) + alpha*thickness_function_points_initial(1, j+1);
    }
    //thickness_function_points = THcurve.eval(parameter_points_thickness);

    // Applying thickness function to the camber line to obtain points on the pressure and suction side of the blade profile
    camber_line_tangents = CLcurve.deriv(parameter_points);
    for (i = 0; i <= num_samples; i++) {
        camber_line_normals(0, i) = - camber_line_tangents(1, i)/camber_line_tangents.col(i).norm();
        camber_line_normals(1, i) = camber_line_tangents(0, i)/camber_line_tangents.col(i).norm();
        suction_side_points(0, i) = camber_line_points(0, i) + thickness_function_points(1, i)*camber_line_normals(0, i);
        suction_side_points(1, i) = camber_line_points(1, i) + thickness_function_points(1, i)*camber_line_normals(1, i);
        pressure_side_points(0, i) = camber_line_points(0, i) - thickness_function_points(1, i)*camber_line_normals(0, i);
        pressure_side_points(1, i) = camber_line_points(1, i) - thickness_function_points(1, i)*camber_line_normals(1, i);
    }
    camber_line_normal_leading.col(0) = camber_line_normals.col(0);
    camber_line_normal_leading2.col(0) = -camber_line_normals.col(0);

    // Fitting B-spline on suction side
    suction_side_points = suction_side_points.transpose();
    parameter_points_suction = centripetalParameterization(suction_side_points);
    suctionSide = curveFittingWithBoundaryAndInputTangent(suction_side_points, parameter_points_suction, kvfit, camber_line_normal_leading);
    control_points_suction = suctionSide.coefs();

    // Fitting B-spline on pressure side
    pressure_side_points = pressure_side_points.transpose();
    parameter_points_pressure = centripetalParameterization(pressure_side_points);
    //gsBSpline<T> pressure_side_curve = curveFittingWithBoundary(pressure_side_points, parameter_points, kvfit);
    pressureSide = curveFittingWithBoundaryAndInputTangent(pressure_side_points, parameter_points_pressure, kvfit, camber_line_normal_leading2);
    control_points_pressure = pressureSide.coefs();

    // Making a joint curve of suction and pressure sides C1 continuous (taking smaller tangent vector and correcting the other curve)x
    T l1 = sqrt(pow(control_points_suction(1,0)-control_points_suction(0,0),2) + pow(control_points_suction(1,1)-control_points_suction(0,1),2));
    T l2 = sqrt(pow(control_points_pressure(1,0)-control_points_pressure(0,0),2) + pow(control_points_pressure(1,1)-control_points_pressure(0,1),2));
    //gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction(0,0) << "\t\t\t" << control_points_suction(0,1) << "\n" <<
    //          control_points_suction(1,0) << "\t\t" << control_points_suction(1,1) << "\n" <<
    //          control_points_pressure(0,0) << "\t\t\t" << control_points_pressure(0,1) << "\n" <<
    //          control_points_pressure(1,0) << "\t\t" << control_points_pressure(1,1) << "\n";
    if (l1 > l2) {
        control_points_suction(1,0) = control_points_suction(0,0) + l2*(control_points_suction(1,0) - control_points_suction(0,0))/l1;
        control_points_suction(1,1) = control_points_suction(0,1) + l2*(control_points_suction(1,1) - control_points_suction(0,1))/l1;
    //    gsInfo << "Upravene body:\n" << control_points_suction << "\n";
        suctionSide = gsBSpline<T>( kvfit, control_points_suction);
    }
    else {
        control_points_pressure(1,0) = control_points_pressure(0,0) + l1/l2*(control_points_pressure(1,0) - control_points_pressure(0,0));
        control_points_pressure(1,1) = control_points_pressure(0,1) + l1/l2*(control_points_pressure(1,1) - control_points_pressure(0,1));
    //    gsInfo << "Upravene body:\n" << control_points_pressure << "\n";
        pressureSide = gsBSpline<T>( kvfit, control_points_pressure);
    }

    l1 = sqrt(pow(control_points_suction(1,0)-control_points_suction(0,0),2) + pow(control_points_suction(1,1)-control_points_suction(0,1),2));
    l2 = sqrt(pow(control_points_pressure(1,0)-control_points_pressure(0,0),2) + pow(control_points_pressure(1,1)-control_points_pressure(0,1),2));
    //gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction(0,0) << "\t\t\t" << control_points_suction(0,1) << "\n" <<
    //         control_points_suction(1,0) << "\t\t" << control_points_suction(1,1) << "\n" <<
    //         control_points_pressure(0,0) << "\t\t\t" << control_points_pressure(0,1) << "\n" <<
    //         control_points_pressure(1,0) << "\t\t" << control_points_pressure(1,1) << "\n";

    myBladeProfile2DPressureSide = pressureSide;
    myBladeProfile2DSuctionSide = suctionSide;


    // Output
    // Point-wise output of pressure and suction side to Paraview
    gsWriteParaviewPoints<T>( suction_side_points.transpose(), "suction_side_points");
    gsWriteParaviewPoints<T>( pressure_side_points.transpose(), "pressure_side_points");

    return 0;
}

template <class T>
int BladeProfile<T>::compute2D(gsBSpline<T> & suctionSide, gsBSpline<T> & pressureSide, gsKnotVector<T> kvfit, int num_samples, int deg, int cnt)
{
    setKnotVector(kvfit);

    compute2D(suctionSide, pressureSide, num_samples, deg, cnt);

    return 0;
}

template <class T>
int BladeProfile<T>::computeOffset(T const offsetdist, gsBSpline<T> & suctionSideOffset, gsBSpline<T> & pressureSideOffset, int num_samples)
{
    //int num_samples_fortf = 1000;

    gsKnotVector<T> kvfit = myKnotVector;
    gsBSpline<T> pressureSide;
    gsBSpline<T> suctionSide;
    gsMatrix<T> parameter_points(1, num_samples+1);
    gsMatrix<T> parameter_points_suction_offset(1, num_samples+1);
    gsMatrix<T> parameter_points_pressure_offset(1, num_samples+1);
    gsMatrix<T> suction_side_points(2, num_samples+1);
    gsMatrix<T> pressure_side_points(2, num_samples+1);
    gsMatrix<T> suction_side_tangents(2, num_samples+1);
    gsMatrix<T> pressure_side_tangents(2, num_samples+1);
    gsMatrix<T> suction_side_normals(2, num_samples+1);
    gsMatrix<T> pressure_side_normals(2, num_samples+1);
    gsMatrix<T> suction_side_offset_points(2, num_samples+1);
    gsMatrix<T> pressure_side_offset_points(2, num_samples+1);

    gsMatrix<T> control_points_suction_offset(kvfit.size()-kvfit.degree()-1, 2);
    gsMatrix<T> control_points_pressure_offset(kvfit.size()-kvfit.degree()-1, 2);
    gsMatrix<T> control_points_pom(2*(kvfit.size()-kvfit.degree()-1), 2);

    gsMatrix<T> pressure_side_tangent_leading(2,1);
    gsMatrix<T> suction_side_tangent_leading(2,1);
    int i;

    // Sampling points in the parameter domain
    for (i = 0; i < num_samples; i++) {
        parameter_points(0,i) = (1-cos(pi/2*i/num_samples));
    }
    parameter_points(0,num_samples) = 1;

    pressureSide = myBladeProfile2DPressureSide;
    suctionSide = myBladeProfile2DSuctionSide;
    pressure_side_points = pressureSide.eval(parameter_points);
    suction_side_points = suctionSide.eval(parameter_points);
    pressure_side_tangents = pressureSide.deriv(parameter_points);
    suction_side_tangents = suctionSide.deriv(parameter_points);

    // Applying thickness function to the camber line to obtain points on the pressure and suction side of the blade profile
    for (i = 0; i <= num_samples; i++) {
        pressure_side_normals(0, i) = - pressure_side_tangents(1, i)/pressure_side_tangents.col(i).norm();
        pressure_side_normals(1, i) = pressure_side_tangents(0, i)/pressure_side_tangents.col(i).norm();
        pressure_side_offset_points(0, i) = pressure_side_points(0, i) + offsetdist * pressure_side_normals(0, i);
        pressure_side_offset_points(1, i) = pressure_side_points(1, i) + offsetdist * pressure_side_normals(1, i);
        suction_side_normals(0, i) = - suction_side_tangents(1, i)/suction_side_tangents.col(i).norm();
        suction_side_normals(1, i) = suction_side_tangents(0, i)/suction_side_tangents.col(i).norm();
        suction_side_offset_points(0, i) = suction_side_points(0, i) - offsetdist * suction_side_normals(0, i);
        suction_side_offset_points(1, i) = suction_side_points(1, i) - offsetdist * suction_side_normals(1, i);
    }
    pressure_side_tangent_leading.col(0) = pressure_side_tangents.col(0);
    suction_side_tangent_leading.col(0) = suction_side_tangents.col(0);

    gsInfo << "suction points:\n" << pressure_side_points << "\n";
    gsInfo << "suction offset points:\n" << pressure_side_offset_points << "\n";


    // Fitting B-spline on suction side
    suction_side_points = suction_side_points.transpose();
    suction_side_offset_points = suction_side_offset_points.transpose();
    parameter_points_suction_offset = centripetalParameterization(suction_side_points);
    suctionSideOffset = curveFittingWithBoundaryAndInputTangent(suction_side_offset_points, parameter_points_suction_offset, kvfit, suction_side_tangent_leading);
    control_points_suction_offset = suctionSideOffset.coefs();

    // Fitting B-spline on pressure side
    pressure_side_points = pressure_side_points.transpose();
    pressure_side_offset_points = pressure_side_offset_points.transpose();
    parameter_points_pressure_offset = centripetalParameterization(pressure_side_points);
    //gsBSpline<T> pressure_side_curve = curveFittingWithBoundary(pressure_side_points, parameter_points, kvfit);
    pressureSideOffset = curveFittingWithBoundaryAndInputTangent(pressure_side_offset_points, parameter_points_pressure_offset, kvfit, pressure_side_tangent_leading);
    control_points_pressure_offset = pressureSideOffset.coefs();

    gsInfo << suctionSideOffset << "\n";
    gsInfo << pressureSideOffset << "\n";

    // Making a joint curve of suction and pressure sides C1 continuous (taking smaller tangent vector and correcting the other curve)x
    T l1 = sqrt(pow(control_points_suction_offset(1,0)-control_points_suction_offset(0,0),2) + pow(control_points_suction_offset(1,1)-control_points_suction_offset(0,1),2));
    T l2 = sqrt(pow(control_points_pressure_offset(1,0)-control_points_pressure_offset(0,0),2) + pow(control_points_pressure_offset(1,1)-control_points_pressure_offset(0,1),2));
    //gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction(0,0) << "\t\t\t" << control_points_suction(0,1) << "\n" <<
    //          control_points_suction(1,0) << "\t\t" << control_points_suction(1,1) << "\n" <<
    //          control_points_pressure(0,0) << "\t\t\t" << control_points_pressure(0,1) << "\n" <<
    //          control_points_pressure(1,0) << "\t\t" << control_points_pressure(1,1) << "\n";
    if (l1 > l2) {
        control_points_suction_offset(1,0) = control_points_suction_offset(0,0) + l2*(control_points_suction_offset(1,0) - control_points_suction_offset(0,0))/l1;
        control_points_suction_offset(1,1) = control_points_suction_offset(0,1) + l2*(control_points_suction_offset(1,1) - control_points_suction_offset(0,1))/l1;
    //    gsInfo << "Upravene body:\n" << control_points_suction << "\n";
        suctionSideOffset = gsBSpline<T>( kvfit, control_points_suction_offset);
    }
    else {
        control_points_pressure_offset(1,0) = control_points_pressure_offset(0,0) + l1/l2*(control_points_pressure_offset(1,0) - control_points_pressure_offset(0,0));
        control_points_pressure_offset(1,1) = control_points_pressure_offset(0,1) + l1/l2*(control_points_pressure_offset(1,1) - control_points_pressure_offset(0,1));
    //    gsInfo << "Upravene body:\n" << control_points_pressure << "\n";
        pressureSideOffset = gsBSpline<T>( kvfit, control_points_pressure_offset);
    }

    l1 = sqrt(pow(control_points_suction_offset(1,0)-control_points_suction_offset(0,0),2) + pow(control_points_suction_offset(1,1)-control_points_suction_offset(0,1),2));
    l2 = sqrt(pow(control_points_pressure_offset(1,0)-control_points_pressure_offset(0,0),2) + pow(control_points_pressure_offset(1,1)-control_points_pressure_offset(0,1),2));
    //gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction(0,0) << "\t\t\t" << control_points_suction(0,1) << "\n" <<
    //         control_points_suction(1,0) << "\t\t" << control_points_suction(1,1) << "\n" <<
    //         control_points_pressure(0,0) << "\t\t\t" << control_points_pressure(0,1) << "\n" <<
    //         control_points_pressure(1,0) << "\t\t" << control_points_pressure(1,1) << "\n";

    myBladeProfile2DPressureSideOffset = pressureSideOffset;
    myBladeProfile2DSuctionSideOffset = suctionSideOffset;


    // Output
    // Point-wise output of pressure and suction side to Paraview
    gsWriteParaviewPoints<real_t>( suction_side_points.transpose(), "suction_side_points");
    gsWriteParaviewPoints<real_t>( pressure_side_points.transpose(), "pressure_side_points");
    gsWriteParaviewPoints<real_t>( suction_side_offset_points.transpose(), "suction_side_offset_points");
    gsWriteParaviewPoints<real_t>( pressure_side_offset_points.transpose(), "pressure_side_offset_points");

    return 0;
}

template <class T>
int BladeProfile<T>::computeOffset(T const offsetdist, gsBSpline<T> & suctionSideOffset, gsBSpline<T> & pressureSideOffset, gsKnotVector<T> kvfit, int num_samples)
{
    setKnotVector(kvfit);

    computeOffset(offsetdist, suctionSideOffset, pressureSideOffset, num_samples);

    return 0;
}

template <class T>
int BladeProfile<T>::compute3D(gsBSpline<T> & suctionSide3d, gsBSpline<T> & pressureSide3d, gsKnotVector<T> kvfit, int num_samples3d)
{
    T rr = myCylinderRadius;
    gsVector<T> vec(2);
    gsMatrix<T> parameter_points(1, num_samples3d+1);
    gsMatrix<T> suction_side_curve_points(2, num_samples3d+1);
    gsMatrix<T> pressure_side_curve_points(2, num_samples3d+1);
    gsMatrix<T> suction_side_curve_points3d(num_samples3d+1, 3);
    gsMatrix<T> pressure_side_curve_points3d(num_samples3d+1, 3);
    gsMatrix<T> parameter_points_suction(1, num_samples3d+1);
    gsMatrix<T> parameter_points_pressure(1, num_samples3d+1);

    // Transformation of profile to blade net
    vec(0) = myRotationCenterX;
    vec(1) = myRotationCenterY;
    myBladeProfile2DSuctionSide.translate(-vec);
    myBladeProfile2DSuctionSide.scale(myChordLength);
    myBladeProfile2DSuctionSide.rotate(myAngle);
    myBladeProfile2DPressureSide.translate(-vec);
    myBladeProfile2DPressureSide.scale(myChordLength);
    myBladeProfile2DPressureSide.rotate(myAngle);


    // Computation of point defining 3D blade profile
    for (index_t i = 0; i < num_samples3d; i++) {
        //parameter_points(0,i) = (1-cos(pi/2*i/num_samples3d));
        parameter_points(0,i) = pow((real_t) i/num_samples3d,2);
    }
    parameter_points(0,num_samples3d) = 1;
    //gsInfo << parameter_points << "\n";
    suction_side_curve_points = myBladeProfile2DSuctionSide.eval(parameter_points);
    suction_side_curve_points = suction_side_curve_points.transpose();
    pressure_side_curve_points = myBladeProfile2DPressureSide.eval(parameter_points);
    pressure_side_curve_points = pressure_side_curve_points.transpose();
    //gsInfo << "pressure side curve points: " << pressure_side_curve_points << "\n";
    for (index_t i = 0; i < num_samples3d+1; i++) {
        suction_side_curve_points3d(i,0) = rr*sin(suction_side_curve_points(i,0)/rr);
        suction_side_curve_points3d(i,1) = rr*cos(suction_side_curve_points(i,0)/rr);
        suction_side_curve_points3d(i,2) = suction_side_curve_points(i,1);
        pressure_side_curve_points3d(i,0) = rr*sin(pressure_side_curve_points(i,0)/rr);
        pressure_side_curve_points3d(i,1) = rr*cos(pressure_side_curve_points(i,0)/rr);
        pressure_side_curve_points3d(i,2) = pressure_side_curve_points(i,1);
    }
    //gsInfo << "pressure side curve points 3D: " << pressure_side_curve_points3d << "\n";

    // Computation of tangent vector of 3D blade profile at leading point
    //gsInfo << myBladeProfile2DSuctionSide.coefs() << "\n";
    gsMatrix<T> suction_side_tangent = myBladeProfile2DSuctionSide.deriv(parameter_points);
    //gsInfo << "suction_side_tangent: " << suction_side_tangent.col(0);
    //gsInfo << "Rozmery suction_side_tangent: " << suction_side_tangent.rows() << ", " << suction_side_tangent.cols() << "\n";
    gsMatrix<T> suction_side_leading_tangent3d(1,3);
    suction_side_leading_tangent3d(0) = suction_side_tangent(0, 0);
    suction_side_leading_tangent3d(1) = 0;
    suction_side_leading_tangent3d(2) = suction_side_tangent(1, 0);
    //gsInfo << suction_side_leading_tangent3d << "\n";
    T alfa = acos(suction_side_curve_points3d(0,1)/sqrt(pow(suction_side_curve_points3d(0,0), 2) + pow(suction_side_curve_points3d(0,1), 2)));
    //gsInfo << alfa*180/pi;
    gsMatrix<T> suction_leading_direction(3,1);
    gsMatrix<T> suction_leading_direction2(3,1);
    if (suction_side_curve_points3d(0,1) < 0) {
        suction_leading_direction(0) = -(suction_side_leading_tangent3d(0)*cos(-alfa)-suction_side_leading_tangent3d(1)*sin(-alfa));
        suction_leading_direction(1) = -(suction_side_leading_tangent3d(1)*cos(-alfa)+suction_side_leading_tangent3d(0)*sin(-alfa));
        suction_leading_direction(2) = suction_side_leading_tangent3d(2);
    }
    else {
        suction_leading_direction(0) = suction_side_leading_tangent3d(0)*cos(alfa)-suction_side_leading_tangent3d(1)*sin(alfa);
        suction_leading_direction(1) = suction_side_leading_tangent3d(1)*cos(alfa)+suction_side_leading_tangent3d(0)*sin(alfa);
        suction_leading_direction(2) = suction_side_leading_tangent3d(2);

    }
    suction_leading_direction2 = -suction_leading_direction;
    //gsInfo << "Suction - pocatecni smer na valci: " << suction_leading_direction << "\n";

    // Approximation of 3D blade profiles
    parameter_points_suction = centripetalParameterization(suction_side_curve_points3d);
    //suctionSide3d = curveFittingWithBoundary(suction_side_curve_points3d, parameter_points_suction, kvfit);
    suctionSide3d = curveFittingWithBoundaryAndInputTangent(suction_side_curve_points3d, parameter_points_suction, kvfit, suction_leading_direction);
    gsMatrix<T> control_points_suction = suctionSide3d.coefs();
    parameter_points_pressure = centripetalParameterization(pressure_side_curve_points3d);
    //pressureSide3d = curveFittingWithBoundary(pressure_side_curve_points3d, parameter_points_pressure, kvfit);
    pressureSide3d = curveFittingWithBoundaryAndInputTangent(pressure_side_curve_points3d, parameter_points_pressure, kvfit, suction_leading_direction2);
    gsMatrix<T> control_points_pressure = pressureSide3d.coefs();

    //gsInfo << "Vysledne ridici body - suction:\n" << control_points_suction << "\n";
    //gsInfo << "Vysledne ridici body - pressure:\n" << control_points_pressure << "\n";

    // Making 3D blade profile C1 at leading point
    T l1 = sqrt(pow(control_points_suction(1,0)-control_points_suction(0,0),2) + pow(control_points_suction(1,1)-control_points_suction(0,1),2));
    T l2 = sqrt(pow(control_points_pressure(1,0)-control_points_pressure(0,0),2) + pow(control_points_pressure(1,1)-control_points_pressure(0,1),2));
    //gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction(0,0) << "\t\t\t" << control_points_suction(0,1) << "\n" <<
    //          control_points_suction(1,0) << "\t\t" << control_points_suction(1,1) << "\n" <<
    //          control_points_pressure(0,0) << "\t\t\t" << control_points_pressure(0,1) << "\n" <<
    //          control_points_pressure(1,0) << "\t\t" << control_points_pressure(1,1) << "\n";
    if (l2 > l1) {
        control_points_suction(1,0) = control_points_suction(0,0) + l2*(control_points_suction(1,0) - control_points_suction(0,0))/l1;
        control_points_suction(1,1) = control_points_suction(0,1) + l2*(control_points_suction(1,1) - control_points_suction(0,1))/l1;
    //    gsInfo << "Upravene body:\n" << control_points_suction << "\n";
        suctionSide3d = gsBSpline<T>( kvfit, control_points_suction);
    }
    else {
        control_points_pressure(1,0) = control_points_pressure(0,0) + l1/l2*(control_points_pressure(1,0) - control_points_pressure(0,0));
        control_points_pressure(1,1) = control_points_pressure(0,1) + l1/l2*(control_points_pressure(1,1) - control_points_pressure(0,1));
    //    gsInfo << "Upravene body:\n" << control_points_pressure << "\n";
        pressureSide3d = gsBSpline<T>( kvfit, control_points_pressure);
    }

    l1 = sqrt(pow(control_points_suction(1,0)-control_points_suction(0,0),2) + pow(control_points_suction(1,1)-control_points_suction(0,1),2));
    l2 = sqrt(pow(control_points_pressure(1,0)-control_points_pressure(0,0),2) + pow(control_points_pressure(1,1)-control_points_pressure(0,1),2));
    //gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction(0,0) << "\t\t\t" << control_points_suction(0,1) << "\n" <<
    //         control_points_suction(1,0) << "\t\t" << control_points_suction(1,1) << "\n" <<
    //         control_points_pressure(0,0) << "\t\t\t" << control_points_pressure(0,1) << "\n" <<
    //         control_points_pressure(1,0) << "\t\t" << control_points_pressure(1,1) << "\n";

    //gsInfo << "Vysledne ridici body po uprave - suction:\n" << control_points_suction << "\n";
    //gsInfo << "Vysledne ridici body po uprave - pressure:\n" << control_points_pressure << "\n";

    myBladeProfile3DSuctionSide = suctionSide3d;
    myBladeProfile3DPressureSide = pressureSide3d;

    //gsInfo << "pressure side 3d coeffs: " << pressureSide3d.coefs() << "\n";

    return 0;
}

template <class T>
int BladeProfile<T>::compute3DGuideVane(gsBSpline<T> & suctionSide3d, gsBSpline<T> & pressureSide3d, gsKnotVector<T> kvfit, int num_samples3d)
{

    gsVector<T> vec(2);
    gsMatrix<T> parameter_points(1, num_samples3d+1);
    gsMatrix<T> suction_side_curve_points(2, num_samples3d+1);
    gsMatrix<T> pressure_side_curve_points(2, num_samples3d+1);
    gsMatrix<T> suction_side_curve_points_rot(num_samples3d+1,2);
    gsMatrix<T> pressure_side_curve_points_rot(num_samples3d+1,2);
    gsMatrix<T> suction_side_curve_points_rot2(num_samples3d+1,2);
    gsMatrix<T> pressure_side_curve_points_rot2(num_samples3d+1,2);
    gsMatrix<T> suction_side_curve_points3d(num_samples3d+1, 3);
    gsMatrix<T> pressure_side_curve_points3d(num_samples3d+1, 3);
    gsMatrix<T> suction_side_curve_points3d_cone(num_samples3d+1, 3);
    gsMatrix<T> pressure_side_curve_points3d_cone(num_samples3d+1, 3);
    gsMatrix<T> parameter_points_suction(1, num_samples3d+1);
    gsMatrix<T> parameter_points_pressure(1, num_samples3d+1);

    // Transformation of profile to blade net

    gsVector<T> transvector(2);
    transvector << myCircularMeshSetting(1,0),0;

    //Scaling chord, rotate angleinMesh around origin, translate lengthinMesh
    myBladeProfile2DSuctionSide.scale(myChordLength);
    myBladeProfile2DPressureSide.scale(myChordLength);
    myBladeProfile2DSuctionSide.rotate(myCircularMeshSetting(0,0));
    myBladeProfile2DPressureSide.rotate(myCircularMeshSetting(0,0));
    myBladeProfile2DSuctionSide.translate(transvector);
    myBladeProfile2DPressureSide.translate(transvector);

    //gsInfo << "myBladeProfile2DSuctionSide \n" << myBladeProfile2DSuctionSide.coefs() << "\n";
    //gsInfo << "myBladeProfile2DPressureSide \n" << myBladeProfile2DPressureSide.coefs() << "\n";

    //transformations for blade centre of rotation


    vec(0)=  myRotationCenterX*myChordLength*(math::cos(myCircularMeshSetting(0,0)))-myRotationCenterY*myChordLength*(math::sin(myCircularMeshSetting(0,0)))+ myCircularMeshSetting(1,0);
    vec(1)= myRotationCenterY*myChordLength*(math::cos(myCircularMeshSetting(0,0)))+ myRotationCenterX*myChordLength*(math::sin(myCircularMeshSetting(0,0)));



    //tangent in the leading edge
     gsMatrix<T> suction_side_tangent = myBladeProfile2DSuctionSide.deriv(parameter_points_suction);
     gsVector<T> suction_side_leading_tangent(2);
     suction_side_leading_tangent = suction_side_tangent.row(0);


    // Computation of point defining 3D blade profile
    for (index_t i = 0; i < num_samples3d; i++) {
        //parameter_points(0,i) = (1-cos(pi/2*i/num_samples3d));
        parameter_points(0,i) = pow((real_t) i/num_samples3d,2);
    }
    parameter_points(0,num_samples3d) = 1;



    suction_side_curve_points = myBladeProfile2DSuctionSide.eval(parameter_points);
    suction_side_curve_points = suction_side_curve_points.transpose();
    pressure_side_curve_points = myBladeProfile2DPressureSide.eval(parameter_points);
    pressure_side_curve_points = pressure_side_curve_points.transpose();



    //transformations to circular mesh
    const real_t circleK = math::log(myCircularMeshSetting(2,0)/myCircularMeshSetting(3,0))/myCircularMeshSetting(1,0);
    const real_t E = 2.71828;
    for (index_t i = 0; i < num_samples3d+1; i++){

         suction_side_curve_points_rot(i,0) = myCircularMeshSetting(3,0) * math::pow(E, circleK * suction_side_curve_points(i,0))  *math::cos(circleK *suction_side_curve_points(i,1) );
         pressure_side_curve_points_rot(i,0) = myCircularMeshSetting(3,0) * math::pow(E,circleK * pressure_side_curve_points(i,0))  *math::cos(circleK *pressure_side_curve_points(i,1) );
         suction_side_curve_points_rot(i,1) = myCircularMeshSetting(3,0) * math::pow(E, circleK * suction_side_curve_points(i,0)) *math::sin(circleK *suction_side_curve_points(i,1) );
         pressure_side_curve_points_rot(i,1) = myCircularMeshSetting(3,0) * math::pow(E,circleK * pressure_side_curve_points(i,0))  *math::sin(circleK *pressure_side_curve_points(i,1) );


    }



  //  gsInfo<< "suction_side_before approx \n" << suction_side_curve_points << "\n";
  //   gsInfo<< "pressure_side_before approx \n" << pressure_side_curve_points << "\n";




     real_t vec_norm = math::sqrt(math::pow(myCircularMeshSetting(3,0) * math::pow(E, circleK * vec(0))  *math::cos(circleK *vec(1) ),2)+math::pow(myCircularMeshSetting(3,0) * math::pow(E, circleK * vec(0)) *math::sin(circleK *vec(1) ),2));
//gsInfo<< vec_norm << "\n";
     real_t alpha = - math::acos(myCircularMeshSetting(3,0) * math::pow(E, circleK * vec(0))  *math::cos(circleK *vec(1) )/vec_norm);
//gsInfo<< alpha << "\n";
gsInfo<<  ((myCircularMeshSetting(6,0)+myCircularMeshSetting(7,0))*(math::tan(myCircularMeshSetting(4,0))) /(2))<< "\n";
    for (index_t i = 0; i < num_samples3d+1; i++){
        suction_side_curve_points_rot2(i,0) =  (suction_side_curve_points_rot(i,0)*math::cos(alpha) - suction_side_curve_points_rot(i,1)*math::sin(alpha))  ;
        pressure_side_curve_points_rot2(i,0) =  (pressure_side_curve_points_rot(i,0)*math::cos(alpha) - pressure_side_curve_points_rot(i,1)*math::sin(alpha)) ;
        suction_side_curve_points_rot2(i,1) =  (suction_side_curve_points_rot(i,1)*math::cos(alpha) + suction_side_curve_points_rot(i,0)*math::sin(alpha)) ;
        pressure_side_curve_points_rot2(i,1) =  (pressure_side_curve_points_rot(i,1)*math::cos(alpha) + pressure_side_curve_points_rot(i,0)*math::sin(alpha));
    }
//       suction_side_leading_tangent(0) = - (suction_side_leading_tangent(0)*math::cos(alpha) - suction_side_leading_tangent(1)*math::sin(alpha));
//        suction_side_leading_tangent(1) = - (suction_side_leading_tangent(1)*math::cos(alpha) + suction_side_leading_tangent(0)*math::sin(alpha));

//    real_t beta = 3.14/2 - myCircularMeshSetting(4,0);
//    for (index_t i = 0; i < num_samples3d+1; i++) {
//        suction_side_curve_points3d(i,0) = suction_side_curve_points(i,0);
//        suction_side_curve_points3d(i,1) = suction_side_curve_points(i,1);
//        suction_side_curve_points3d(i,2) = 0;
//        pressure_side_curve_points3d(i,0) = pressure_side_curve_points(i,0);
//        pressure_side_curve_points3d(i,1) = pressure_side_curve_points(i,1);
//        pressure_side_curve_points3d(i,2) = 0;
//    }
    const real_t PI=3.14159;
    real_t beta = -((PI/2) - myCircularMeshSetting(4,0));

    for (index_t i = 0; i < num_samples3d+1; i++) {
        suction_side_curve_points3d(i,1) = -suction_side_curve_points_rot2(i,1);
        suction_side_curve_points3d(i,0) = -suction_side_curve_points_rot2(i,0)*math::sin(beta);
        suction_side_curve_points3d(i,2) =-suction_side_curve_points_rot2(i,0)*math::cos(beta) + myConeSetting(0,2);
        pressure_side_curve_points3d(i,1) = -pressure_side_curve_points_rot2(i,1);
        pressure_side_curve_points3d(i,0) = -pressure_side_curve_points_rot2(i,0)*math::sin(beta);
        pressure_side_curve_points3d(i,2) = -pressure_side_curve_points_rot2(i,0)*math::cos(beta) + myConeSetting(0,2);
    }

    //gsInfo<< "suction_side_before approx \n" << suction_side_curve_points3d << "\n";
    //gsInfo<< "pressure_side_before approx \n" << pressure_side_curve_points3d << "\n";

    real_t Sz = myCircularMeshSetting(5,0);
    real_t z0 = myConeSetting(0,2);
    real_t c = math::tan(math::abs(beta));
    real_t Ax,Ay,Az,frac;

    for (index_t i = 0; i < num_samples3d+1; i++) {


        Ax = suction_side_curve_points3d(i,0);
        Ay = suction_side_curve_points3d(i,1);
        Az = suction_side_curve_points3d(i,2);
        frac =(c*(Sz-z0)*(Sz-z0))/(-c*(Az-Sz)*(Sz-z0)+math::sqrt((Ax*Ax+Ay*Ay)*(Sz-z0)*(Sz-z0)));
        suction_side_curve_points3d_cone(i,0) = Ax*frac;
        suction_side_curve_points3d_cone(i,1) = Ay*frac;
        suction_side_curve_points3d_cone(i,2) = Sz + (Az-Sz)*frac;

        Ax = pressure_side_curve_points3d(i,0);
        Ay = pressure_side_curve_points3d(i,1);
        Az = pressure_side_curve_points3d(i,2);
       frac =(c*(Sz-z0)*(Sz-z0))/(-c*(Az-Sz)*(Sz-z0)+math::sqrt((Ax*Ax+Ay*Ay)*(Sz-z0)*(Sz-z0)));
        pressure_side_curve_points3d_cone(i,0) = Ax*frac;
        pressure_side_curve_points3d_cone(i,1) =  Ay*frac;
        pressure_side_curve_points3d_cone(i,2) =  Sz + (Az-Sz)*frac;

    }


    gsMatrix<T> suction_side_leading_tangent3d(1,3);
    suction_side_leading_tangent3d(0) = suction_side_leading_tangent(0)*math::sin(beta);
    suction_side_leading_tangent3d(1) = suction_side_leading_tangent(1);
    suction_side_leading_tangent3d(2) = suction_side_leading_tangent(0)*math::cos(beta);

gsMatrix<T> suction_side_leading2_tangent3d(1,3);
suction_side_leading2_tangent3d = -suction_side_leading_tangent3d;
    // Approximation of 3D blade profiles

    parameter_points_suction = centripetalParameterization(suction_side_curve_points3d_cone);

  suctionSide3d = curveFittingWithBoundary(suction_side_curve_points3d_cone, parameter_points_suction, kvfit);
    //suctionSide3d = curveFittingWithBoundaryAndInputTangent(suction_side_curve_points3d, parameter_points_suction, kvfit, suction_side_leading_tangent3d);
    gsMatrix<T> control_points_suction = suctionSide3d.coefs();
    parameter_points_pressure = centripetalParameterization(pressure_side_curve_points3d_cone);
    pressureSide3d = curveFittingWithBoundary(pressure_side_curve_points3d_cone, parameter_points_pressure, kvfit);
    //pressureSide3d = curveFittingWithBoundaryAndInputTangent(pressure_side_curve_points3d, parameter_points_pressure, kvfit, suction_side_leading2_tangent3d);
    gsMatrix<T> control_points_pressure = pressureSide3d.coefs();

    gsInfo << "Vysledne ridici body - suction:\n" << control_points_suction << "\n";
    gsInfo << "Vysledne ridici body - pressure:\n" << control_points_pressure << "\n";

    // Making 3D blade profile C1 at leading point
    T l1 = sqrt(pow(control_points_suction(1,0)-control_points_suction(0,0),2) + pow(control_points_suction(1,1)-control_points_suction(0,1),2));
    T l2 = sqrt(pow(control_points_pressure(1,0)-control_points_pressure(0,0),2) + pow(control_points_pressure(1,1)-control_points_pressure(0,1),2));
    gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction(0,0) << "\t\t\t" << control_points_suction(0,1) << "\n" <<
              control_points_suction(1,0) << "\t\t" << control_points_suction(1,1) << "\n" <<
              control_points_pressure(0,0) << "\t\t\t" << control_points_pressure(0,1) << "\n" <<
              control_points_pressure(1,0) << "\t\t" << control_points_pressure(1,1) << "\n";
    if (l2 > l1) {
        control_points_suction(1,0) = control_points_suction(0,0) + l2*(control_points_suction(1,0) - control_points_suction(0,0))/l1;
        control_points_suction(1,1) = control_points_suction(0,1) + l2*(control_points_suction(1,1) - control_points_suction(0,1))/l1;
        gsInfo << "Upravene body:\n" << control_points_suction << "\n";
        suctionSide3d = gsBSpline<T>( kvfit, control_points_suction);
    }
    else {
        control_points_pressure(1,0) = control_points_pressure(0,0) + l1/l2*(control_points_pressure(1,0) - control_points_pressure(0,0));
        control_points_pressure(1,1) = control_points_pressure(0,1) + l1/l2*(control_points_pressure(1,1) - control_points_pressure(0,1));
        gsInfo << "Upravene body:\n" << control_points_pressure << "\n";
        pressureSide3d = gsBSpline<T>( kvfit, control_points_pressure);
    }

    l1 = sqrt(pow(control_points_suction(1,0)-control_points_suction(0,0),2) + pow(control_points_suction(1,1)-control_points_suction(0,1),2));
    l2 = sqrt(pow(control_points_pressure(1,0)-control_points_pressure(0,0),2) + pow(control_points_pressure(1,1)-control_points_pressure(0,1),2));
    gsInfo << "Delky tecnych vektoru a body:\n" << l1 << "\n" << l2 << "\n" << control_points_suction(0,0) << "\t\t\t" << control_points_suction(0,1) << "\n" <<
             control_points_suction(1,0) << "\t\t" << control_points_suction(1,1) << "\n" <<
             control_points_pressure(0,0) << "\t\t\t" << control_points_pressure(0,1) << "\n" <<
             control_points_pressure(1,0) << "\t\t" << control_points_pressure(1,1) << "\n";

    gsInfo << "Vysledne ridici body po uprave - suction:\n" << control_points_suction << "\n";
    gsInfo << "Vysledne ridici body po uprave - pressure:\n" << control_points_pressure << "\n";


    myBladeProfile3DSuctionSide = suctionSide3d;
    myBladeProfile3DPressureSide = pressureSide3d;

    return 0;
}
#endif // BLADEPROFILE_H


