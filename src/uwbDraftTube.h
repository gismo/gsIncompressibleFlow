/** @file uwbDraftTube.h
 * compute Straight draft tube for given parameters

    Author(s): B. Bastl, K. Michalkova
*/

#ifndef UWBDRAFTTUBE_H
#define UWBDRAFTTUBE_H

template <class T>
class DraftTube
{
public:

    DraftTube(int thisNumConeSections, int thisNumDiffuserSections, gsMatrix<T> thisWidthDraftTubeSection, gsMatrix<T> thisHeightDraftTubeSection, gsMatrix<T> thisRadiusDraftTubeSection, gsMatrix<T> thisCentreDraftTubeSection) : myNumConeSections(thisNumConeSections), myNumDiffuserSections(thisNumDiffuserSections),myWidthDraftTubeSection(thisWidthDraftTubeSection), myHeightDraftTubeSection(thisHeightDraftTubeSection), myRadiusDraftTubeSection(thisRadiusDraftTubeSection),myCentreDraftTubeSection(thisCentreDraftTubeSection) {}
    ~DraftTube() {}

    int getNumConeSections() const{ return myNumConeSections; }
    int getNumDiffuserSections() const { return myNumDiffuserSections; }
    gsMatrix<T> getWidthDraftTubeSection() const { return myWidthDraftTubeSection; }
    gsMatrix<T> getHeightDraftTubeSection() const { return myHeightDraftTubeSection; }
    gsMatrix<T> getRadiusDraftTubeSection() const { return myRadiusDraftTubeSection; }
    gsMatrix<T> getCentreDraftTubeSection() const { return myCentreDraftTubeSection; }
    gsBSpline<T> getDraftTubeSection() const {return myDraftTubeSection;}
   gsTensorBSpline<2, T>  getDraftTubeCone() const { return myDraftTubeCone; }
    gsTensorBSpline<2, T> getDraftTubeDiffuser() const { return myDraftTubeDiffuser; }
    gsBSpline<> getDraftTube() const { return myDraftTube; }
    void setNumConeSections(int thisNumConeSections) {myNumConeSections=thisNumConeSections;}
    void setNumDiffuserSections( int thisNumDiffuserSections) {myNumConeSections=thisNumDiffuserSections;}
    void setWidthDraftTubeSection( gsMatrix<T> thisWidthDraftTubeSection) { myWidthDraftTubeSection = thisWidthDraftTubeSection; }
    void setHeightDraftTubeSection( gsMatrix<T> thisHeightDraftTubeSection) { myHeightDraftTubeSection = thisHeightDraftTubeSection; }
    void setRadiusDraftTubeSection( gsMatrix<T> thisRadiusDraftTubeSection) { myRadiusDraftTubeSection = thisRadiusDraftTubeSection; }
    void setCentreDraftTubeSection( gsMatrix<T> thisCentreDraftTubeSection) { myCentreDraftTubeSection = thisCentreDraftTubeSection; }
    void setDraftTubeSection (gsBSpline<T> section) { myDraftTubeSection = section; }
    void setDraftTubeCone(gsTensorBSpline<2, T> section) { myDraftTubeCone = section; }
    void setDraftTubeDiffuser(gsTensorBSpline<2, T> section) { myDraftTubeDiffuser = section; }
    void setDraftTube(gsBSpline<T> section) { myDraftTube = section; }

    int computeDraftTubeCircle(T radius, T centre);
    int computeDraftTubeDiffuser(T  width, T  height, T  radius, T centre);
    int computeDraftTube();


private:
    int myNumConeSections;
    int myNumDiffuserSections;
    gsMatrix<T> myWidthDraftTubeSection;
    gsMatrix<T> myHeightDraftTubeSection;
    gsMatrix<T> myRadiusDraftTubeSection;
    gsMatrix<T> myCentreDraftTubeSection;
    gsBSpline<T> myDraftTubeSection;
    gsTensorBSpline<2, T>  myDraftTubeCone;
    gsTensorBSpline<2, T>  myDraftTubeDiffuser;
    gsBSpline<T> myDraftTube;
};

template<class T>
int DraftTube<T>::computeDraftTubeCircle(T radius, T centre)
{
    real_t const k = 0.55228475; //constant for circle section


        gsKnotVector<T> kvclp(0, 1, 3, 4);//start,end,interior knots, start/end multiplicites of knots1
        kvclp.insert(0.25, 2);
        kvclp.insert(0.5, 2);
        kvclp.insert(0.75, 2);

        gsMatrix<T> coefsclp(13, 3);
        coefsclp << radius, 0, centre,
                    radius, k * radius, centre,
                    k*radius, radius, centre,
                    0, radius , centre,
                    -k*radius, radius, centre,
                    -radius, k * radius, centre,
                    -radius, 0, centre,
                    -radius, - k * radius, centre,
                    -k*radius, - radius, centre,
                    0, - radius, centre,
                    k*radius, - radius, centre,
                    radius, - k * radius, centre,
                    radius, 0, centre;

        gsBSpline<T> section(kvclp, coefsclp);
        setDraftTubeSection(section);

        //gsWriteParaview(section, "section", number_of_points);
    return 0;
}

template<class T>
int DraftTube<T>::computeDraftTubeDiffuser(T  width, T height, T radius, T centre)
{

    real_t k = 0.55228475; //constant for circle section


    gsKnotVector<T> kvclp(0, 1, 7, 4);//start,end,interior knots, start/end multiplicites of knots1

   for(real_t i = 0.125; i <= 0.875; i += 0.125){
        kvclp.insert(i, 2);
    }

        gsMatrix<T> coefsclp(25, 3);
        coefsclp <<  height + radius, - width,   centre,
                     height + radius, - width/3,  centre,
                     height + radius, width/3,   centre,
                     height + radius, width, centre,
                     height + radius, width + k * radius,centre,
                     height + k*radius, width + radius, centre,
                      height, width + radius,  centre,
                      height/3, width + radius, centre,
                     - height/3, width + radius, centre,
                     - height, width + radius, centre,
                     - height-k*radius, width + radius, centre,
                     - height-radius, width + k * radius, centre,
                     - height-radius, width, centre,
                     - height-radius, width/3, centre,
                     - height-radius, -width/3,  centre,
                     - height-radius, -width, centre,
                     - height-radius, -width - k * radius,centre,
                     - height-k*radius, -width - radius,centre,
                     - height, -width - radius,centre,
                     - height/3, -width - radius,centre,
                      height/3, -width - radius, centre,
                      height, -width - radius, centre,
                      height+k*radius, -width - radius,centre,
                     height+radius, -width - k * radius,centre,
                     height+radius, - width,  centre;

        gsBSpline<T> section(kvclp, coefsclp);

        setDraftTubeSection(section);
    return 0;
}

template<class T>
int DraftTube<T>::computeDraftTube(){
    int sections = this->getNumConeSections()+this->getNumDiffuserSections();
    gsBSpline<T>  sectionCurve[sections];
    gsInfo << "Cone initialization \n";


    for(int i=0;i<myNumConeSections;i++){
        computeDraftTubeCircle(myRadiusDraftTubeSection(i),myCentreDraftTubeSection(i));
        sectionCurve[i]=this->getDraftTubeSection();
    }
    gsInfo << "Diffuser initialization \n";

    for(int i=myNumConeSections-1;i<myNumConeSections+myNumDiffuserSections-1;i++){
        computeDraftTubeDiffuser(myWidthDraftTubeSection(i),myHeightDraftTubeSection(i),myRadiusDraftTubeSection(i),myCentreDraftTubeSection(i));
        sectionCurve[i+1]=this->getDraftTubeSection(); //include parameters of last circle in Diffuser
    }
    //CONE
    gsMatrix<T> coefs ((sectionCurve[0].basis()).size() * myNumConeSections, 3);

    for (int i = 0; i< myNumConeSections; i++){
        for (int j = 0; j<sectionCurve[0].basis().size(); j++){
            coefs(i*(sectionCurve[0].basis().size())+j,0)=sectionCurve[i].coef(j,0);
            coefs(i*(sectionCurve[0].basis().size())+j,1)=sectionCurve[i].coef(j,1);
            coefs(i*(sectionCurve[0].basis().size())+j,2)=sectionCurve[i].coef(j,2);
        }
    }

    /*for (int i = 0; i< myNumConeSections; i++){
        for (int j = 0; j<sectionCurve[0].basis().size(); j++){
            gsInfo << i*(sectionCurve[0].basis().size())+j  << coefs(i*(sectionCurve[0].basis().size())+j,0) << " ";
            gsInfo << coefs(i*(sectionCurve[0].basis().size())+j,1) << " ";
            gsInfo << coefs(i*(sectionCurve[0].basis().size())+j,2) << "\n ";
        }
    }*/
    gsKnotVector<T> kv1(0, 1, myNumConeSections-2, 2);
    gsTensorBSplineBasis<2, T> basisCone(sectionCurve[0].basis().knots(),kv1);
    gsTensorBSpline<2, T>  surfaceCone(basisCone, coefs);
    this->setDraftTubeCone(surfaceCone);


    //DIFFUSER
    gsMatrix<T> coefs2(sectionCurve[myNumConeSections+myNumDiffuserSections-2].basis().size() * myNumDiffuserSections, 3);

    for (int i = myNumConeSections; i< myNumConeSections+myNumDiffuserSections; i++){
        for (int j = 0; j<sectionCurve[myNumConeSections+myNumDiffuserSections-2].basis().size(); j++){
            coefs2((i-myNumConeSections)*(sectionCurve[myNumConeSections+myNumDiffuserSections-2].basis().size())+j,0)=sectionCurve[i].coef(j,0);
            coefs2((i-myNumConeSections)*(sectionCurve[myNumConeSections+myNumDiffuserSections-2].basis().size())+j,1)=sectionCurve[i].coef(j,1);
            coefs2((i-myNumConeSections)*(sectionCurve[myNumConeSections+myNumDiffuserSections-2].basis().size())+j,2)=sectionCurve[i].coef(j,2);
        }
    }

    gsKnotVector<T> kv2(0, 1, myNumDiffuserSections-2, 2);
    gsTensorBSplineBasis<2, T> basisDiffuser(sectionCurve[myNumConeSections+myNumDiffuserSections-2].basis().knots(),kv2);
    gsTensorBSpline<2, T>  surfaceDiffuser(basisDiffuser, coefs2);
     this->setDraftTubeDiffuser(surfaceDiffuser);

    return 0;
}



#endif // UWBDRAFTTUBE_H

