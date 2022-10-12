/** @file uwbHydraulicProfile.h
 * compute Hydraulic Profile for given parameters
    Author(s): B. Bastl, K. Michalkova
*/


#ifndef UWBHYDRAULICPROFILE_H
#define UWBHYDRAULICPROFILE_H

#include <iostream>
#include <gismo.h>
#include <math.h>

#include "uwbTurbineUtils.h"
#include "uwbBladeProfile.h"
#include "uwbDraftTube.h"
#include "uwbProfileOptimization.h"

template <class T>
class HydraulicProfile
{
public:

    HydraulicProfile() {}
    HydraulicProfile(gsMatrix<T> thisParametersInner, gsMatrix<T> thisParametersOuter, bool thisSimplify) : myParametersInner(thisParametersInner), myParametersOuter(thisParametersOuter), mySimplify(thisSimplify){}
    ~HydraulicProfile() {}


    gsMatrix<T> getParametersInner() const { return myParametersInner; }
    gsMatrix<T> getParametersOuter() const { return myParametersOuter; }
    bool getSimplify() const { return mySimplify; }
    gsTensorBSpline<2, T> getHydraulicProfileOuter() const { return myHydraulicProfileOuter; }
    gsTensorBSpline<2, T> getHydraulicProfileInner() const { return myHydraulicProfileInner; }
    void setParametersInner( gsMatrix<T> thisParametersInner) { myParametersInner = thisParametersInner; }
    void setParametersOuter( gsMatrix<T> thisParametersOuter) { myParametersOuter = thisParametersOuter; }
    void setSimplify( bool thisSimplify) { mySimplify = thisSimplify; }
    void setHydraulicProfileInner(gsTensorBSpline<2, T> surface) { myHydraulicProfileInner = surface; }
    void setHydraulicProfileOuter(gsTensorBSpline<2, T> surface) { myHydraulicProfileOuter = surface; }

    int computeHydraulicProfile();
    int computeHydraulicProfileInner( gsMatrix<T> params, bool simplify);
    int computeHydraulicProfileOuter( gsMatrix<T> params, bool simplify);
    gsMatrix<T> computeCircleArc(gsMatrix<T> start, T centreX, T centreY, T alpha);
    gsMatrix<T> computeCircle(T centre, T radius);
    gsMatrix<T> computeLine(gsMatrix<T> start, T length, T height, T degree);
    gsMatrix<T> fitLeastSquares(gsMatrix<T> data);



private:

    gsMatrix<T> myParametersInner;
    gsMatrix<T> myParametersOuter;
    bool mySimplify;
    gsTensorBSpline<2, T>  myHydraulicProfileInner;
    gsTensorBSpline<2, T>  myHydraulicProfileOuter;

};

//least squares fitting for given data
template<class T>
 gsMatrix<T> HydraulicProfile<T>::fitLeastSquares(gsMatrix<T> data)
 {

    gsInfo << data << "\n";

    int i,j,k;
    int degree = data.rows();
    gsMatrix<T> sigmaX(2*degree + 1,1);
    sigmaX.setZero();

    for (i=0;i<2*degree+1;i++)
    {
        sigmaX(i,0)=0;
        for (j=0;j<data.rows();j++)
            sigmaX(i,0)=sigmaX(i,0)+math::pow(data(j,0),i);
    }
    gsInfo << "sigmax:\n" << sigmaX << "\n";

    gsMatrix<T> normal(degree+1,degree+2);
    normal.setZero();
    gsMatrix<T> coef(degree+1,1);
    coef.setZero();

    for (i=0;i<=degree;i++)
    {
        for (j=0;j<=degree;j++)
        {
            normal(i,j)=sigmaX(i+j,0);
        }
    }

    gsInfo << "normal:\n" << normal << "\n";

    gsMatrix<T> sigmaY(degree+1,1);
    sigmaY.setZero();

    for (i=0;i<degree+1;i++)
    {
        sigmaY(i,0)=0;
        for (j=0;j<data.rows();j++)
        {
        sigmaY(i,0)=sigmaY(i,0)+math::pow(data(j,0),i)*data(j,1);
        }
    }
    gsInfo << "sigmay:\n" << sigmaY << "\n";

    for (i=0;i<=degree;i++)
    {
        normal(i,degree+1)=sigmaY(i,0);
    }
    degree=degree+1;
    gsInfo << "normal:\n" << normal << "\n";


    for (i=0;i<degree;i++)
    {
        for (k=i+1;k<degree;k++)
        {
            if (normal(i,i)<normal(k,i))
            {
                for (j=0;j<=degree;j++)
                {
                    real_t temp=normal(i,j);
                    normal(i,j)=normal(k,j);
                    normal(k,j)=temp;
                }
            }
        }
    }
    gsInfo << "normal:\n" << normal << "\n";

    gsMatrix<T> m_A(degree, degree);
    m_A.setZero();
    gsMatrix<T> m_B(degree, 1);
    m_B.setZero();

    for (i = 0; i < degree; i++) {
        m_B(i) = normal(i, degree);
        for (j = 0; j < degree; j++) {
            m_A(i, j) = normal(i, j);
        }
    }
    gsInfo << "m_A:\n" << m_A << "\n";
    gsInfo << "m_B:\n" << m_B << "\n";
    //gsMatrix<T> coef (m_B.rows(), m_B.cols());
    coef=m_A.fullPivHouseholderQr().solve( m_B);

    /*
    for (i=0;i<degree-1;i++)
    {
        for (k=i+1;k<degree;k++)
            {
                real_t t=normal(k,i)/normal(i,i);
                for (j=0;j<=degree;j++)
                {
                    //normal(k,j)=normal(k,j)-t*normal(i,j);
                    normal(k,j) -= t*normal(i,j);
                }
            }
    }
    gsInfo << "normal:\n" << normal << "\n";

    for (i=degree-1;i>=0;i--)
    {

        coef(i,0)=normal(i,degree);
        for (j=0;j<degree;j++)
        {
            if (j!=i)
            {
                coef(i,0)=coef(i,0)-normal(i,j)*coef(j,0);
            }
        }
        coef(i,0)=coef(i,0)/normal(i,i);

    }
    */
    gsInfo << "koef fitleastsquares: " << coef << "\n";

    /*gsInfo<<"\n--------------------------------------\n";
    gsInfo<<"leastsquere \n";
    gsInfo<<coef;
    gsInfo<<"\n--------------------------------------\n";*/
    return coef;
    }

template<class T>
 gsMatrix<T> HydraulicProfile<T>::computeCircle(T centre, T radius)
{
    real_t const k = 0.55228475; //constant for circle section

    gsMatrix<T> coefsclp(13, 3);
    coefsclp << centre, 0, radius,
                centre, k * radius, radius,
                centre, radius, k * radius,
                centre, radius , 0,
                centre, radius, - k * radius,
                centre, k * radius, - radius,
                centre, 0, -radius,
                centre, - k * radius, - radius,
                centre, - radius, - k * radius,
                centre, - radius, 0,
                centre, - radius, k * radius,
                centre, - k * radius, radius,
                centre, 0, radius;

    /*gsInfo<<"\n--------------------------------------\n";
    gsInfo<<"circle \n";
    gsInfo<<coefsclp;
    gsInfo<<"\n--------------------------------------\n";*/

    return coefsclp;
}

 template<class T>
  gsMatrix<T> HydraulicProfile<T>::computeLine(gsMatrix<T> start, T length, T height, T degree)
 {
     gsMatrix<T> cp(degree+1,2);
     for(int i=0; i<cp.rows(); i++)
     {
         cp(i,0)= (1.0-i/(degree+0.0))*start(0) + (i/(degree+0.0))*(start(0)+length);
         cp(i,1)=(1.0-i/(degree+0.0))*start(1) + (i/(degree+0.0))*height;
     }

     /*gsInfo<<"\n--------------------------------------\n";
     gsInfo<<"line \n";
     gsInfo<<cp.block(1,0,cp.rows()-1,cp.cols());
     gsInfo<<"\n--------------------------------------\n";*/

    return cp.block(1,0,cp.rows()-1,cp.cols());
 }

  template<class T>
   gsMatrix<T> HydraulicProfile<T>::computeCircleArc(gsMatrix<T> start, T centreX, T centreY, T alpha)
  {
       //constant for circle section

      real_t radius = math::sqrt(math::pow(start(0)-centreX,2)+math::pow(start(1)-centreY,2));
      real_t norm_koef =  - radius * (4.0/3.0)*(math::tan(alpha/4));

      gsMatrix<T> end(2,1);
      end << centreX + (-centreX + start(0))*( math::cos(alpha)) - (-centreY + start(1)) * (math:: sin(alpha)),
             centreY + (-centreY + start(1))*(math::cos(alpha)) + (-centreX + start(0)) *(math::sin(alpha));

      gsMatrix<T> coefsclp(4, 2);
      coefsclp << start(0), start(1),
                  start(0) + (norm_koef*(start(1)-centreY))/(math::sqrt(math::pow(-start(0)+centreX,2)+math::pow(start(1)-centreY,2))), start(1) + (norm_koef*(-start(0)+centreX))/(math::sqrt(math::pow(-start(0)+centreX,2)+math::pow(start(1)-centreY,2))),
                  end(0) + (norm_koef*(-end(1)+centreY))/(math::sqrt(math::pow(-end(1)+centreY,2)+math::pow(end(0)-centreX,2))), end(1) + (norm_koef*(+end(0)-centreX))/(math::sqrt(math::pow(-end(1)+centreY,2)+math::pow(end(0)-centreX,2))),
                  end(0), end(1);

      /*gsInfo<<"\n--------------------------------------\n";
      gsInfo<<"arc \n";
      gsInfo<< coefsclp.block(1,0,coefsclp.rows()-1,coefsclp.cols());
      gsInfo<<"\n--------------------------------------\n";*/

      return coefsclp.block(1,0,coefsclp.rows()-1,coefsclp.cols());
  }



template<class T>
int HydraulicProfile<T>::computeHydraulicProfileInner(gsMatrix<T> params, bool simplify)
{
    gsKnotVector<T> kv;
    gsMatrix<T> coefnew;

    //not simplified case of hydraulic profile
    unsigned degree = 3; //in case of change trouble for the circular parts in generatrix;
    gsMatrix<T> generatrix((degree+1)*2+(params.rows()-2)*degree-(degree+1),3);

    for(int i = 0; i < generatrix.rows();i++)
    {
       for(int j = 0; j < generatrix.cols();j++)
       {
          generatrix(i,j)=0.0;
       }
    }

    generatrix(0,0)=params(0,0);
    generatrix(0,1)=params(0,1);
    gsMatrix<T> start(1,2);
    start = generatrix.block(0,0,1,2);

    for (int i=0; i < params.rows()-1; i++ )
    {
        int part_degree = params(i+1,3);
        switch (part_degree)
        {
            case 1:   //line
                  generatrix.block(i*(degree)+1,0,degree,generatrix.cols()-1) = computeLine(start, params(i+1,0), params(i+1,2), degree);
            break;

            case 2:     //circle
                  generatrix.block(i*(degree)+1,0,degree,generatrix.cols()-1) = computeCircleArc(start, params(i+1,0), params(i+1,1), params(i+1,2));
            break;

            default:
                    gsWarn << "Unexpected part in generatrix.";
            break;
        }

        start = generatrix.block(i*(degree)+degree,0,1,generatrix.cols()-1);
      }

     gsKnotVector<T> kvclp(0, 1, 3, 4); //start,end,interior knots, start/end multiplicites of knots1

     for(real_t i = 0.25; i <= 0.75; i += 0.25)
     {
          kvclp.insert(i, 2);
     }

     gsMatrix<T> par(1,params.rows());
     gsMatrix<T> points_parts(params.rows(),generatrix.cols());

     for (int i=0; i < points_parts.rows(); i++ )
     {
         points_parts.row(i)=generatrix.row(i*degree);
     }

     par = centripetalParameterization(points_parts);
     std::vector<real_t> knots;

     for (int i = 0; i < par.cols();i++)
     {
         knots.push_back(par(0,i));
     }

     gsKnotVector<> kvline(knots, degree, 0); // knots, degree, regularity
     kv = kvline;
     gsMatrix<T> coefs(generatrix.rows()*13,3);

     for(unsigned i = 0; i < kvclp.size() - (degree+1); i++)
     {
        for(int j = 0; j< generatrix.rows(); j++ )
        {
                 coefs(i*(generatrix.rows())+j,2)=computeCircle(generatrix(j,0), generatrix(j,1))(i,0);
                 coefs(i*(generatrix.rows())+j,1)=computeCircle(generatrix(j,0), generatrix(j,1))(i,1);
                 coefs(i*(generatrix.rows())+j,0)=computeCircle(generatrix(j,0), generatrix(j,1))(i,2);
         }
     }

     coefnew=coefs;
     if(simplify)
     {
         gsMatrix<real_t> params_in(7,2);
         params_in << generatrix(15,0),generatrix(15,1),
                    (1060-1200)/1000.0,159.0/1000,
                    (1080-1200)/1000.0,163.0/1000,
                    0.0, generatrix(21,1),
                    (1330-1200)/1000.0,154.0/1000,
                    (1400-1200)/1000.0,130.0/1000,
                    generatrix(30,0),generatrix(30,1);


         params_in << (1010.29 - 1200)/1000.0, 155.0/1000,
                      (1060-1200)/1000.0,159.0/1000,
                      (1080-1200)/1000.0,163.0/1000,
                      0.0, 187.5/1000,
                      (1330-1200)/1000.0,154.0/1000,
                      (1400-1200)/1000.0,130.0/1000,
                      (1449.12-1200)/1000.0, 122.926/1000;

         unsigned degree = params_in.rows();
         gsMatrix<T> generatrixnew(params_in.rows()+1,3);
         generatrixnew.setZero(params_in.rows()+1,3);

         //find coefficients with least square method for the given points
         gsMatrix<T> fit_coef(params_in.rows()+1,1);
         fit_coef.setZero(params_in.rows()+1,1);
         fit_coef = HydraulicProfile::fitLeastSquares(params_in);

         //  fit_coef << 0.1875, 0.028807, -2.43946, -4.12466, 43.538, 54.6389, -264.984, -224.392;

         //reparameterization matrix for given points
         gsMatrix<T> coef_repar(params_in.rows()+1,params_in.rows()+1);
         coef_repar.setZero(params_in.rows()+1,params_in.rows()+1);

         for (int i = 0; i < coef_repar.rows() ;i++ )
         {
             for (int j = 0; j < coef_repar.cols() ;j++ )
             {
                  if(i>j)
                  {
                      coef_repar(i,j) = 0.0;
                  }
                  else
                  {
                      coef_repar(i,j) = (factorial(j)/(factorial(i)*factorial(j-i))) *math::pow(params_in(0,0),j-i)*math::pow((params_in(params_in.rows()-1,0)-params_in(0,0)),i);
                  }

              }
          }

          //repar. coefficients
          gsMatrix<T> fit_coef_repar(params_in.rows()+1,1);
          fit_coef_repar.setZero(params_in.rows()+1,1);
          for(int i = 0; i < coef_repar.rows(); i++)
          {
             for(int j = 0; j < fit_coef.cols(); j++)
             {
                for(int k = 0; k < coef_repar.cols(); k++)
                {
                    fit_coef_repar(i,j) += coef_repar(i,k) * fit_coef(k,j);
                }
             }
          }

          for (int i = 0; i < generatrixnew.rows() ;i++ )
          {
             generatrixnew(i,0) = params_in(0,0) + (i/(generatrixnew.rows()-1.0))*(params_in(params_in.rows()-1,0)-params_in(0,0));
          }

          //coefficient matrix for repar. coefficients
          gsMatrix<T> pom(params_in.rows()+1,params_in.rows()+1);
          for(int i = 0; i < pom.rows(); i++)
          {
              for(int j = 0; j < pom.cols(); j++)
              {
                           if (j>i)
                           {
                              pom(i,j) = 0.0;
                           }
                           else
                           {
                               pom(i,j) = (factorial(i) * factorial(pom.cols()-1-j)) / (factorial(pom.cols()-1) * factorial(i-j));
                           }
              }
          }

          pom << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.142857, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
                 1.0, 0.285714, 0.047619, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.428571, 0.142857, \
                 0.0285714, 0.0, 0.0, 0.0, 0.0, 1.0, 0.571429, 0.285714, 0.114286, \
                 0.0285714, 0.0, 0.0, 0.0, 1.0, 0.714286, 0.47619, 0.285714, 0.142857, \
                 0.047619, 0.0, 0.0, 1.0, 0.857143, 0.714286, 0.571429, 0.428571, \
                 0.285714, 0.142857, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;


         for(int i = 0; i < pom.rows(); i++)
         {
             for(int j = 0; j < fit_coef_repar.cols(); j++)
             {
                 for(int k = 0; k < pom.cols(); k++)
                 {
                     generatrixnew(i,1) += (pom(i,k) * fit_coef_repar(k,j));
                 }
             }
         }

         gsKnotVector<T> kvclp(0, 1, 3, 4);//start,end,interior knots, start/end multiplicites of knots1
         for(real_t i = 0.25; i <= 0.75; i += 0.25)
         {
            kvclp.insert(i, 2);
         }

         gsKnotVector<T> kvline(0, 1, 0, degree + 1);
         kvline.insert(0.25,degree);
         kv = kvline;

         gsMatrix<T> generatrix_in(2*degree+1,2);
         generatrix_in.setZero(2*degree+1,2);
         generatrix_in(0,0)=generatrix(12,0);
         generatrix_in(0,1)=generatrix(12,1);
         generatrix_in.block(1,0,degree,2)=computeLine(generatrix_in.row(0),generatrix(15,0)-generatrix(12,0),generatrix(15,1),degree);
         generatrix_in.block(1+degree,0,degree,2)=generatrixnew.block(1,0,degree,2);

         gsMatrix<T> coefs(15*13,3);

         for(unsigned i = 0; i < kvclp.size() - (3+1); i++)
         {
             for(int j = 0; j< generatrix_in.rows(); j++ )
             {
                coefs(i*(generatrix_in.rows())+j,2)=computeCircle(generatrix_in(j,0), generatrix_in(j,1))(i,0);
                coefs(i*(generatrix_in.rows())+j,1)=computeCircle(generatrix_in(j,0), generatrix_in(j,1))(i,1);
                coefs(i*(generatrix_in.rows())+j,0)=computeCircle(generatrix_in(j,0), generatrix_in(j,1))(i,2);
             }
         }

         coefnew=coefs;
    }//end if simplify

    gsTensorBSplineBasis<2, T> basisInner(kv,kvclp);
    gsTensorBSpline<2, T>  surfaceInner(basisInner, coefnew);
    this->setHydraulicProfileInner(surfaceInner);

    /*gsInfo<<"\n--------------------------------------\n";
    gsInfo<<"profileInner \n";
    gsInfo<<coefnew;
    gsInfo<<"\n--------------------------------------\n";*/

    return 0;
}

template<class T>
int HydraulicProfile<T>::computeHydraulicProfileOuter(gsMatrix<T> params, bool simplify)
{
    gsKnotVector<T> kv;
    gsMatrix<T> coefnew;
    //not simplified case of hydraulic profile
    unsigned degree = 3; //in case of change trouble for the circular parts in generatrix;
    gsMatrix<T> generatrix((degree+1)*2+(params.rows()-2)*degree-(degree+1),3);

    for(int i = 0; i < generatrix.rows();i++)
    {
       for(int j = 0; j < generatrix.cols();j++)
       {
           generatrix(i,j)=0.0;
       }
    }

    generatrix(0,0)=params(0,0);
    generatrix(0,1)=params(0,1);
    gsMatrix<T> start(1,2);
    start = generatrix.block(0,0,1,2);

    for (int i=0; i < params.rows()-1; i++ )
    {
         int part_degree = params(i+1,3);

         switch (part_degree)
         {
             case 1:   //line
                  generatrix.block(i*(degree)+1,0,degree,generatrix.cols()-1) = computeLine(start, params(i+1,0), params(i+1,2), degree);
             break;

             case 2:     //circle
                  generatrix.block(i*(degree)+1,0,degree,generatrix.cols()-1) = computeCircleArc(start, params(i+1,0), params(i+1,1), params(i+1,2));
             break;

             default:
                  gsWarn << "Unexpected part in generatrix.";
             break;
         }

         start = generatrix.block(i*(degree)+degree,0,1,generatrix.cols()-1);
    }

    gsKnotVector<T> kvclp(0, 1, 3, 4); //start,end,interior knots, start/end multiplicites of knots1
    for(real_t i = 0.25; i <= 0.75; i += 0.25)
    {
         kvclp.insert(i, 2);
    }

    gsMatrix<T> par(1,params.rows());
    gsMatrix<T> points_parts(params.rows(),generatrix.cols());
    for (int i=0; i < points_parts.rows(); i++ )
    {
         points_parts.row(i)=generatrix.row(i*degree);
    }
    par = centripetalParameterization(points_parts);

    std::vector<real_t> knots;
    for (int i = 0; i < par.cols();i++)
    {
        knots.push_back(par(0,i));
    }

    gsKnotVector<> kvline(knots, degree, 0); // knots, degree, regularity
    kv = kvline;
    gsMatrix<T> coefs(generatrix.rows()*13,3);
    for(unsigned i = 0; i < kvclp.size() - (degree+1); i++)
    {
       for(int j = 0; j< generatrix.rows(); j++ )
       {
           coefs(i*(generatrix.rows())+j,2)=computeCircle(generatrix(j,0), generatrix(j,1))(i,0);
           coefs(i*(generatrix.rows())+j,1)=computeCircle(generatrix(j,0), generatrix(j,1))(i,1);
           coefs(i*(generatrix.rows())+j,0)=computeCircle(generatrix(j,0), generatrix(j,1))(i,2);
       }
    }

    coefnew = coefs;

    if (simplify)
    {
        unsigned degree = 3;

        gsMatrix<T> line(2*degree+1,2);
        line.block(0,0,degree+1,2) = generatrix.block(12,0,degree+1,2);

        gsMatrix<T>  paramsnew(2,2);
        paramsnew.row(0) = generatrix.block(15,0,1,2);
        paramsnew.row(1) = generatrix.block(27,0,1,2);
        paramsnew(1,1)=0.24912;
        line.block(degree+1,0,degree,2)=computeLine(line.row(degree),paramsnew(1,0),paramsnew(0,1),degree);

        gsKnotVector<T> kvline(0, 1, 0, degree + 1);
        kvline.insert(0.25, 3);
        kv = kvline;

        gsMatrix<T> coefs(7*13,3);
        for(int i = 0; i < coefs.rows();i++)
        {
             for(int j = 0; j < coefs.cols();j++)
             {
                coefs(i,j)=0.0;
             }
        }

        for(unsigned i = 0; i < kvclp.size() - (degree+1); i++)
        {
             for(int j = 0; j< line.rows(); j++ )
             {
                coefs(i*(line.rows())+j,2)=computeCircle(line(j,0), line(j,1))(i,0);
                coefs(i*(line.rows())+j,1)=computeCircle(line(j,0), line(j,1))(i,1);
                coefs(i*(line.rows())+j,0)=computeCircle(line(j,0), line(j,1))(i,2);
             }
        }

        coefnew = coefs;
    }

    gsTensorBSplineBasis<2, T> basisOuter(kv,kvclp);
    gsTensorBSpline<2, T>  surfaceOuter(basisOuter, coefnew);
    this->setHydraulicProfileOuter(surfaceOuter);

    /*gsInfo<<"\n--------------------------------------\n";
    gsInfo<<"profileOuter \n";
    gsInfo<<coefnew;
    gsInfo<<"\n--------------------------------------\n";*/

    return 0;
}

template<class T>
int HydraulicProfile<T>::computeHydraulicProfile()
{
    computeHydraulicProfileInner(myParametersInner, mySimplify);
    computeHydraulicProfileOuter(myParametersOuter, mySimplify);

    return 0;
}



#endif // UWBHYDRAULICPROFILE_H
