#ifndef _uwbSwarm_
#define _uwbSswarm_

//#include<vector>
#include <math.h>
#include <gismo.h>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace gismo;

typedef enum {const_w1, const_w2, const_c1, const_c2, const_K} constant;

template<class T>
class opening {

public:
    gsMatrix<T> f_parts, previousSolution;
    gsVector<T> CH, CH_f; // objective function
    gsVector<int> N_iter;
    T f, f_eff;

opening()
{
}

opening(int nProfiles)
{
    initialize(nProfiles);
}

~opening()
{
}

void initialize(int nProfiles){
    f = 1e+19;
    f_eff = 1e+19;
    f_parts.resize(3,nProfiles); f_parts.setZero();
    N_iter.resize(nProfiles);    //N_iter.setZero(nProfiles);
    CH.resize(nProfiles);        CH.setZero(nProfiles);
    CH_f.resize(nProfiles);      CH_f.setZero(nProfiles);
}    
	
};

//=====================================================================================================================================

template <class T>
class swarmMember {
	
private:
    int NParams; // number of parameters
    gsVector<T> p; // parameters with the best value of f
    T fp; // best value of f
	
public:
    gsVector<opening<T>> open;
    gsVector<T> x, v; // parameters, velocities of x
    T f, penaltyF;
    int nP, nO; //number of profiles, number of opens

swarmMember(gsVector<T> xx, int nnP, int nnO)
{
    NParams = xx.size();
    nP = nnP;
    nO = nnO;
    x = xx;
    gsVector<T> vv (NParams); vv.setZero(NParams);
    v = vv; 
    p = x;
    initialize();
}

~swarmMember()
{
}

void initialize(){
    open.resize(nO);
    for (int i = 0; i < nO; i++){
        open(i).initialize(nP);
    }
    f = 1e+19;
    penaltyF = 0.0; 
    fp = f;
}    
	
void actual_best(){ // evaluation of possibly new f and p
    if(f < fp){
        fp = f;
        p = x;
    }
}	 

int get_NParams(){
    return NParams;
}

void set_position(gsVector<T> xx){
    x = xx;
} 
void set_velocity(gsVector<T> vv){
    v = vv;
} 
gsVector<T> get_position(){
    return x;
} 
gsVector<T> get_velocity(){
    return v;
} 

T get_bestValue(){
    return fp; // best value of objective function
}

gsVector<T> get_bestPosition(){
    return p; // parameters of best member
}
};

//========================================== swarm ==========================================

template <class T>
class swarm {
	
private:
	int NParams;
	gsVector<T> g; // best parameters
	T fg; // best value of objective function
	int gi; // index of member with fg
	T w1, w2, c1, c2, K; // constants
        gsVector<T> upperBound, lowerBound;
	T PI = 3.141592653589793238463;

public:
	int NMembers; // number of members of the sworm
	std::vector<swarmMember<T>> member; // pointers to members
	
swarm(){
	NMembers = 0;
	w1 = 0.4;
	w2 = 0.6;
	c1 = 2.05;
	c2 = 2.05;
	T psi = c1 + c2;
	K = 2./fabs(2. - psi - sqrt(psi*psi - 4.*psi));
}

swarm(gsVector<T> w){
	NMembers = 0;
	w1 = w[0];
	w2 = w[1];
	c1 = 2.05;
	c2 = 2.05;
	T psi = c1 + c2;
	K = 2./fabs(2. - psi - sqrt(psi*psi - 4.*psi));
}

swarm(gsVector<T> w, gsVector<T> c){
	NMembers = 0;
	w1 = w[0];
	w2 = w[1];
	c1 = c[0];
	c2 = c[1];
	T psi = c1 + c2;
	K = 2./fabs(2. - psi - sqrt(psi*psi - 4.*psi));
}

swarm(gsVector<T> w, gsVector<T> c, T KK){
	NMembers = 0;
	w1 = w[0];
	w2 = w[1];
	c1 = c[0];
	c2 = c[1];
	K = KK;
}

~swarm(){
//	for(int i = 0; i < NMembers; i++)
//	delete[] member[i];
}
	
void addMember(swarmMember<T> newmember){  // add new member
	member.push_back(newmember);
	NMembers = member.size();

	if(NMembers == 1){ // fisrt member -> initialize g and fg
		fg = newmember.get_bestValue();
		NParams = newmember.get_NParams();
		g = newmember.get_bestPosition();
		gi = 0;
		}
	else if(NParams != newmember.get_NParams()){
		gsInfo << "New swarmMember has different number of parameters." << "\n";
		while (getchar() != '\n');
		exit(0);
	}
}

void bestMember(){  // evaluation of possibly new fg and g
    for(int i = 0; i < NMembers; i ++){
        if(member[i].get_bestValue() < fg){
	    fg = member[i].get_bestValue();
            g = member[i].get_bestPosition();
            gi = i;
        }
    }
}

void next_gen(){  // evolution step
    for(int i = 0; i < NMembers; i++){
        member[i].actual_best();  // actualizaton of p and fp of every member
    }
    bestMember(); // actualizaton of g and fg
    T possibleV;
    for(int i = 0; i < NMembers; i++){
        for(int j = 0; j < NParams; j++){
            possibleV = K*(member[i].v[j] + c1*w1*(  member[i].get_bestPosition()[j] - member[i].x[j]  ) + c2*w2*(  g[j] - member[i].x[j] )  );
            member[i].v[j] = std::max(std::min(possibleV,upperBound[j] - member[i].x[j]),lowerBound[j] - member[i].x[j]);
            member[i].x[j] = member[i].x[j] + member[i].v[j];
        }
    }
}	

T get_bestValue(){
	return fg; // best value of objective function
}

gsVector<T> get_bestPosition(){
	return g; // parameters of best member
}

int get_bestIndex(){
	return gi; // index of best member
}

void set_constant(constant con, T value){
	switch(con){
		case const_w1:
			w1 = value; break;
		case const_w2:
			w2 = value; break;
		case const_c1:
			c1 = value; break;
		case const_c2:
			c2 = value; break;
		case const_K:
			K = value; break;
		default: 
			gsInfo << "Class swarm does not contain this constant. Possibilities are: const_w1, const_w2, const_w3, const_w4 or const_K" << "\n";
			while (getchar() != '\n');
			exit(0);
	}
}
	
T get_constant(constant con){
	switch(con){
		case const_w1:
			return w1;
		case const_w2:
			return w2;
		case const_c1:
			return c1;
		case const_c2:
			return c2;
		case const_K:
			return K;
		default: 
			gsInfo << "Class swarm does not contain this constant. Possibilities are: const_w1, const_w2, const_w3, const_w4 or const_K" << "\n";
			while (getchar() != '\n');
			exit(0);
	}
}

void initBounds(int N){
    gsVector<T> initVectorL(N), initVectorU(N);
    for(int i = 0; i < N; i++){
        initVectorL[i] = -1e+300;
        initVectorU[i] = 1e+300;
    }
    upperBound = initVectorU;
    lowerBound = initVectorL;
}

void set_upperBound(int i, T value){
    upperBound[i] = value;
}

void set_lowerBound(int i, T value){
    lowerBound[i] = value;
}

T get_upperBound(int i){
    return upperBound[i];
}

T get_lowerBound(int i){
    return lowerBound[i];
}

T randomNumber(T lower, T upper)
{
    return (upper - lower) * ( (double)rand() / (double)RAND_MAX ) + lower;
}

gsVector<T> randVelocity(gsVector<T> designP) // setting of random velocity
{
    int nPar = designP.size();
    gsVector<T> velocity(nPar), geometryParams(nPar); 
    T lB, uB, lV, uV; //lower and upper bounds of design parameter and velocity
    for(int i=0; i<nPar; i++){
        lB = lowerBound(i);
        uB = upperBound(i);
        lV = lB - designP[i];
        uV = uB - designP[i];
        velocity[i] = 0.5*randomNumber(lV, uV);
    }    
    return velocity;
}

gsVector<T> randMember(gsVector<T> designP, T toler) // random member of the swarm
{
    int nPar = designP.size();
    gsVector<T> geometryParams(nPar);
    //toler relative tolerance to the boundary, values in <0, 0.5>
    T lB, uB, sizeInter; //lower and upper bounds of design parameter
    for(int i=0; i<nPar; i++){
        lB = lowerBound(i);
        uB = upperBound(i);
        sizeInter = uB - lB;
        lB += toler*sizeInter;
        uB -= toler*sizeInter;
        geometryParams[i] = randomNumber(lB, uB);
    }
    return geometryParams;
}

void trialMembers(int nM, gsVector<T> x0, T limit, const int nP, const int nO)
{
    addMember(swarmMember<T>(x0, nP, nO));
    member[0].penaltyF = penaltyConstraints7(x0, nP);
    member[0].v = randVelocity(x0);

    int nTrial = 0; //counter of trial members
    for (int i = 1; i < nM; i++){    
        T H = limit + 1.0;
        gsVector<T> x (x0.size());
        while ((H > limit) && (nTrial < 1e+4)){  // boundary of penalty constraints to choose trial member as swarm.member
            x = randMember(x0,0.01);
            H = penaltyConstraints7(x, nP);
            nTrial += 1;
        }
        if (H <= limit){
            gsInfo << "trial member number: " << nTrial << "\n";
            addMember(swarmMember<T>(x, nP, nO));
            member[i].penaltyF = H;
            member[i].v = randVelocity(x);
        } else {
            limit *= 2.;
            nTrial = 0;
            i -= 1;
        }
    }
}

T penaltyConstraints7(const gsVector<T> u, const int nP)
{ 
    gsVector<T> constr;
    if (nP == 7){ 
        constr.resize(5*nP+9*(nP-1));
    } else {
        constr.resize(5*nP);
    }
    for (int IOP = 0; IOP < nP; IOP++){
        int j = 0;
        T camber_x = u(j+IOP); j += nP;
        T camber_y = u(j+IOP); j += nP;
        T leading_angle = u(j+IOP); j += nP;
        T trailing_angle = u(j+IOP); j += nP;
        T thickness_x = u(j+IOP); j += nP;
        T thickness_y = u(j+IOP); j += nP;
        T output_angle = u(j+IOP); j += nP;
        T radius = u(j+IOP); j += 2*nP;
        T endingOffset = u(j+IOP);
        constr[0+IOP*5] = math::tan(leading_angle) - camber_y/camber_x;
        constr[1+IOP*5] = math::tan(trailing_angle) - camber_y/(1.0-camber_x);
        constr[2+IOP*5] = math::tan(output_angle) - (thickness_y - endingOffset)/(1.0-thickness_x);
        constr[3+IOP*5] = thickness_y/2.0 - radius; 
        constr[4+IOP*5] = -leading_angle - trailing_angle + PI/2.0;// constr[4] = leading_angle + trailing_angle;
    }
    
    if (nP == 7){
        int i = 0;
        int k = 0; //camber_x
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        k = 1; // camber_y
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;    
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        k = 2; //leading_angle
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        k = 3; //trailing_angle
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;    
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        k = 4; // thickness_x    
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        k = 5; // thickeness_y
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;    
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        k = 6; // output_angle
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;    
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        k = 7; // radius    
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        constr[i+5*nP] = u(i+k) - u(i+k+1); i++;
        k = 8; // angle
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
        constr[i+5*nP] = -u(i+k) + u(i+k+1); i++;
    }

    int gamma;
    T theta, q;
    T H = 0.0;

    for (int i = 0; i < constr.size(); i++){
        q = std::max(0.0,-constr[i]);
        gamma = 1;
        if (-constr[i] < 0.01){
            theta = 1.0;
        } else if (-constr[i] < 0.1){
            theta = 10.0;
        } else {
            theta = 100;
        }            
        H += theta*std::pow(q,gamma);
    }
    return H;
}

};	


#endif
