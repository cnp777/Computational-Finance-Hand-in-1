#include "CIR.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

// CIR class
CIR::CIR(double r_, double kappa_, double theta_, double Sigma_, double h_): r(r_), kappa(kappa_), theta(theta_),Sigma(Sigma_), h(h_)
{
}
double CIR::getKappa() const {return kappa;}
double CIR::getR() const {return r;}
double CIR::getSigma() const {return Sigma;}
double CIR::getTheta() const {return theta;}
double CIR::getH() const {return h;}

double CirODE1(double Bt, double kappa, double Sigma){return kappa * Bt - (Sigma * Sigma * Bt * Bt)/2 + 1;}
double CirODE2(double Bt, double kappa, double theta)
{return -kappa * theta * Bt ;}

double RungeKuttaCIR(double r, double T, double theta, double kappa, double Sigma, double h)
{
    std::vector<double> X(2, 0.0); // Vector holding starting B and A values
    double t = T;                  // Starting point
    
    do {
        double K1 = -h*CirODE1(X[0], kappa, Sigma);
        double L1 = -h*CirODE2(X[0], kappa, theta);
        
        double K2 = -h*CirODE1(X[0]+0.5*K1, kappa, Sigma);
        double L2 = -h*CirODE2(X[0]+0.5*L1, kappa, theta);
        
        double K3 = -h*CirODE1(X[0]+0.5*K2, kappa, Sigma);
        double L3 = -h*CirODE2(X[0]+0.5*L2, kappa, theta);
        
        double K4 = -h*CirODE1(X[0]+K3, kappa, Sigma);
        double L4 = -h*CirODE2(X[0]+L3, kappa, theta);
        
        // Update B-position
        X[0]=X[0]+(K1+2*K2+2*K3+K4)/6;
        // Update A-position
        X[1]=X[1]+(L1+2*L2+2*L3+L4)/6;
        // Update timestep
        t -= h;
        
    } while (t>=0);
    
    return exp(X[1]+X[0]*r);
}

void CIR::solveODE(unsigned int T)
{
    std::vector<double> X(2, 0.0);  // Starting B and A values
    solutions.push_back(X);
    
    double t = T;                   // Starting timepoint
    unsigned int i = 0;             // Starting indexpoint
    
    do {
        std::vector<double> Xnew;
        
        double K1 = -h*CirODE1(solutions[i][0], kappa, Sigma);
        double L1 = -h*CirODE2(solutions[i][0], kappa, theta);
        
        double K2 = -h*CirODE1(solutions[i][0]+0.5*K1, kappa, Sigma);
        double L2 = -h*CirODE2(solutions[i][0]+0.5*L1, kappa, theta);
        
        double K3 = -h*CirODE1(solutions[i][0]+0.5*K2, kappa, Sigma);
        double L3 = -h*CirODE2(solutions[i][0]+0.5*L2, kappa, theta);
        
        double K4 = -h*CirODE1(solutions[i][0]+K3, kappa, Sigma);
        double L4 = -h*CirODE2(solutions[i][0]+L3, kappa, theta);
        
        // Update B-position
        Xnew.push_back(solutions[i][0]+(K1+2*K2+2*K3+K4)/6);
        // Update A-position
        Xnew.push_back(solutions[i][1]+(L1+2*L2+2*L3+L4)/6);
        // Update solutions matrix
        solutions.push_back(Xnew);
        // Update timestep and index
        i +=1;
        t -= h;
        
    } while (t>=0);
}

std::vector<double> CIR::getODE(unsigned int T)
{
    std::vector<double> solution;
    
    if (solutions.size()>0) {
    double i =  T/h;
    
    for (unsigned int j=0; j < 2; j++){
        
        solution.push_back(solutions[i][j]);
    }
    }
    else{
        std::cout << "Solutions matrix not initialized yet" << std::endl;
    }
    return solution;
}

/////////////////////////////////////////////////////
// CIR bond class
CIRbond::CIRbond(CIR cir_, ZCB zcb_): cir(cir_), zcb(zcb_)
{
    
}

double CIRbond::getPrice()
{    
    double r = zcb.getR();
    double T = zcb.getT();

    std::vector<double> solution = cir.getODE(T);
    
    return exp(solution[1]+solution[0]*r);
}
