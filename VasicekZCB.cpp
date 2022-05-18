#include "VasicekZCB.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

// Vasicek class
Vasicek::Vasicek(double r_, double kappa_, double theta_, double Sigma_, double h_): r(r_), kappa(kappa_), theta(theta_),Sigma(Sigma_), h(h_)
{
}
double Vasicek::getKappa() const {return kappa;}
double Vasicek::getR() const {return r;}
double Vasicek::getSigma() const {return Sigma;}
double Vasicek::getTheta() const {return theta;}
double Vasicek::getH() const {return h;}

// Vasicek ODEs
double VasicekODE1(double Bt, double kappa){return kappa * Bt + 1;}
double VasicekODE2(double Bt, double kappa, double theta, double Sigma)
{return - kappa * theta * Bt - 0.5 * pow(Sigma,2)*pow(Bt,2);}

// Runge-Kutta methods
double RungeKuttaVas(double r, double T, double theta, double kappa, double Sigma, double h)
{
    std::vector<double> X(2, 0.0); // Vector holding starting B and A values
    double t = T;                  // Starting point
    
    do {
        double K1 = -h*VasicekODE1(X[0], kappa);
        double L1 = -h*VasicekODE2(X[0], kappa, theta, Sigma);
        
        double K2 = -h*VasicekODE1(X[0]+0.5*K1, kappa);
        double L2 = -h*VasicekODE2(X[0]+0.5*L1,kappa, theta, Sigma);
        
        double K3 = -h*VasicekODE1(X[0]+0.5*K2, kappa);
        double L3 = -h*VasicekODE2(X[0]+0.5*L2,kappa, theta, Sigma);
        
        double K4 = -h*VasicekODE1(X[0]+K3, kappa);
        double L4 = -h*VasicekODE2(X[0]+L3,kappa, theta, Sigma);
        
        // Update B-position
        X[0]=X[0]+(K1+2*K2+2*K3+K4)/6;
        // Update A-position
        X[1]=X[1]+(L1+2*L2+2*L3+L4)/6;
        // Update timestep
        t -= h;
        
    } while (t>=0);
    
    return exp(X[1]+X[0]*r);
}

double RungeKuttaVas(Vasicek vasicek, double T)
{
    double kappa = vasicek.getKappa();
    double theta = vasicek.getTheta();
    double Sigma = vasicek.getSigma();
    double r = vasicek.getR();
    double h = vasicek.getH();
    
    std::vector<double> X(2, 0.0); // Vector holding starting B and A values
    double t = T;                  // Starting point
    
    do {
        double K1 = -h*VasicekODE1(X[0], kappa);
        double L1 = -h*VasicekODE2(X[0], kappa, theta, Sigma);
        
        double K2 = -h*VasicekODE1(X[0]+0.5*K1, kappa);
        double L2 = -h*VasicekODE2(X[0]+0.5*L1,kappa, theta, Sigma);
        
        double K3 = -h*VasicekODE1(X[0]+0.5*K2, kappa);
        double L3 = -h*VasicekODE2(X[0]+0.5*L2,kappa, theta, Sigma);
        
        double K4 = -h*VasicekODE1(X[0]+K3, kappa);
        double L4 = -h*VasicekODE2(X[0]+L3,kappa, theta, Sigma);
        
        // Update B-position
        X[0]=X[0]+(K1+2*K2+2*K3+K4)/6;
        // Update A-position
        X[1]=X[1]+(L1+2*L2+2*L3+L4)/6;
        // Update timestep
        t -= h;
        
    } while (t>=0);
    
    return exp(X[1]+X[0]*r);
}

double RungeKuttaVas(Vasicek vasicek, ZCB zcb)
{
    double kappa = vasicek.getKappa();
    double theta = vasicek.getTheta();
    double Sigma = vasicek.getSigma();
    double r = vasicek.getR();
    double h = vasicek.getH();
    double T = zcb.getT();
    
    std::vector<double> X(2, 0.0); // Vector holding starting B and A values
    double t = T;                  // Starting point
    
    do {
        double K1 = -h*VasicekODE1(X[0], kappa);
        double L1 = -h*VasicekODE2(X[0], kappa, theta, Sigma);
        
        double K2 = -h*VasicekODE1(X[0]+0.5*K1, kappa);
        double L2 = -h*VasicekODE2(X[0]+0.5*L1,kappa, theta, Sigma);
        
        double K3 = -h*VasicekODE1(X[0]+0.5*K2, kappa);
        double L3 = -h*VasicekODE2(X[0]+0.5*L2,kappa, theta, Sigma);
        
        double K4 = -h*VasicekODE1(X[0]+K3, kappa);
        double L4 = -h*VasicekODE2(X[0]+L3,kappa, theta, Sigma);
        
        // Update B-position
        X[0]=X[0]+(K1+2*K2+2*K3+K4)/6;
        // Update A-position
        X[1]=X[1]+(L1+2*L2+2*L3+L4)/6;
        // Update timestep
        t -= h;
        
    } while (t>=0);
    
    return exp(X[1]+X[0]*r);
}

// ODE methods 
void Vasicek::solveODE(double T)
{
    std::vector<double> X(2, 0.0);  // Starting B and A values
    solutions.push_back(X);
    
    double t = T;                   // Starting timepoint
    double i = 0;                   // Starting indexpoint
    
    do {
        std::vector<double> Xnew;
        
        double K1 = -h*VasicekODE1(solutions[i][0], kappa); // solutions[i][0] retrieves the i'th B value
        double L1 = -h*VasicekODE2(solutions[i][0], kappa, theta, Sigma);
        
        double K2 = -h*VasicekODE1(solutions[i][0]+0.5*K1, kappa);
        double L2 = -h*VasicekODE2(solutions[i][0]+0.5*L1, kappa, theta, Sigma);
        
        double K3 = -h*VasicekODE1(solutions[i][0]+0.5*K2, kappa);
        double L3 = -h*VasicekODE2(solutions[i][0]+0.5*L2, kappa, theta, Sigma);
        
        double K4 = -h*VasicekODE1(solutions[i][0]+K3, kappa);
        double L4 = -h*VasicekODE2(solutions[i][0]+L3, kappa, theta, Sigma);
        
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
};

std::vector<double> Vasicek::getODE(double T)
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
};


/////////////////////////////////////////////////////
// Zero Coupon Bond class

ZCB::ZCB(double r_, double T_): r(r_), T(T_)
{
}
double ZCB::getR() const {return r;}
double ZCB::getT() const {return T;}


/////////////////////////////////////////////////////
// Vasicek Bond class
VasicekBond::VasicekBond(Vasicek vasicek_, ZCB zcb_): vasicek(vasicek_), zcb(zcb_)
{
}

 double VasicekBond::getPrice()
{
    double r = zcb.getR();
    double T = zcb.getT();
    
    std::vector<double> solution = vasicek.getODE(T);
    
    return exp(solution[1]+solution[0]*r);
}
