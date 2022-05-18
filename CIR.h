#pragma once
#include "ShortRateBond.h"
#include "VasicekZCB.h"
#include <cmath>
#include <vector>

class CIR // Class which holds the relevant field variables and methods in regards to the CIR Model 
{
public:
    CIR(double r_, double kappa_, double theta_, double Sigma, double h_); // Constructor
    double getKappa() const;
    double getTheta() const;
    double getSigma() const;
    double getR() const;
    double getH() const;
    void solveODE(unsigned int T);
    std::vector<double> getODE(unsigned int T);
    
private:
    double r;
    double kappa;
    double theta;
    double Sigma;
    double h;
    std::vector<std::vector<double> > solutions;
    
};

class CIRbond: public ShortRateBond
{
public:
    CIRbond(CIR cir_, ZCB zcb_); // Constructor
    double getPrice();
    
private:
    CIR cir;
    ZCB zcb;
};

double CirODE1(double Bt, double kappa, double Sigma);
double CirODE2(double Bt, double kappa, double theta);
double RungeKuttaCIR(double r, double T, double theta, double kappa, double Sigma, double h);

