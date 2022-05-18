#pragma once
#include <cmath>
#include <vector>
#include "ShortRateBond.h"

class Vasicek // Class which holds the relevant field variables and methods in regards to the Vasicek Model
{
public:
    Vasicek(double r_, double kappa_, double theta_, double Sigma_, double h_);    // Constructor
    
    double getKappa() const;
    double getTheta() const;
    double getSigma() const;
    double getR() const;
    double getH() const;
    
    void solveODE(double T);
    std::vector<double> getODE(double T);
    
    
private:
    double r;
    double kappa;
    double theta;
    double Sigma;
    double h;
    std::vector<std::vector<double>> solutions;  // Matrix which holds solutions to solveODE
};

class ZCB // Class which holds the relevant field variables and methods in regards to a Zero Coupon Bond
{
public:
    ZCB(double r, double T); // Constructor
    double getR() const;
    double getT() const;
    
private:
    double r;
    double T;
    
};

class VasicekBond: public ShortRateBond
{
public:
    VasicekBond(Vasicek vasicek_, ZCB zcb_); // Constructor
    double getPrice();
    
private:
    Vasicek vasicek;
    ZCB zcb;
};


double VasicekODE1(double Bt); // dBt/dt = kappa * Bt + 1
double VasicekODE2(double Bt); // dAt/dt = -kappa * theta * Bt - 1/2 * Sigma*Sigma * Bt*Bt
double RungeKuttaVas(double r, double T, double theta, double kappa, double Sigma, double h);
double RungeKuttaVas(Vasicek vasicek, double T); // Overload of pricing function to work of Vasicek model object
double RungeKuttaVas(Vasicek vasicek, ZCB zcb);  // Overload of pricing function to work of Vasicek model object

