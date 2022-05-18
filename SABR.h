#pragma once
#include <vector>
#include "MTRNG.h"
#include "CallOption.h"

class SABR
{
public:
    SABR(double S0_, double K_, double sigma0_, double alpha_, double beta_, double rho_, double TT_); // Constructor
    double getTT() const;                                                                            
    double impVol();
    double BS_CallPrice();
    double genPath(std::vector<double>& Normals, double TT) const;
    
private:
    double S0;
    double K;
    double sigma0;
    double alpha;
    double beta;
    double rho;
    unsigned long TT;
};

double MC_SABR(SABR sabr, CallOption calloption, MTRNG MT, double nSteps, double nPaths);  // Monte Carlo call pricing function
double PathsRequired(SABR sabr, CallOption calloption, MTRNG MT, double nSteps, double epsilon);
