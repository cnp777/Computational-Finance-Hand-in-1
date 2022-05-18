#include "SABR.h"
#include "gaussian_header.hpp"
#include <cmath>
#include <vector>
#include <iostream>

SABR::SABR(double S0_, double K_, double sigma0_, double alpha_, double beta_, double rho_, double TT_): S0(S0_), K(K_), sigma0(sigma0_), alpha(alpha_), beta(beta_), rho(rho_), TT(TT_)
{
}
double SABR::getTT() const
{
    return TT;
}

double SABR::impVol()
{
    if (alpha == 0) // To avoid division with zero
    {
        std::cout << "Error. Did you mean to call the CEV class instead?" << std::endl;
        return 0; }
    else
    {
    double Fmid = sqrt(S0*K);
    double CF = pow(Fmid, beta);
    double gamma1 = beta/Fmid;
    double gamma2 = -beta*(1.0-beta)/pow(Fmid,2.0);
    double epsilon = TT*pow(alpha, 2.0);
    double zeta = (alpha/(sigma0*(1.0-beta)))*(pow(S0,(1.0-beta))-pow(K,(1.0-beta)));
    double D = log( (sqrt(1.0-2.0*rho*zeta + pow(zeta,2.0) )+zeta-rho) / (1.0-rho) );
    
    return alpha*(log(S0/K)/D) * (1.0 + ( (2.0*gamma2 - pow(gamma1,2.0)+1/pow(Fmid,2.0))/24.0 * pow(sigma0*CF/alpha,2.0) + rho*gamma1/4*sigma0*CF/alpha + (2.0-3.0*pow(rho, 2.0))/24.0)*epsilon);

    }
}

double SABR::BS_CallPrice() 
{
    double impVol = SABR::impVol();
    double d1 = ( log(S0/K) + ( pow(impVol,2.0) / 2.0)*TT) / (impVol*sqrt(TT));
    double d2 = d1-impVol*sqrt(TT);
    
    return S0*normCdf(d1)-K*normCdf(d2);
}

double SABR::genPath(std::vector<double>& Normals, double TT) const
{
    double size = Normals.size()/(double)(2);
    double StepSize = TT/size;
    double NewS = S0;
    double NewVol = sigma0;
    
    for (unsigned int i=0; i < size; i++){
        NewS   = fmax(NewS,0) + NewVol*pow(fmax(NewS,0), beta)*sqrt(StepSize)*Normals[2*i];
        NewVol = NewVol*exp(-0.5*pow(alpha, 2.0)*StepSize +alpha*sqrt(StepSize)*(rho*Normals[2*i]+sqrt(1-pow(rho, 2.0))*Normals[2*i+1]) );
    }
    return NewS;
}

double MC_SABR(SABR sabr, CallOption calloption, MTRNG MT, double nSteps, double nPaths)
{
    std::vector<double> Normals(2*nSteps, 0.0);
    std::vector<double> MCprice(nPaths, 0.0);
    unsigned long TT = sabr.getTT();
    
    for (auto &x : MCprice)
    {
        MT.genNormal(Normals);                                    // Generates 200 new iid normals each loop
        x =  calloption.CallPayoff( sabr.genPath(Normals, TT) );  // Evaluates the payoff of one path with 100 steps
    }
    
    return std::accumulate(MCprice.begin(), MCprice.end()-1, 0.0)/nPaths; // Note: as r=0 we don't discount
}

double PathsRequired(SABR sabr, CallOption calloption, MTRNG MT, double nSteps, double epsilon)
{
    double nPaths = 0;
    
    do {
        nPaths += 1;
        
    } while ( abs( sabr.BS_CallPrice() - MC_SABR(sabr, calloption, MT, nSteps, nPaths) ) / sabr.BS_CallPrice() > epsilon );
    
    return nPaths;
}
