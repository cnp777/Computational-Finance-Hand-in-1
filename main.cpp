// HandIn 1 by Clara Tørsløv
// Student id :  cnp777
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "MTRNG.h"
#include "Portfolio.h"
#include "SABR.h"

// Somehow gaussian_header.hpp keeps deleting itself. However it is included in the handin1 folder and can be added such that the code runs:
// go to --> File --> Add files to "handin1" --> choose the folder "handin1" --> choose the file "gaussian_header.hpp" --> the code should run


int main() {
    
    std::cout << " Question 1 " << std::endl;
    std::cout << "-------------------------- " << std::endl;
    
    // Problem 1.1
    std::cout << "---  Q1.P1  --- \n";
    
    std::mt19937 MT(1); // A MT object is instantiated with seed 1
    
    std::cout << "Outputs using instantiation of MT object with seed 1" << std::endl;
    std::cout << "Random " << MT() << std::endl;  // outputs one random value
    
    // Problem 1.2
    std::cout << "\n" << "---  Q1.P2  --- \n";
    
    // A MTRNG object is instantiated using the newly defined class
    uint_fast32_t seed = 1;
    MTRNG MTRNG_(seed);
    
    // Problem 1.3
    std::cout << "\n" << "---  Q1.P3  --- \n";

    std::cout << "Realization of a U(0,1) random variable: " << MTRNG_() << std::endl;
    
    // Problem 1.4
    std::cout << "\n" << "---  Q1.P4  --- \n";

    MTRNG_.setSeed(1); // Sets seed = 1
    std::cout << "Realization of U(0,1) with seed 1: " << MTRNG_() << std::endl;
    MTRNG_.setSeed(1); // Resets seed to 1
    std::cout << "Realization of U(0,1) after resetting seed to 1: " << MTRNG_() << std::endl;

    // Problem 1.5
    std::cout << "\n" << "---  Q1.P5  --- \n";
    
    MTRNG_.setSeed(1);
    std::cout << "Realization of U(1,2) with seed 1: " << MTRNG_(1,2) << std::endl;
    
    // Problem 1.6
    std::cout << "\n" << "---  Q1.P6  --- \n";

    // Problem 1.7
    std::cout << "\n" << "---  Q1.P7  --- \n";
    
    MTRNG_.setSeed(1);
    std::cout << "Realization of N(0,1) with seed 1: " << MTRNG_.genNormal() << std::endl;
    
    // Problem 1.8
    std::cout << "\n" << "---  Q1.P8  --- \n";

    std::mt19937 MTRNG10(10);     // Instantiating a mt19937 object to seed 10
    MTRNG MTRNGn(MTRNG10);        // Instantiate class using the mt19937 object
    std::cout << "Realization of N(0,1) with seed 10: " << MTRNGn.genNormal() << std::endl;
    
    // Problem 1.9
    std::cout << "\n" << "---  Q1.P9  --- \n";

    std::vector<double> inputVector(5,0.0); // Creates a vector with dimension 5 and 0 in every place
    MTRNG MTRNGvec(inputVector);            // Creates MTRNG object using the vector constructor
    MTRNGvec.setSeed(1);                    // Sets the seed to 1
    MTRNGvec.genNormal(inputVector);        // Fills 'inputVector' with standard normals
    
    std::cout << "Vector of dimension 5 holding realizations of N(0,1)" << std::endl;
    for (int i=0; i<inputVector.size() ; i++){
        std::cout << inputVector[i] << std::endl;
    }
        
    // Problem 1.10
    std::cout << "\n" << "---  Q1.P10  --- \n";

    MTRNG_.setSeed(1);       // Sets the seed to 1 in a previous instantiation
    MTRNG MTRNGcopy(MTRNG_); // Uses copy constructor

    std::cout << "Realization of U(0,1) using copied object: " << MTRNGcopy() << std::endl;
    
    // Problem 1.11
    std::cout << "\n" << "---  Q1.P11  --- \n";
    
    tempClass<std::mt19937> tClassObj(1);      // Instantiates tempClass object using Mersenne Twister generator and seed 1
    std::cout << "Realization of U(0,1) (using templated class and MT) " << tClassObj.genUniform() << std::endl;
    
    std::mt19937 mt1(1);                       // Instantiate mt19937 object with seed 1
    tempClass<std::mt19937> tClassObj2(mt1);   // Instantiates tempClass object using overloaded constructor on MT object
    std::cout << "Realization of U(0,1) (using overloaded constructor and MT) " << tClassObj2.genUniform() << std::endl;
    
    tempClass<std::minstd_rand> tClassObj3(1); // Instantiates tempClass object using minstd random generator and seed 1
    std::cout << "Realization of U(0,1) (using templated class and minstd) " << tClassObj3.genUniform() << std::endl;
    
    std::minstd_rand minstd1(1);                     // Instantiate minstd_rand object with seed 1
    tempClass<std::minstd_rand> tClassObj4(minstd1); // Instantiates tempClass object using overloaded constructor on minstd_rand object
    std::cout << "Realization of U(0,1) (using overloaded constructor and minstd) " << tClassObj4.genUniform() << std::endl;



    std::cout << "\n" << std::endl;
    std::cout << "\n" << std::endl;

    std::cout << " Question 2 " << std::endl;
    std::cout << "-------------------------- " << std::endl;
    
    // Problem 2.1
    std::cout << "\n" << "---  Q2.P1 --- \n";
    
    // Problem 2.2
    std::cout << "\n" << "---  Q2.P2 --- \n";
    
    // Problem 2.3
    std::cout << "\n" << "---  Q2.P3 --- \n";

    double r = 0.02;
    double kappa = 0.3;
    double theta = 0.03;
    double Sigma = 0.01;
    double T = 1;
    double h = 0.01;

    std::cout << "Bond price using Runge-Kutta function = " << RungeKuttaVas(r, T, theta, kappa, Sigma, h) << std::endl;
    
    // Problem 2.4
    std::cout << "\n" << "---  Q2.P4 --- \n";
    
    Vasicek vasicek(r,kappa,theta,Sigma,h);  // Initialize Vasicek class object using the given parameter values
    
    std::cout << "kappa = " << vasicek.getKappa() << std::endl;
    std::cout << "r = " << vasicek.getR() << std::endl;
    std::cout << "theta = " << vasicek.getTheta() << std::endl;
    std::cout << "Sigma = " << vasicek.getSigma() << std::endl;
    std::cout << "h = " << vasicek.getH() << std::endl;
    
    // Problem 2.5
    std::cout << "\n" << "---  Q2.P5 --- \n";

    std::cout << "Bond price using Runge-Kutta function and Vasicek object = " << RungeKuttaVas(vasicek, T) << std::endl;

    // Problem 2.6
    std::cout << "\n" << "---  Q2.P6 --- \n";
    
    ZCB zcb(r,T); // Initialize ZCB class object
    
    std::cout << "Bond price using Runge-Kutta function, Vasicek object and ZCB object = " << RungeKuttaVas(vasicek, zcb) << std::endl;
    
    // Problem 2.7
    std::cout << "\n" << "---  Q2.P7 --- \n";
    
    vasicek.solveODE(10); // Fills the solutionsmatrix

    std::cout << "Solution to the ODE for T = 10 using getODE() from Vasicek class " << "\n" << "B = " << vasicek.getODE(10)[0] << " A = " << vasicek.getODE(10)[1] << std::endl;
    
    // Problem 2.8
    std::cout << "\n" << "---  Q2.P8 --- \n";
        
    ZCB zcb3(r,3);                          // Instantiate ZCB with 3 years to maturity
    VasicekBond vasicekBond3(vasicek,zcb3); // Instantiate Vasicek Bond object with T = 3
    VasicekBond vasicekBond(vasicek,zcb);   // Instantiate Vasicek Bond object with T = 1
    
    std::cout << "Value of a Vasicek bond 1. 1 year to maturity using the getPrice method " << vasicekBond.getPrice() <<  std::endl;
    std::cout << "Value of a Vasicek bond w. 3 years to maturity using the getPrice method " << vasicekBond3.getPrice() << std::endl;

    // Problem 2.9
    std::cout << "\n" << "---  Q2.P9 --- \n";
    
    // Problem 2.10
    std::cout << "\n" << "---  Q2.P10 --- \n";
     
    std::cout << "CIR bond value using Runge-Kutta method = " << RungeKuttaCIR(r, T, theta, kappa, Sigma, h) << std::endl;
    
    // Problem 2.11
    std::cout << "\n" << "---  Q2.P11 --- \n";

    CIR cir(r, kappa, theta, Sigma, h); // Instantiate CIR object
    cir.solveODE(10);                   // Fill solutionsmatrix
    
    std::cout << "Solution to the ODE for T = 10 using getODE() from CIR class " << "\n" << "B = " << cir.getODE(10)[0] << " A = " << cir.getODE(10)[1] << "\n" << std::endl; 
        
    CIRbond cirBond(cir, zcb);         // Instantiate Cir Bond class object w. T = 1
    CIRbond cirBond3(cir, zcb3);       // Instantiate Cir Bond class object w. T = 3
    
    std::cout << "Value of a CIR bond w. 1 year to maturity using the getPrice method " <<  cirBond.getPrice() << std::endl;
    std::cout << "Value of a CIR bond w. 3 years to maturity using the getPrice method " << cirBond3.getPrice() << std::endl;

    // Problem 2.12
    std::cout <<  "\n" << "---  Q2.P12 --- " << std::endl;
    
    ShortRateBond* const ptrVas = &vasicekBond; // Instantiate an object which refers to a Vasicek Bond object w. T = 1
    ShortRateBond* const ptrCir = &cirBond;     // Instantiate an object which refers to a CIR Bond object w. T = 1
    
    std::vector<ShortRateBond*> ShortRateBondVec({ptrVas,ptrCir});  // Instantiate vector of pointers
    std::vector<double>  NotionalValue(2,100);                      // Instantiate vector of notional values
    Portfolio pf(ShortRateBondVec, NotionalValue);                  // Instantiate portfolio

    std::cout << "The value of the portfolio holding 1 Vasicek Bond and 1 CIR Bond, both with 1 year to maturity and notional value 100, is " << pf.getPrice() << std::endl;

    // Problem 2.13
    std::cout <<  "\n" << "---  Q2.P13 --- " << std::endl;
    
    Portfolio newPF;                              // Instantiate an empty portfolio using the portfolio constructor
    std::cout << "The value of the empty portfolio is " << newPF.getPrice() << std::endl;
    
    ShortRateBond* const ptrVas3 = &vasicekBond3; // Instantiate an object which refers to a Vasicek Bond object w. T = 3
    
    newPF.addBond(ptrCir, 100);                   // Add a 100 notional value CIR bond with T = 1 to the portfolio
    newPF.addBond(ptrVas3, 100);                  // Add a 100 notional value Vasicek bond with T = 3 to the portfolio

    std::cout << "The value of the portfolio holding 1 Vasicek Bond w. T=3 and 1 CIR Bond w. T=1, both with notional value 100, is " << newPF.getPrice() << std::endl;
 
    newPF.removeBond(2);                          // Removes the Vasicek bond
    std::cout << "The value of the portfolio holding 1 CIR Bond w. T=1 and notional value 100 is " << newPF.getPrice() << std::endl;
    
    
    
    std::cout << "\n" << std::endl;
    std::cout << "\n" << std::endl;
    
    std::cout << " Question 3 " << std::endl;
    std::cout << "-------------------------- " << std::endl;
    
    // Problem 3.1
    std::cout << "\n" << "---  Q3.P1 --- \n";

    std::cout << "   x   | normCdf(x) |  theoretical" << std::endl;
    std::cout << "------------------------------------  " << std::endl;
    std::cout << " -1.96 | " << normCdf(-1.96) << "  |  0.025 " << std::endl;
    std::cout << "  0    | " << normCdf(0) <<  "        |  0.5 " << std::endl;
    std::cout << "  1.96 | " << normCdf(1.96) << "   |  0.975 " << "\n\n";
    
    // Problem 3.2
    std::cout << "\n" << "---  Q3.P2 --- \n";
    
    double S0 = 90;
    double K = 100;
    double sigma0 = 1.3;
    double alpha = 0.5;
    double beta = 0.5;
    double rho = -0.5;
    double TT = 3;

    SABR sabr(S0,K,sigma0,alpha,beta,rho,TT); // Instantiate SABR model object 
    
    std::cout << "The value of the call option in the SABR model is " << sabr.BS_CallPrice() << std::endl;
    
    // Problem 3.3
    std::cout << "\n" << "---  Q3.P3 --- \n";

    // Problem 3.4
    std::cout << "\n" << "---  Q3.P4 --- \n";
    
    std::vector<double> Normals(200,0); // Creates a vector with dimension 100 and 0 in every place
    MTRNG Normals1(Normals);            // Instantiate MTRNG object with the right dimensions
    Normals1.setSeed(1);                // Sets the seed to 1
    Normals1.genNormal(Normals);        // Fills the vector with iid realizations of N(0,1)
    
    std::cout << "The value at maturity of a single path with 100 steps using seed 1 is " << sabr.genPath(Normals, TT) << std::endl;
    
    // Problem 3.5
    std::cout << "\n" << "---  Q3.P5 --- \n";
    
    CallOption calloption(K, TT);   // Instantiates Call Option object
    
    // Problem 3.6
    std::cout << "\n" << "---  Q3.P6 --- \n";
    
    MTRNG_.setSeed(1); // Sets the seed of the random number generator to 1
    std::cout << "The payoff of a call option, using seed 1 and MC_SABR with 100 steps and 10.000 paths, is " << MC_SABR(sabr, calloption, MTRNG_, 100, 10000) << std::endl;
    
    // Problem 3.7
    std::cout << "\n" << "---  Q3.P7 --- \n";
    
    std::cout << "Number of paths required, using seed = 1, 100 steps per path, for relative error < 0.001 is " << PathsRequired(sabr, calloption, MTRNG_, 100, 0.001) << std::endl;
    
    std::cout << "\n" << "-------------------------- \n";
    return 0;
}

