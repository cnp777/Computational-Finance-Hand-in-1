#pragma once
#include <cmath>

//  Zelen and Severo (1964) approximation
// See https://en.wikipedia.org/wiki/Normal_distribution#Numerical_approximations_for_the_normal_CDF
inline double normCdf( double x) // inline means the function can be implemented in a header and still work.
{
    // static constexpr means it will only be set once at compile time. Don't worry about it.
    static constexpr double p = 0.2316419;
    static constexpr double b1 = 0.319381530;
    static constexpr double b2 = -0.356563782;
    static constexpr double b3 = 1.781477937;
    static constexpr double b4 = -1.821255978;
    static constexpr double b5 = 1.330274429;
    
    double y = abs(x);
    double t = 1/(1+p*y);
    double result = 1-1/(2.506628274631)*exp(-pow(y,2)/2)*(b1*t+b2*pow(t,2)+b3*pow(t,3)+b4*pow(t,4)+b5*pow(t,5));
    
    if (x<0){
        result = 1-result;
    }
    
    return result;
}

//    Inverse CDF approximation
//  See Glasserman, Monte Carlo Methods in Financial Engineering, p 68
inline double invNormCdf( double u) // inline means the function can be implemented in a header and still work.
{
    // static constexpr means it will only be set once at compile time. Don't worry about it.
    static constexpr double a0 = 2.50662823884;
    static constexpr double a1 = -18.61500062529;
    static constexpr double a2 = 41.39119773534;
    static constexpr double a3 = -25.44106049637;
    
    static constexpr double b0 = -8.47351093090;
    static constexpr double b1 = 23.08336743743;
    static constexpr double b2 = -21.06224101826;
    static constexpr double b3 = 3.13082909833;
    
    static constexpr double c0 = 0.3374754822726147;
    static constexpr double c1 = 0.9761690190917186;
    static constexpr double c2 = 0.1607979714918209;
    static constexpr double c3 = 0.0276438810333863;
    static constexpr double c4 = 0.0038405729373609;
    static constexpr double c5 = 0.0003951896511919;
    static constexpr double c6 = 0.0000321767881768;
    static constexpr double c7 = 0.0000002888167364;
    static constexpr double c8 = 0.0000003960315187;
    
    double x;
    double r;
    double y = u-0.5;
    if (abs(y)<0.42){
        r = y*y;
        x = y*(((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1);
    }
    else
    {
        r = u;
        if (y>0){r = 1 - u;}
        r = log(-log(r));
            x = c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))));
        if (y<0){x=-x;}
    }
    
    return x;
}


