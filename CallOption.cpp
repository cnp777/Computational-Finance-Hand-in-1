#include "CallOption.h"

CallOption::CallOption(double strike_, double maturity_): strike(strike_),maturity(maturity_)
{
}

double CallOption::getStrike() const
{
    return strike;
}

double CallOption::getMaturity() const
{
    return maturity;
}

double CallOption::CallPayoff(double spot) const
{
    double payoff = 0;
    
    if (spot > strike) {
        payoff = spot - strike;
    }
    
    return payoff;
}
