#pragma once

class CallOption{
public:
    CallOption(double strike_, double maturity_); // Constructor
    double getStrike() const;
    double getMaturity() const;
    double CallPayoff(double spot) const;         // Evaluates the payoff of a call option
    
private:
    double strike;
    double maturity;
};
