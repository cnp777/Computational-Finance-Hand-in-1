#pragma once
#include <vector>
#include "VasicekZCB.h"
#include "CIR.h"

class Portfolio // Class for constructing and pricing a bond portfolio
{
public:
    Portfolio(std::vector<ShortRateBond*> BondObject_, std::vector<double> NotionalValue_); // Constructor
    Portfolio();                                                                            // Constructs empty portfolio
    double getPrice();                                                                      // Returns value of the portfolio
    void addBond(ShortRateBond* BondObject_, double NotionalValue_);                        // Adds a bond to the portfolio
    void removeBond(unsigned int n);                                                  // Removes the n'th bond in the portfolio
    
private:
    unsigned long dim;                                                                      // Dimension of portfolio vector
    std::vector<ShortRateBond*> BondObject;                                                 // Types of bonds in the portfolio
    std::vector<double> NotionalValue;                                                      // Notional value of each bond in the portfolio
};

