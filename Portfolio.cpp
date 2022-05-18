#include "Portfolio.h"
#include <cmath>
#include <vector>
#include <iostream>

Portfolio::Portfolio(std::vector<ShortRateBond*> BondObject_, std::vector<double> NotionalValue_):
BondObject(BondObject_), NotionalValue(NotionalValue_){
    
    dim = BondObject_.size();
    
    if (BondObject_.size() != NotionalValue_.size()) {
        std::cout << "Dimensions don't match" << std::endl;
    }
}

double Portfolio::getPrice(){
    double value = 0;
    for (unsigned int i = 0 ; i < dim; i++){
        value += (*BondObject[i]).getPrice()*NotionalValue[i];
    }
    return value;
}

Portfolio::Portfolio()
{
    BondObject = {};
    NotionalValue = {};
    dim = {};
}

void Portfolio::addBond(ShortRateBond* BondObject_, double NotionalValue_)
{
    BondObject.push_back(BondObject_);
    NotionalValue.push_back(NotionalValue_);
    dim += 1;
}

void Portfolio::removeBond(unsigned int n)
{
    if (dim > 0) {
        BondObject.erase(BondObject.begin() + n -1);
        NotionalValue.erase(NotionalValue.begin()+n -1);
        dim -= 1;
    }
    else {
        std::cout << "No elements to remove" << std::endl;
    }
}
