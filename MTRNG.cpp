#include "MTRNG.h"
#include <iostream>
#include <cmath>

MTRNG::MTRNG(uint32_t seed_)
{
    MT.seed(seed_);
}

double MTRNG::operator()()
{
    return MT() / (double)(MT.max());
}

void MTRNG::setSeed(unsigned int n)
{
    MT.seed(n);
}

double MTRNG::operator()(double a, double b)
{
    return a + (b-a)*(MT() / (double)(MT.max()));
}

double MTRNG::genNormal()
{
    double u = operator()();
    return invNormCdf(u);
}

MTRNG::MTRNG(std::mt19937 MT_): MT(MT_)
{
}

MTRNG::MTRNG(std::vector<double> inputVector_): dim(inputVector_.size())
{
}

size_t MTRNG::getDim() const{
    return dim;
}

void MTRNG::genNormal(std::vector<double>& inputVector)
{
    for (auto& x : inputVector){
        x = genNormal();
    }
}

MTRNG::MTRNG(const MTRNG& originalObject): dim(originalObject.dim), MT(originalObject.MT)
{
  
}
