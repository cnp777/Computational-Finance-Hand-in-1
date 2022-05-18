#pragma once
#include "gaussian_header.hpp"
#include <random>

class MTRNG // Random number class
{
public:
    MTRNG(uint_fast32_t seed);                          // Constructor which takes in seed
    MTRNG(std::mt19937 MT_);                            // Constructor which takes in a mt19937 object
    MTRNG(std::vector<double> inputVector_);            // Constructor which takes in a vector
    MTRNG(const MTRNG& originalObject);                 // Copy constructor
    
    double operator()();                                // ()-operator to return a realization of a U(0,1) r.v.
    double operator()(double a, double b);              // ()-operator to return a realization of a U(a,b) r.v
    
    void setSeed(unsigned int n);                       // Sets the seed of the underlying generator

    double genNormal();                                 // Generates one N(0,1) variable
    void genNormal(std::vector<double>& inputVector);   // Fills the input vector with N(0,1)

    size_t getDim() const;                                    // Returns the dimension
    
    
private:
    std::mt19937 MT;    // Random generator
    size_t dim;         // Dimension
};


// Templated Class
template <typename T>
class tempClass
{
public:
    tempClass(unsigned int seed_){randomGenerator.seed(seed_); }             // Constructor which instantiates the generator with the given seed
    tempClass(T rnGenerator):randomGenerator(rnGenerator){ }                 // Constructor which takes in a random number generator
    
    void setSeed(unsigned int seed_){ randomGenerator.seed(seed_); }         // Changes the seed of the underlying generator

    double genUniform(){ return randomGenerator() / (double)(randomGenerator.max()) ;}  // Generates one U(0,1) variable
    
private:
    T randomGenerator;  // Templated random number generator
};

