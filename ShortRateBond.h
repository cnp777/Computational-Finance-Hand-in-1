#pragma once

class ShortRateBond  // Abstract base class for bond pricing
{
public:
    virtual double getPrice() =0;
private:
};

