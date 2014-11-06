#ifndef PANACHE_FOURINDEXTENSOR_H
#define PANACHE_FOURINDEXTENSOR_H

#include <array>


namespace panache {

struct FourIndexIntegral
{
    int i, j, k, l;
    int perm;
    double integral;

    FourIndexIntegral(int i, int j, int k, int l, int perm, int integral)
                   : i(i), j(j), k(k), l(l), perm(perm), integral(integral)
    { }

    FourIndexIntegral(const std::array<int, 4> & indices, int perm, int integral)
                   : i(indices[0]), j(indices[1]), 
                     k(indices[2]), l(indices[3]),
                     perm(perm), integral(integral)
    { }
};



class FourIndexTensor
{

public:
    typedef std::array<int, 4> IndexArray;

    FourIndexTensor();

    // wrappers for derived class virtual functions
    int GetNBatches(void) const;
    int GetNLocalIntegrals(void) const;
    bool GetNextBatch(void);
    FourIndexIntegral LocalIntegral(int index);


protected:
    virtual int GetNBatches(void) const = 0;
    virtual int GetNLocalIntegrals_(void) const = 0;
    virtual bool GetNextBatch_(void) = 0;
    virtual FourIndexIntegral LocalIntegral_(int index) = 0;

};

}

#endif

