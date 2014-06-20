#ifndef PANACHE_SLOWERI_H
#define PANACHE_SLOWERI_H

#include "SlowTwoElectronInt.h"

namespace panache
{

class BasisSet;

class SlowERI : public SlowTwoElectronInt
{
public:
    SlowERI(const SharedBasisSet bs1,
           const SharedBasisSet bs2,
           const SharedBasisSet bs3,
           const SharedBasisSet bs4);

    virtual ~SlowERI();
};

} // close namespace panache

#endif // PANACHE_SLOWERI_H 

