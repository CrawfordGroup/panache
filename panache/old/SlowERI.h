#ifndef PANACHE_SLOWERI_H
#define PANACHE_SLOWERI_H

#include "SlowTwoElectronInt.h"

namespace panache
{

class BasisSet;

class SlowERI : public SlowTwoElectronInt
{
public:
    SlowERI(const std::shared_ptr<BasisSet> bs1,
           const std::shared_ptr<BasisSet> bs2,
           const std::shared_ptr<BasisSet> bs3,
           const std::shared_ptr<BasisSet> bs4);

    virtual ~SlowERI();
};

} // close namespace panache

#endif // PANACHE_SlowERI_H 

