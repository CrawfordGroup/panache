#ifndef PANACHE_ERDERI_H
#define PANACHE_ERDERI_H

#include "ERDTwoElectronInt.h"

namespace panache
{

class BasisSet;

class ERDERI : public ERDTwoElectronInt
{
public:
    ERDERI(const std::shared_ptr<BasisSet> bs1,
           const std::shared_ptr<BasisSet> bs2,
           const std::shared_ptr<BasisSet> bs3,
           const std::shared_ptr<BasisSet> bs4,
           bool use_shell_pairs=false);

    virtual ~ERDERI();
};

} // close namespace panache

#endif // PANACHE_ERDERI_H 

