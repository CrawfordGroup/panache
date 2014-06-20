#ifndef PANACHE_ERDERI_H
#define PANACHE_ERDERI_H

#include "ERDTwoElectronInt.h"

namespace panache
{

class BasisSet;

class ERDERI : public ERDTwoElectronInt
{
public:
    ERDERI(const SharedBasisSet bs1,
           const SharedBasisSet bs2,
           const SharedBasisSet bs3,
           const SharedBasisSet bs4,
           bool use_shell_pairs=false);

    virtual ~ERDERI();
};

} // close namespace panache

#endif // PANACHE_ERDERI_H 

