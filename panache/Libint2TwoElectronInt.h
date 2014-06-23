#ifndef PANACHE_LIBINT2TWOELECTRONINT_H
#define PANACHE_LIBINT2TWOELECTRONINT_H

#include <libint2.h>
#include "TwoBodyAOInt.h"

namespace panache {


class Fjt;
class AOShellCombinationsIterator;

class Libint2TwoElectronInt : public TwoBodyAOInt
{
protected:
    //! Maximum cartesian class size.
    int max_cart_;

    //! Computes the fundamental
    Fjt *fjt_;

    //! Computes the ERIs between four shells.
    size_t compute_quartet(int, int, int, int);

    //! Original shell index requested
    int osh1_, osh2_, osh3_, osh4_;

    //! Were the indices permuted?
    bool p13p24_, p12_, p34_;

    Libint_t * erival_;

public:
    //! Constructor
    Libint2TwoElectronInt(const SharedBasisSet bs1,
                   const SharedBasisSet bs2,
                   const SharedBasisSet bs3,
                   const SharedBasisSet bs4);

    virtual ~Libint2TwoElectronInt();

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator&);

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual size_t compute_shell(int, int, int, int);

};

}

#endif //PANACHE_LIBINT2TWOELECTRONINT_H
