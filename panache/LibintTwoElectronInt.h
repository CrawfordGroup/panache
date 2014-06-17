#ifndef PANACHE_LIBINTTWOELECTRONINT_H
#define PANACHE_LIBINTTWOELECTRONINT_H

#include <libint/libint.h>
#include "TwoBodyAOInt.h"

namespace panache {


class Fjt;
class AOShellCombinationsIterator;

/*! \ingroup MINTS
 *  \class ERI
 *  \brief Capable of computing two-electron repulsion integrals.
 */
class LibintTwoElectronInt : public TwoBodyAOInt
{
protected:
    //! Libint object.
    Libint_t libint_;

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


public:
    //! Constructor
    LibintTwoElectronInt(const std::shared_ptr<BasisSet> bs1,
                         const std::shared_ptr<BasisSet> bs2,
                         const std::shared_ptr<BasisSet> bs3,
                         const std::shared_ptr<BasisSet> bs4);

    virtual ~LibintTwoElectronInt();

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator&);

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual size_t compute_shell(int, int, int, int);

};

}

#endif //PANACHE_LIBINTTWOELECTRONINT_H
