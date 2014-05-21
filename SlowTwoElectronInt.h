#ifndef PANACHE_SLOWTWOELECTRONINT_H
#define PANACHE_SLOWTWOELECTRONINT_H

#include "TwoBodyAOInt.h"
#include "SlowERIBase.h"

namespace panache
{


/*! \ingroup MINTS
 *  \class ERI
 *  \brief Capable of computing two-electron repulsion integrals.
 */
class SlowTwoElectronInt : public TwoBodyAOInt
{
private:
    SlowERIBase sloweri_;

    size_t compute_quartet(int sh1, int sh2, int sh3, int sh4);

public:
    //! Constructor
    SlowTwoElectronInt(const std::shared_ptr<BasisSet> bs1,
                       const std::shared_ptr<BasisSet> bs2,
                       const std::shared_ptr<BasisSet> bs3,
                       const std::shared_ptr<BasisSet> bs4);

    virtual ~SlowTwoElectronInt();

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator&);

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual size_t compute_shell(int, int, int, int);
};

}

#endif //PANACHE_TWOELECTRONINT_H


