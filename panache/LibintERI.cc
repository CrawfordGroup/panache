/*! \file
 * \brief Class for calculating libint1-based ERI (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/LibintERI.h"
#include "panache/Fjt.h"
#include "panache/BasisSet.h"

namespace panache
{

/////////
// Normal two-electron repulsion integrals
/////////

LibintERI::LibintERI(const SharedBasisSet bs1,
         const SharedBasisSet bs2,
         const SharedBasisSet bs3,
         const SharedBasisSet bs4)
    : LibintTwoElectronInt(bs1, bs2, bs3, bs4)
{
    fjt_ = new Taylor_Fjt(basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          1, 1e-15);
}

LibintERI::~LibintERI()
{
    delete fjt_;
}

/////////
// LibintErfERI
/////////

LibintErfERI::LibintErfERI(double omega,
               const SharedBasisSet bs1,
               const SharedBasisSet bs2,
               const SharedBasisSet bs3,
               const SharedBasisSet bs4)
    : LibintTwoElectronInt(bs1, bs2, bs3, bs4)
{
    fjt_ = new ErfFundamental(omega,
                              basis1()->max_am() +
                              basis2()->max_am() +
                              basis3()->max_am() +
                              basis4()->max_am() +
                              1);
}

LibintErfERI::~LibintErfERI()
{
    delete fjt_;
}

void LibintErfERI::setOmega(double omega)
{
    (static_cast<ErfFundamental*>(fjt_))->setOmega(omega);
}

} // close namespace panache
