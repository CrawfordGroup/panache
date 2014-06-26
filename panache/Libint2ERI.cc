/*! \file
 * \brief Class for calculating libint2-based ERI (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "Libint2ERI.h"
#include "Fjt.h"
#include "BasisSet.h"

namespace panache
{

/////////
// Normal two-electron repulsion integrals
/////////

Libint2ERI::Libint2ERI(const SharedBasisSet bs1,
         const SharedBasisSet bs2,
         const SharedBasisSet bs3,
         const SharedBasisSet bs4)
    : Libint2TwoElectronInt(bs1, bs2, bs3, bs4)
{
    fjt_ = new Taylor_Fjt(basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          1, 1e-15);
}

Libint2ERI::~Libint2ERI()
{
    delete fjt_;
}

/////////
// ErfLibint2ERI
/////////

Libint2ErfERI::Libint2ErfERI(double omega,
               const SharedBasisSet bs1,
               const SharedBasisSet bs2,
               const SharedBasisSet bs3,
               const SharedBasisSet bs4)
    : Libint2TwoElectronInt(bs1, bs2, bs3, bs4)
{
    fjt_ = new ErfFundamental(omega,
                              basis1()->max_am() +
                              basis2()->max_am() +
                              basis3()->max_am() +
                              basis4()->max_am() +
                              1);
}

Libint2ErfERI::~Libint2ErfERI()
{
    delete fjt_;
}

void Libint2ErfERI::setOmega(double omega)
{
    (static_cast<ErfFundamental*>(fjt_))->setOmega(omega);
}

} // close namespace panache
