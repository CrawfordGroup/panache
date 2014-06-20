#include "ERDERI.h"

namespace panache {

ERDERI::ERDERI(const SharedBasisSet bs1,
         const SharedBasisSet bs2,
         const SharedBasisSet bs3,
         const SharedBasisSet bs4,
         bool use_shell_pairs)
    : ERDTwoElectronInt(bs1, bs2, bs3, bs4, use_shell_pairs)
{
}

ERDERI::~ERDERI()
{
}

} // Namespace
