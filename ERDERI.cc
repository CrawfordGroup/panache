#include "ERDERI.h"

namespace panache {

ERDERI::ERDERI(const std::shared_ptr<BasisSet> bs1,
         const std::shared_ptr<BasisSet> bs2,
         const std::shared_ptr<BasisSet> bs3,
         const std::shared_ptr<BasisSet> bs4,
         bool use_shell_pairs)
    : ERDTwoElectronInt(bs1, bs2, bs3, bs4, use_shell_pairs)
{
}

ERDERI::~ERDERI()
{
}

} // Namespace
