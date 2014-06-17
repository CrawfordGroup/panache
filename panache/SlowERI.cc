#include "SlowERI.h"

namespace panache {

SlowERI::SlowERI(const std::shared_ptr<BasisSet> bs1,
         const std::shared_ptr<BasisSet> bs2,
         const std::shared_ptr<BasisSet> bs3,
         const std::shared_ptr<BasisSet> bs4)
    : SlowTwoElectronInt(bs1, bs2, bs3, bs4)
{
}

SlowERI::~SlowERI()
{
}

} // Namespace
