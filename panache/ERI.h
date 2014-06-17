#ifndef PANACHE_ERI_H
#define PANACHE_ERI_H

#include <memory>

#ifdef USE_LIBINT
#include "LibintERI.h"
#endif

#ifdef USE_LIBINT2
#include "Libint2ERI.h"
#endif

#ifdef USE_SLOWERI
#include "SlowERI.h"
#endif

namespace panache {

inline std::shared_ptr<TwoBodyAOInt> GetERI(shared_ptr<BasisSet> & bs1, shared_ptr<BasisSet> & bs2,
                                            shared_ptr<BasisSet> & bs3, shared_ptr<BasisSet> & bs4)
{
            #ifdef USE_LIBINT
              return std::shared_ptr<TwoBodyAOInt>(new LibintERI(bs1, bs2, bs3, bs4));
            #endif
            #ifdef USE_LIBINT2
              return std::shared_ptr<TwoBodyAOInt>(new Libint2ERI(bs1, bs2, bs3, bs4));
            #endif
            #ifdef USE_SLOWERI
              return std::shared_ptr<TwoBodyAOInt>(new SlowERI(bs1, bs2, bs3, bs4));
            #endif
}

inline std::shared_ptr<TwoBodyAOInt> GetErfERI(double omega,
                                               shared_ptr<BasisSet> & bs1, shared_ptr<BasisSet> & bs2,
                                               shared_ptr<BasisSet> & bs3, shared_ptr<BasisSet> & bs4)
{
            #ifdef USE_LIBINT
              return std::shared_ptr<TwoBodyAOInt>(new LibintErfERI(omega, bs1, bs2, bs3, bs4));
            #endif
            #ifdef USE_LIBINT2
              return std::shared_ptr<TwoBodyAOInt>(new Libint2ErfERI(omega, bs1, bs2, bs3, bs4));
            #endif
            #ifdef USE_SLOWERI
              throw RuntimeError("ErfERI for SlowERI not implemented!");
            #endif
}


} // close namespace panache


#endif
