#ifndef PANACHE_AOSHELLCOMBINATIONSITERATOR_H
#define PANACHE_AOSHELLCOMBINATIONSITERATOR_H

#include <memory>
#include "AOIntegralsIterator.h"

namespace panache {

class BasisSet;

/*! \ingroup MINTS */
class AOShellCombinationsIterator
{
private:
    struct ShellQuartet {
        int P;
        int Q;
        int R;
        int S;
        bool end_of_PK;
    };

    ShellQuartet current;
    int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
    int usii, usjj, uskk, usll, upk;

    int num_unique_pk;

    bool done;

    std::shared_ptr<BasisSet> bs1_;
    std::shared_ptr<BasisSet> bs2_;
    std::shared_ptr<BasisSet> bs3_;
    std::shared_ptr<BasisSet> bs4_;

public:
    AOShellCombinationsIterator(std::shared_ptr<BasisSet>bs1, std::shared_ptr<BasisSet>bs2,
                              std::shared_ptr<BasisSet>bs3, std::shared_ptr<BasisSet>bs4);
    AOShellCombinationsIterator();
    void init(std::shared_ptr<BasisSet>bs1, std::shared_ptr<BasisSet>bs2,
            std::shared_ptr<BasisSet>bs3, std::shared_ptr<BasisSet>bs4);

    void first();
    void next();
    bool is_done() { return done; }

    int p() const { return current.P; }
    int q() const { return current.Q; }
    int r() const { return current.R; }
    int s() const { return current.S; }
    int end_of_PK() const { return current.end_of_PK; }

    AOIntegralsIterator integrals_iterator();
};


}

#endif //PANACHE_AOSHELLCOMBINATIONSITERATOR_H
