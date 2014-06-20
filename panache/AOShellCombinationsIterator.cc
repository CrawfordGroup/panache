
#include "BasisSet.h"
#include "AOShellCombinationsIterator.h"

namespace panache {

// ===========================================================================
//  AOShellCombinationsIterator
// ===========================================================================
AOShellCombinationsIterator::AOShellCombinationsIterator(SharedBasisSet bs1, SharedBasisSet bs2,
                                                         SharedBasisSet bs3, SharedBasisSet bs4) :
    bs1_(bs1), bs2_(bs2), bs3_(bs3), bs4_(bs4)
{

}

AOShellCombinationsIterator::AOShellCombinationsIterator()
{

}

AOIntegralsIterator AOShellCombinationsIterator::integrals_iterator()
{
    return AOIntegralsIterator(bs1_->shell(p()), bs2_->shell(q()), bs3_->shell(r()), bs4_->shell(s()));
}

void AOShellCombinationsIterator::init(SharedBasisSet bs1, SharedBasisSet bs2,
                                     SharedBasisSet bs3, SharedBasisSet bs4)
{
    bs1_=bs1;
    bs2_=bs2;
    bs3_=bs3;
    bs4_=bs4;
}

void AOShellCombinationsIterator::first()
{
    usii = usjj = uskk = usll = upk = 0;
    done = false;

    num_unique_pk = 1;
    usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;

    int usi, usj, usk, usl;
    usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];

    // Sort shells based on AM, save ERI some work doing permutation resorting.
    if (bs1_->shell(usi).am() < bs2_->shell(usj).am()) {
        std::swap(usi, usj);
    }
    if (bs3_->shell(usk).am() < bs4_->shell(usl).am()) {
        std::swap(usk, usl);
    }
    if (bs1_->shell(usi).am() + bs2_->shell(usj).am() >
            bs3_->shell(usk).am() + bs4_->shell(usl).am()) {
        std::swap(usi, usk);
        std::swap(usj, usl);
    }

    current.P = usi; current.Q = usj; current.R = usk; current.S = usl; current.end_of_PK = false;

    if (upk == num_unique_pk - 1) {
        // If this is the last unique shell flag it as end of a pk block.
        current.end_of_PK = true;
    }
    else{
        current.end_of_PK = false;
    }

}

void AOShellCombinationsIterator::next()
{
    ++upk;
    if(upk >= num_unique_pk){
        upk = 0;
        ++usll;
        if (usll > uskk){
            ++uskk;
            usll = 0;
            if(uskk > usjj){
                ++usjj;
                uskk = 0;
                if(usjj > usii){
                    ++usii;
                    usjj = 0;
                    if(usii >= bs1_->nshell()){
                        done = true;
                        return;
                    }
                }
            }
        }
        usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;
        if ((usii == usjj && usii == uskk) || (usjj == uskk && usjj == usll))
            num_unique_pk = 1;
        else if (usii == uskk || usjj == usll) {
            num_unique_pk = 2;
            usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
        }
        else if (usjj == uskk) {
            num_unique_pk = 2;
            usi_arr[1] = usii; usj_arr[1] = usll; usk_arr[1] = usjj; usl_arr[1] = uskk;
        }
        else if (usii == usjj || uskk == usll) {
            num_unique_pk = 2;
            usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
        }
        else {
            num_unique_pk = 3;
            usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
            usi_arr[2] = usii; usj_arr[2] = usll; usk_arr[2] = usjj; usl_arr[2] = uskk;
        }
    }



    int usi, usj, usk, usl;
    usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];


    // Sort shells based on AM, save ERI some work doing permutation resorting.
    if (bs1_->shell(usi).am() < bs2_->shell(usj).am()) {
        std::swap(usi, usj);
    }
    if (bs3_->shell(usk).am() < bs4_->shell(usl).am()) {
        std::swap(usk, usl);
    }
    if (bs1_->shell(usi).am() + bs2_->shell(usj).am() >
            bs3_->shell(usk).am() + bs4_->shell(usl).am()) {
        std::swap(usi, usk);
        std::swap(usj, usl);
    }

    current.P = usi; current.Q = usj; current.R = usk; current.S = usl; current.end_of_PK = false;

    if (upk == num_unique_pk - 1) {
        // If this is the last unique shell flag it as end of a pk block.
        current.end_of_PK = true;
    }
    else{
        current.end_of_PK = false;
    }

}

} // close namespace panache
