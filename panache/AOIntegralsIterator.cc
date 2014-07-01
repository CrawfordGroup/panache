
#include <utility> // for std::swap
#include "panache/AOIntegralsIterator.h"
#include "panache/GaussianShell.h"

namespace panache {


AOIntegralsIterator::AOIntegralsIterator(const GaussianShell& s1, const GaussianShell& s2,
                                     const GaussianShell& s3, const GaussianShell& s4)
    : usi(s1), usj(s2), usk(s3), usl(s4)
{
    done = false;
    ni = usi.nfunction();
    nj = usj.nfunction();
    nk = usk.nfunction();
    nl = usl.nfunction();

    fii = usi.function_index();
    fij = usj.function_index();
    fik = usk.function_index();
    fil = usl.function_index();

    iimax = ni - 1;
    if (&usi == &usj && &usk == &usl && &usi == &usk) {
        kkmax = 0;
        llmax = 0;
        jjmax = 0;
    }
    else if(&usi == &usk && &usj == &usl){
        kkmax = 0;
        llmax = 0;
        jjmax = nj - 1;
    }
    else{
        kkmax = nk - 1;
        jjmax = (&usi == &usj) ? 0 : nj - 1;
        llmax = (&usk == &usl) ? 0 : nl - 1;
    }

    ii = 0;
    jj = 0;
    kk = 0;
    ll = 0;
}

void AOIntegralsIterator::first()
{
    current.i = 0 + fii;
    current.j = 0 + fij;
    current.k = 0 + fik;
    current.l = 0 + fil;
    current.index = 0;
    if (&usi == &usj && &usk == &usl && &usi == &usk) {     // (aa|aa) case
    }
    else if(&usi== &usk && &usj == &usl){
        if (current.i < current.j) {
            std::swap(current.i, current.j);
            std::swap(current.k, current.l);
        }
        if (current.i < current.k) {
            std::swap(current.i, current.k);
            std::swap(current.j, current.l);
        }
    }
    else{
        if (current.i < current.j) {
            std::swap(current.i, current.j);
        }
        if (current.k < current.l) {
            std::swap(current.k, current.l);
        }
        if ((current.i < current.k) || (current.i == current.k && current.j < current.l)) {
            std::swap(current.i, current.k);
            std::swap(current.j, current.l);
        }
    }
}

void AOIntegralsIterator::next()
{
    if (&usi == &usj && &usk == &usl && &usi == &usk) {
        ++ll;
        if(ll > llmax){
            ++kk;
            ll = 0;
            if(kk > kkmax){
                kk = 0;
                ++jj;
                if(jj > jjmax){
                    jj = 0;
                    ++ii;
                    if(ii > iimax){
                        done = true;
                    }
                    jjmax = ii;
                }
                kkmax = ii;

            }
            llmax = (kk==ii) ? jj : kk;
        }
        current.i = ii + fii;
        current.j = jj + fij;
        current.k = kk + fik;
        current.l = ll + fil;
        current.index = ll+nl*(kk+nk*(jj+nj*ii));

    }
    else if(&usi == &usk && &usj == &usl){ //(ab|ab)
        ++ll;
        if(ll > llmax){
            ++kk;
            ll = 0;
            if(kk > kkmax){
                kk = 0;
                ++jj;
                if(jj > jjmax){
                    jj = 0;
                    ++ii;
                    if(ii > iimax){
                        done = true;
                    }
                }
                kkmax = ii;
            }
            llmax = (kk == ii) ? jj : nl - 1;
        }
        current.i = ii + fii;
        current.j = jj + fij;
        current.k = kk + fik;
        current.l = ll + fil;
        current.index = ll+nl*(kk+nk*(jj+nj*ii));
        if (current.i < current.j) {
            std::swap(current.i, current.j);
            std::swap(current.k, current.l);
        }
        if (current.i < current.k) {
            std::swap(current.i, current.k);
            std::swap(current.j, current.l);
        }
    }
    else{
        ++ll;
        if(ll > llmax){
            ++kk;
            ll = 0;
            if(kk > kkmax){
                kk = 0;
                ++jj;
                if(jj > jjmax){
                    jj = 0;
                    ++ii;
                    if(ii > iimax){
                        done = true;
                    }
                    jjmax = (&usi == &usj) ? ii : nj - 1;
                }
            }
            llmax = (&usk == &usl) ? kk : nl - 1;
        }
        current.i = ii + fii;
        current.j = jj + fij;
        current.k = kk + fik;
        current.l = ll + fil;
        current.index = ll+nl*(kk+nk*(jj+nj*ii));
        if (current.i < current.j) {
            std::swap(current.i, current.j);
        }
        if (current.k < current.l) {
            std::swap(current.k, current.l);
        }
        if ((current.i < current.k) || (current.i == current.k && current.j < current.l)) {
            std::swap(current.i, current.k);
            std::swap(current.j, current.l);
        }
    }
}

} // close namespace panache

