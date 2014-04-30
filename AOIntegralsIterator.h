#ifndef PANACHE_AOINTEGRALSITERATOR_H
#define PANACHE_AOINTEGRALSITERATOR_H

namespace panache {

class GaussianShell;

/*! \ingroup MINTS */
class AOIntegralsIterator
{
private:
    struct Integral {
        int i;
        int j;
        int k;
        int l;
        unsigned int index;
    };

    Integral current;
    const GaussianShell& usi;
    const GaussianShell& usj;
    const GaussianShell& usk;
    const GaussianShell& usl;

    bool done;

    int ii, iimax, jj, jjmax, kk, kkmax, ll, llmax;
    int ni, nj, nk, nl, fii, fij, fik, fil;

public:
    AOIntegralsIterator(const GaussianShell& s1, const GaussianShell& s2,
                      const GaussianShell& s3, const GaussianShell& s4);

    void first();
    void next();
    bool is_done() { return done; }

    int i() const { return current.i; }
    int j() const { return current.j; }
    int k() const { return current.k; }
    int l() const { return current.l; }
    int index() const { return current.index;}
};

}

#endif //PANACHE_AOINTEGRALSITERATOR_H
