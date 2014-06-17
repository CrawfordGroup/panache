#ifndef PANACHE_SOLIDHARMONIC_H
#define PANACHE_SOLIDHARMONIC_H


//! \todo change to forward declaration after creating/moving matrix class
#include "Matrix.h"

namespace panache {

namespace SolidHarmonic {

void solidharm(unsigned int l, int m, unsigned int r2, Matrix& coefmat);
void solidharmonic(int l, Matrix &coefmat);

// there ordering here is arbitrary and doesn't have to match the
// basis set ordering
inline int ncart(int l) { return (l>=0)?((((l)+2)*((l)+1))>>1):0; }
inline int npure(int l) { return 2*l+1; }
inline int icart(int a, int b, int c) { return (((((a+b+c+1)<<1)-a)*(a+1))>>1)-b-1; }
inline int ipure(int l, int m) { return m<0?2*-m:(m==0?0:2*m-1); }

}
}

#endif //PANACHE_SOLIDHARMONIC_H
