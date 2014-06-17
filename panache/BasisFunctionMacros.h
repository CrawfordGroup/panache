#ifndef PANACHE_BASISFUNCTIONMACROS_H
#define PANACHE_BASISFUNCTIONMACROS_H


/*! \def INT_NCART(am)
    Gives the number of cartesian functions for an angular momentum.
*/
#define INT_NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)
/*! \def INT_PURE(am)
    Gives the number of spherical functions for an angular momentum.
*/
#define INT_NPURE(am) (2*(am)+1)
/*! \def INT_NFUNC(pu,am)
    Gives the number of functions for an angular momentum based on pu.
*/
#define INT_NFUNC(pu,am) ((pu)?INT_NPURE(am):INT_NCART(am))
/*! \def INT_CARTINDEX(am,i,j)
    Computes offset index for cartesian function.
*/
#define INT_CARTINDEX(am,i,j) (((i) == (am))? 0 : (((((am) - (i) + 1)*((am) - (i)))>>1) + (am) - (i) - (j)))
/*! \def INT_ICART(a, b, c)
    Given a, b, and c compute a cartesian offset.
*/
#define INT_ICART(a, b, c) (((((((a)+(b)+(c)+1)<<1)-(a))*((a)+1))>>1)-(b)-1)
/*! \def INT_IPURE(l, m)
    Given l and m compute a pure function offset.
*/
#define INT_IPURE(l, m) ((l)+(m))

#endif //PANACHE_BASISFUNCTIONMACROS_H
