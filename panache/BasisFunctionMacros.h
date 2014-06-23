/*! \file
 * \brief Defines some macros for basis functions
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef PANACHE_BASISFUNCTIONMACROS_H
#define PANACHE_BASISFUNCTIONMACROS_H


/*! 
 * \def INT_NCART(am)
 * \brief Gives the number of cartesian functions for an angular momentum.
 * 
 * \param am Angular momentum of the shell
 */
#define INT_NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)


/*!
 * \def INT_NPURE(am)
 * \brief Gives the number of spherical functions for an angular momentum.
 * 
 * \param am Angular momentum of the shell
 */
#define INT_NPURE(am) (2*(am)+1)


/*!
 * \def INT_NFUNC(pu,am)
 * \brief Gives the number of functions for an angular momentum
 *
 * \param pu Is the shell pure (spherical) or ont
 * \param am Angular momentum of the shell
 */
#define INT_NFUNC(pu,am) ((pu)?INT_NPURE(am):INT_NCART(am))


/*! 
 * \def INT_CARTINDEX(am,i,j)
 * \brief Computes offset index for cartesian function with a shell.
 *
 * The exponent on z in the basis function is inferred.
 * For example, for a d shell (\p am = 2)
 *
 * |Function| i | j | k | INT_CARTINDEX |
 * |--------|---|---|---|---------------|
 * | XX     | 2 | 0 | 0 | 0             |
 * | XY     | 1 | 1 | 0 | 1             |
 * | XZ     | 1 | 0 | 1 | 2             |
 * | YY     | 0 | 2 | 0 | 3             |
 * | YZ     | 0 | 1 | 1 | 4             |
 * | ZZ     | 0 | 0 | 2 | 5             |
 *
 * \param am Angular momentum of the shell
 * \param i The exponent on x in the basis function
 * \param j The exponent on y in the basis function
 */
#define INT_CARTINDEX(am,i,j) (((i) == (am))? 0 : (((((am) - (i) + 1)*((am) - (i)))>>1) + (am) - (i) - (j)))




#endif //PANACHE_BASISFUNCTIONMACROS_H
