/*! \file
 * \brief Holds information about a basis set shell (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_GAUSSIANSHELL_H
#define PANACHE_GAUSSIANSHELL_H

#include "ShellInfo.h"

namespace panache {

/*!
 *  \brief Holds information about a Gaussian orbital shell.
 */
class GaussianShell
{
private:
    
    int l_;              //!< Angular momentum
    int puream_;         //!< Flag for pure angular momentum
    const double* exp_;  //!< Exponents (of length nprimitives_)
    const double* original_coef_;  //!< Original (un-normalized) contraction coefficients (of length nprimitives)
    const double* coef_; //!< Contraction coefficients (of length nprimitives_)

    
    int nc_;  //!< Atom number this shell goes to. Needed when indexing integral derivatives.
    const double *center_;  //!< Atomic coordinates of this center
    int start_;  //!< First basis function in this shell

    int ncartesian_;  //!< How many cartesian functions? (1=s, 3=p, 6=d, ...)
    int nprimitive_;  //!< The number of primitives in this shell
    int nfunction_; //!< Number of basis functions (dependent on puream_, 1=s, 3=p, 5/6=d, ...)

    
    static const char *amtypes;  //!< Lookup array for lowercase letter symbolizing the angular momentup (0=s, 1=p, etc)
    static const char *AMTYPES;  //!< Lookup array for uppercase letter symbolizing the angular momentup (0=S, 1=P, etc)

public:
    // This class is ok to be copy constructed or assigned
    //GaussianShell(const GaussianShell & other) = delete;
    //GaussianShell & operator=(const GaussianShell & other) = delete;
    

    /*!
     * \brief Builds and empty GShell
     */
    GaussianShell() {}


    /*!
     *  \brief Constructor
     * 
     *  \param am Angular momentum.
     *  \param nprimitive Number of primitives in the shell
     *  \param oc An array of contraction coefficients (length \p nprimitive)
     *  \param c An array of normalized contraction coefficients (length \p nprimitive)
     *  \param e An array of exponent values (length \p nprimitive)
     *  \param pure an enum describing whether this shell uses pure or Cartesian functions.
     *  \param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
     *  \param center The x, y, z position of the shell.
     *  \param start The starting index of the first function this shell provides. Used to provide starting positions in matrices.
     */
    GaussianShell(int am, int nprimitive,
                  const double *oc, const double *c,
                  const double *e,
                  ShellInfo::GaussianType pure,
                  int nc, const double* center,
                  int start);


    /*!
     * \brief The number of primitive Gaussians
     */
    int nprimitive() const;


    /*!
     * \brief Total number of basis functions
     */
    int nfunction() const;


    /*!
     * \brief Total number of functions if this shell was Cartesian
     */
    int ncartesian() const          { return ncartesian_; }


    /*!
     * \brief The angular momentum of the given contraction
     */
    int am() const                  { return l_; }


    /*!
     * \brief The character symbol for the angular momentum of the given contraction
     */
    char amchar() const             { return amtypes[l_]; }


    /*!
     * \brief The character symbol for the angular momentum of the given contraction (upper case)
     */
    char AMCHAR() const             { return AMTYPES[l_]; }


    /*!
     * \brief Returns true if contraction is Cartesian
     */
    bool is_cartesian() const       { return !puream_; }


    /*!
     * \brief Returns true if contraction is pure
     */
    bool is_pure() const            { return puream_; }


    /*!
     * \brief Returns the center of the Molecule this shell is on
     */
    const double* center() const;


    /*!
     * \brief Returns the atom number this shell is on. Used by integral derivatives for indexing.
     */
    int ncenter() const             { return nc_; }


    /*!
     * \brief Returns the exponent of the given primitive
     *
     * \param [in] prim Number of the primitive (zero-based)
     * \return Exponent on primitive \p prim
     */
    double exp(int prim) const      { return exp_[prim]; }


    /*!
     * \brief Return coefficient of a primitive
     *
     * \param [in] prim Number of the primitive (zero-based)
     * \return Coefficient on primitive \p prim
     */
    double coef(int prim) const       { return coef_[prim]; }


    /*!
     * \brief Return unnormalized coefficient of a primitive
     *
     * \note Depending on the how the corresponding ShellInfo
     *       was created, this may still be the normalized coefficient
     *       (ie if the basis function is normalized already, this code does
     *       not 'unnormalize' it. 
     * 
     * \param [in] prim Number of the primitive (zero-based)
     * \return Unnormalized coefficient on primitive \p prim
     */
    double original_coef(int prim) const { return original_coef_[prim]; }


    /*!
     * \brief Returns all exponents in this shell
     *
     * \return Exponents for all primitives in this shell (of length nprimitive())
     */
    const double* exps() const { return exp_; }


    /*!
     * \brief Return all coefficients in this shell
     *
     * \return Exponents for all coefficients in this shell (of length nprimitive())
     */
    const double* coefs() const { return coef_; }


    /*!
     * \brief Return all unnormalized coefficients in this shell
     *
     * \note Depending on the how the corresponding ShellInfo
     *       was created, this may still be the normalized coefficient
     *       (ie if the basis function is normalized already, this code does
     *       not 'unnormalize' it. 
     *
     * \return Unnormalized coefficients for all primitives in this shell (of length nprimitive())
     */
    const double* original_coefs() const { return original_coef_; }



    /*!
     * \brief Basis function index where this shell starts.
     * \return Basis function index where this shell starts.
     */
    int function_index() const      { return start_; }



    /*!
     * \brief Print information about this shell
     *
     * Destination of the printing is controlled by SetOutput(). See Output.h
     */
    void print(void) const;
};

}

#endif //PANACHE_GAUSSIANSHELL_H
