/*! \file
 * \brief Holds information about a basis set shell (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_SHELLINFO_H
#define PANACHE_SHELLINFO_H

#include <vector>
#include "Vector3.h"

namespace panache {


/*!
 *  \brief Holds information about a Gaussian orbital shell.
 *
 *  This class has the same behavior as GaussianShell, but implements everything using
 *  slower data structures, which are easier to construct. These are used to build the
 *  basis set, which builds more efficient pointer-based GaussianShell objects.
 */
class ShellInfo
{
private:
    
    int l_;  //!< Angular momentum
    int puream_;  //!< Flag for pure angular momentum
    std::vector<double> exp_;  //!< Exponents (of length nprimitives_)
    std::vector<double> coef_;  //!< Contraction coefficients (of length nprimitives_)
    std::vector<double> original_coef_;  //!< Original (un-normalized) contraction coefficients (of length nprimitives)

    int nc_;  //!< Atom number this shell goes to. Needed when indexing integral derivatives.
    Vector3 center_;  //!< Atomic center number in the Molecule
    int ncartesian_;  //!< How many cartesian functions? (1=s, 3=p, 6=d, ...)
    int nfunction_; //!< Number of basis functions (dependent on puream_, 1=s, 3=p, 5/6=d, ...)

    /*!
     *  \brief Normalizes a single primitive.
     *
     *  \param p The primitive index to normalize.
     *  \return Normalization constant to be applied to the primitive.
     */
    double primitive_normalization(int p);


    /*!
     *  \brief Normalizes the entire contraction set. Applies the normalization to the coefficients
     */
    void contraction_normalization();


    static const char *amtypes;  //!< Lookup array for lowercase letter symbolizing the angular momentup (0=s, 1=p, etc)
    static const char *AMTYPES;  //!< Lookup array for uppercase letter symbolizing the angular momentup (0=S, 1=P, etc)


public:
    // This class is ok to be copy constructed or assigned
    //ShellInfo(const ShellInfo & other) = delete;
    //ShellInfo & operator=(const ShellInfo & other) = delete;

   
    /*!
     * \brief Are the primitives normalized or not
     */
    enum PrimitiveType {
        Normalized, //!< Is normalized
        Unnormalized //!< Is NOT normalized
    };

    
    /*!
     * \brief Are the primitives pure spherical or cartesian basis funtions
     */
    enum GaussianType {
        Cartesian = 0, //!< Cartesian (6d, 10f, etc)
        Pure = 1  //!< Spherical (5d, 7f, etc)
    };


    /*!
     *  Constructor.
     *
     *  \param am Angular momentum.
     *  \param c An array of contraction coefficients.
     *  \param e An array of exponent values.
     *  \param pure Pure spherical harmonics, or Cartesian.
     *  \param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
     *  \param center The x, y, z position of the shell. This is passed to reduce the number of calls to the molecule.
     *  \param pt Is the shell already normalized?
     */
    ShellInfo(int am,
                  const std::vector<double>& c,
                  const std::vector<double>& e,
                  GaussianType pure,
                  int nc,
                  const Vector3& center,
                  PrimitiveType pt = Normalized);

    /*!
     * \brief Normalizes the entire shell
     */
    void normalize_shell();

    /*!
     * \brief Make a copy of the ShellInfo.
     */
    ShellInfo copy();

    /*!
     * \brief Make a copy of the ShellInfo.
     */
    ShellInfo copy(int nc, const Vector3& c);


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
    const Vector3& center() const;



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
     * \note Depending on the how this object
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
    const std::vector<double>& exps() const { return exp_; }


    /*!
     * \brief Return all coefficients in this shell
     *
     * \return Exponents for all coefficients in this shell (of length nprimitive())
     */
    const std::vector<double>& coefs() const { return coef_; }



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
    const std::vector<double>& original_coefs() const { return original_coef_; }



    /*!
     * \brief Print information about this shell
     *
     * Destination of the printing is controlled by SetOutput(). See Output.h
     */
    void print(void) const;
};

} // end namespace panache

#endif //PANACHE_SHELLINFO_H
