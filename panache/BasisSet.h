/*! \file
 * \brief A class to hold basis set information (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_BASISSET_H
#define PANACHE_BASISSET_H

#include <map>
#include <memory>
#include "panache/GaussianShell.h"

namespace panache {

class Molecule;
typedef std::shared_ptr<Molecule> SharedMolecule;

class BasisSetParser;
typedef std::shared_ptr<BasisSetParser> SharedBasisSetParser;


/*! 
 * \brief Basis set container class
 *
 * Can be constructed by passing it raw ShellInfo objects
 * or by specifying a molecule, file, and parser and a file
 *
*/
class BasisSet
{
private:
    bool iszero_;  //!< True if this is a zero basis set


    GaussianShell *shells_;  //!< Holds most of the shell information

    //! Molecule object.
    SharedMolecule molecule_;

    int nao_;            //!< Number of atomic orbitals (Cartesian)
    int nbf_;            //!< Number of basis functions (either cartesian or spherical)
    int n_uprimitive_;   //!< The number of unique primitives
    int n_shells_;       //!< Number of shells
    int nprimitive_;     //!< Total number of primitives
    int max_am_;         //!< Maximum angular momentum
    int max_nprimitive_; //!< Maximum number of primitives in all shells
    bool puream_;        //!< True if the basis set uses spherical basis functions

    
    int *n_prim_per_shell_;     //!< The number of primitives (and exponents) in each shell
    int *shell_first_ao_;       //!< The first (Cartesian) atomic orbital in each shell
    int *shell_first_basis_function_; //!< The first (Cartesian / spherical) basis function in each shell
    int *shell_center_;         //!< Shell number to atomic center.
    int *function_to_shell_;    //!< Which shell does a given (Cartesian / spherical) function belong to?
    int *ao_to_shell_;          //!< Which shell does a given Cartesian function belong to?
    int *function_center_;      //!< Which center is a given function on?
    int *center_to_nshell_;     //!< How many shells are there on each center?
    int *center_to_shell_;      //!< What's the first shell on each center?

    
    double *uexponents_;                //!< The flattened lists of unique exponents
    double *ucoefficients_;             //!< The flattened lists of unique contraction coefficients (normalized)
    double *uoriginal_coefficients_;    //!< The flattened lists of unique contraction coefficients (as provided by the user)
    double *xyz_;                       //!< The flattened list of Cartesian coordinates for each atom


    /*!
     * \brief Constructs the basis set object from ShellInfo objects
     *
     * \warning molecule_ is expected to be set already!
     *
     * \param [in] shellmap A map of centers to vector of shells. That is,
     *                      shellmap[0] contains a vector of all shells for the first center,
     *                      shellmap[1] contains a vector of all shells for the second center, etc.
     */ 
    void construct_(const std::vector<std::vector<ShellInfo>> & shellmap);


public:
    /*!
     * \brief Constructs an empty (zero) basis set
     *
     *  A zero BasisSet object that actually has a single s-function
     *  at the origin with an exponent of 0.0 and contraction of 1.0.
     */
    BasisSet();


    /*!
     * Frees all memory
     */ 
    ~BasisSet();



    /*!
     * \brief Constructs a BasisSet with the given shells and molecule
     *
     * \param [in] mol A molecule (information about the centers for this basis set)
     * \param [in] shellmap A map of centers to vector of shells. That is,
     *                      shellmap[0] contains a vector of all shells for the first center,
     *                      shellmap[1] contains a vector of all shells for the second center, etc.
     */
    BasisSet(SharedMolecule mol, const std::vector<std::vector<ShellInfo>> & shellmap);



    /*!
     * \brief Constructs a BasisSet from a file
     *
     * \param [in] parser An object that is used to parse the basis set file
     * \param [in] mol A molecule (information about the centers for this basis set)
     * \param [in] path Full path to a file containing the basis set in a format that \p parser can read. 
     */
    BasisSet(const SharedBasisSetParser & parser,
             const SharedMolecule & mol,
             const std::string & path);



    /*!
     *  Is this a zero basis set
     *
     *  \return True if this is a zero basis set, false if it is not
     */
    bool iszero(void) const { return iszero_; } 




    /*!
     *  \brief Number of primitives.
     *
     *  \return The total number of primitives in all contractions.
     */
    int nprimitive() const             { return nprimitive_; }


    /*!
     *  \brief Maximum number of primitives in a shell.
     *
     *  Examines each shell and find the shell with the maximum number of primitives returns that
     *  number of primitives.
     *
     *  \return Maximum number of primitives.
     */
    int max_nprimitive() const         { return max_nprimitive_; }


    /*!
     * \brief Number of shells.
     *
     * \return Number of shells in the basis set
     */
    int nshell() const                 { return n_shells_;  }


    /*!
     * \brief Number of atomic orbitals (Cartesian).
     *
     * \return The number of atomic orbitals (Cartesian orbitals, always).
     */
    int nao() const                    { return nao_;         }



    /*!
     *  \brief Number of basis functions (Spherical or cartesian)
     *
     *  \return The number of basis functions (Spherical, if has_puream() == true).
     */
    int nbf() const                    { return nbf_;         }



    /*!
     *  \brief Maximum angular momentum used in the basis set.
     *
     *  \return Maximum angular momentum.
     */
    int max_am() const                 { return max_am_;      }



    /*!
     *  \brief Spherical harmonics?
     *
     *  \return true if using spherical harmonics
     */
    bool has_puream() const            { return puream_;      }



    /*!
     *  \brief Compute the maximum number of basis functions contained in a shell.
     *
     *  \return The max number of basis functions in a shell.
     */
    int max_function_per_shell() const { return (puream_) ? 2*max_am_+1 : (max_am_+1)*(max_am_+2)/2; }



    /*!
     * \brief Molecule this basis is for.
     *
     *  \return Shared pointer to the molecule for this basis set.
     */
    SharedMolecule molecule() const { return molecule_; }



    /*!
     *  \brief Given a shell what is its first AO function
     *
     *  \param [in] i Shell number
     *  \return The function number for the first function for the i'th shell.
     */
    int shell_to_ao_function(int i) const { return shell_first_ao_[i]; }


    /*!
     *  \brief Given a shell what is its atomic center
     *
     *  \param [in] i Shell number
     *  \return The atomic center for the i'th shell.
     */
    int shell_to_center(int i) const { return shell_center_[i]; }



    /*!
     *  \brief Given a shell what is its first basis function (spherical) function
     *
     *  \param [in] i Shell number
     *  \return The function number for the first function for the i'th shell.
     */
    int shell_to_basis_function(int i) const { return shell_first_basis_function_[i]; }



    /*!
     *  \brief Shell corresponding to a given basis function
     *
     *  \param [in] i Number of a basis function
     *  \return Number of the shell which contains function \p i
     */
    int function_to_shell(int i) const { return function_to_shell_[i]; }



    /*!
     *  \brief Given a function what is its atomic center
     *
     *  \param [in] i Number of a basis function
     *  \return Number of the atomic center for function \p i.
     */
    int function_to_center(int i) const { return function_center_[i]; }




    /*!
     *  \brief Given a Cartesian function (AO) number what shell does it correspond to.
     *
     *  \param [in] i Number of a cartesian AO
     *  \return Number of the shell which contains cartesian AO function \p i
     */

    int ao_to_shell(int i) const { return ao_to_shell_[i]; }


    /*!
     *  \brief Get a Gaussian shell
     *
     *  \param [in] si Number of the shell
     *  \return A const reference GaussianShell object for shell \p si
     */
    const GaussianShell& shell(int si) const;



    /*!
     *  \brief Get a Gaussian shell on a given center
     *
     *  \param [in] si Shell number for center \p center
     *  \param [in] center The number of the center
     *  \return A const reference to the GaussianShell object for shell \p i on center \p center
     */
    const GaussianShell& shell(int center, int si) const;


    /*!
     * \brief Get the number of shells on a given center.
     *
     * \param [in] i Number of the center
     * \return Number of shells on center \p i
     */
    int nshell_on_center(int i) const { return center_to_nshell_[i]; }


    /*!
     * \brief Get the overall shell number for a shell on a given a center
     *
     * \param [in] center Center that the shell is on
     * \param [in] shell Number of the shell on that center
     * \return Absolute number for that shell
     */
    int shell_on_center(int center, int shell) const { return center_to_shell_[center] + shell; }


    /*!
     * \brief Print a quick summary about this basis set
     *
     * Information is printed through the global Output interface (see Output.h)
     */
    void print_summary(void) const;


    /*!
     * \brief Print details about this basis set
     *
     * Details include all exponents and coefficients, etc.
     *
     * Information is printed through the global Output interface (see Output.h)
     */
    void print_detail(void) const;
};


/*!
 * \brief A shared pointer to a BasisSet object
 */
typedef std::shared_ptr<BasisSet> SharedBasisSet;

} // close namespace panache

#endif //PANACHE_BASISSET_H
