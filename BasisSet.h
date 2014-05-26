/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef PANACHE_BASISSET_H
#define PANACHE_BASISSET_H

#include <map>
#include <memory>
#include "GaussianShell.h"

using std::shared_ptr;

namespace panache {

class Molecule;
typedef shared_ptr<Molecule> SharedMolecule;

class BasisSetParser;

/*! \ingroup MINTS */

//! Basis set container class
/*! Reads the basis set from a checkpoint file object. Also reads the molecule
    from the checkpoint file storing the information in an internal Molecule class
    which can be accessed using molecule().
*/
class BasisSet
{
    //! Array of gaussian shells
    GaussianShell *shells_;

    //! vector of shells numbers sorted in acending AM order.
    std::vector<int> sorted_ao_shell_list_;

    //! Molecule object.
    shared_ptr<Molecule> molecule_;

    // Has static information been initialized?
    static bool initialized_shared_;

    /*
     * Scalars
     */
    /// Number of atomic orbitals (Cartesian)
    int nao_;
    /// Number of basis functions (either cartesian or spherical)
    int nbf_;
    /// The number of unique primitives
    int n_uprimitive_;
    /// The number of shells
    int n_shells_;
    /// The number of primitives
    int nprimitive_;
    /// The maximum angular momentum
    int max_am_;
    /// The maximum number of primitives in a shell
    int max_nprimitive_;
    /// Whether the basis set is uses spherical basis functions or not
    bool puream_;

    /*
     * Arrays
     */
    /// The number of primitives (and exponents) in each shell
    int *n_prim_per_shell_;
    /// The first (Cartesian) atomic orbital in each shell
    int *shell_first_ao_;
    /// The first (Cartesian / spherical) basis function in each shell
    int *shell_first_basis_function_;
    /// Shell number to atomic center.
    int *shell_center_;
    /// Which shell does a given (Cartesian / spherical) function belong to?
    int *function_to_shell_;
    /// Which shell does a given Cartesian function belong to?
    int *ao_to_shell_;
    /// Which center is a given function on?
    int *function_center_;
    /// How many shells are there on each center?
    int *center_to_nshell_;
    /// What's the first shell on each center?
    int *center_to_shell_;

    /// The flattened lists of unique exponents
    double *uexponents_;
    /// The flattened lists of unique contraction coefficients (normalized)
    double *ucoefficients_;
    /// The flattened lists of unique contraction coefficients (as provided by the user)
    double *uoriginal_coefficients_;
    /// The flattened lists of ERD normalized contraction coefficients
    double *uerd_coefficients_;
    /// The flattened list of Cartesian coordinates for each atom
    double *xyz_;

    void delete_arrays(void);


public:
    BasisSet();

    ~BasisSet();

    BasisSet(SharedMolecule mol, const std::vector<std::vector<ShellInfo>> & shellmap);

    /** Initialize singleton values that are shared by all basis set objects. */
    static void initialize_singletons();

    /** Number of primitives.
     *  @return The total number of primitives in all contractions.
     */
    int nprimitive() const             { return nprimitive_; }
    /** Maximum number of primitives in a shell.
     *  Examines each shell and find the shell with the maximum number of primitives returns that
     *  number of primitives.
     *  @return Maximum number of primitives.
     */
    int max_nprimitive() const         { return max_nprimitive_; }
    /** Number of shells.
     *  @return Number of shells.
     */
    int nshell() const                 { return n_shells_;  }
    /** Number of atomic orbitals (Cartesian).
     * @return The number of atomic orbitals (Cartesian orbitals, always).
     */
    int nao() const                    { return nao_;         }
    /** Number of basis functions (Spherical).
     *  @return The number of basis functions (Spherical, if has_puream() == true).
     */
    int nbf() const                    { return nbf_;         }
    /** Maximum angular momentum used in the basis set.
     *  @return Maximum angular momentum.
     */
    int max_am() const                 { return max_am_;      }
    /** Spherical harmonics?
     *  @return true if using spherical harmonics
     */
    bool has_puream() const            { return puream_;      }
    /** Compute the maximum number of basis functions contained in a shell.
     *  @return The max number of basis functions in a shell.
     */
    int max_function_per_shell() const { return (puream_) ? 2*max_am_+1 : (max_am_+1)*(max_am_+2)/2; }
    /** Molecule this basis is for.
     *  @return Shared pointer to the molecule for this basis set.
     */
    shared_ptr<Molecule> molecule() const;
    /** Given a shell what is its first AO function
     *  @param i Shell number
     *  @return The function number for the first function for the i'th shell.
     */
    int shell_to_ao_function(int i) const { return shell_first_ao_[i]; }
    /** Given a shell what is its atomic center
     *  @param i Shell number
     *  @return The atomic center for the i'th shell.
     */
    int shell_to_center(int i) const { return shell_center_[i]; }
    /** Given a shell what is its first basis function (spherical) function
     *  @param i Shell number
     *  @return The function number for the first function for the i'th shell.
     */
    int shell_to_basis_function(int i) const { return shell_first_basis_function_[i]; }

    /** Given a function number what shell does it correspond to. */
    int function_to_shell(int i) const { return function_to_shell_[i]; }
    /** Given a function what is its atomic center
     *  @param i Function number
     *  @return The atomic center for the i'th function.
     */
    int function_to_center(int i) const { return function_center_[i]; }

    /** Given a Cartesian function (AO) number what shell does it correspond to. */
    int ao_to_shell(int i) const { return ao_to_shell_[i]; }

    /** Return the si'th Gaussian shell
     *  @param i Shell number
     *  @return A shared pointer to the GaussianShell object for the i'th shell.
     */
    const GaussianShell& shell(int si) const;

    /** Return the i'th Gaussian shell on center
     *  @param i Shell number
     *  @return A shared pointer to the GaussianShell object for the i'th shell.
     */
    const GaussianShell& shell(int center, int si) const;

    /// Return the number of shells on a given center.
    int nshell_on_center(int i) const { return center_to_nshell_[i]; }
    /// Return the overall shell number
    int shell_on_center(int center, int shell) const { return center_to_shell_[center] + shell; }

    /** Returns an empty basis set object.
     *
     *  Returns a BasisSet object that actually has a single s-function
     *  at the origin with an exponent of 0.0 and contraction of 1.0.
     *  @return A new empty BasisSet object.
     */
    static shared_ptr<BasisSet> zero_ao_basis_set();

    /// Global arrays of x, y, z exponents
    static std::vector<Vector3> exp_ao[];

    //! Returns the value of the sorted shell list.
    int get_ao_sorted_shell(const int &i) { return sorted_ao_shell_list_[i]; }
    //! Returns the vector of sorted shell list.
    std::vector<int> get_ao_sorted_list() { return sorted_ao_shell_list_; }

    // Returns the values of the basis functions at a point
    void compute_phi(double *phi_ao, double x, double y, double z);

    std::shared_ptr<BasisSet> construct(const std::shared_ptr<BasisSetParser>& parser,
                                        const std::shared_ptr<Molecule>& mol,
                                        const std::string& basisname);


    void print(FILE *out = stdout) const;
    void print_summary(FILE *out = stdout) const;
    void print_detail(FILE *out = stdout) const;
};

inline
bool shell_sorter_ncenter(const GaussianShell& d1, const GaussianShell& d2)
{
    return d1.ncenter() < d2.ncenter();
}

inline
bool shell_sorter_am(const GaussianShell& d1, const GaussianShell& d2)
{
    return d1.am() < d2.am();
}


// Free functions
int **compute_shell_map(int **atom_map, const std::shared_ptr<BasisSet> &);
void delete_shell_map(int **shell_map, const std::shared_ptr<BasisSet> &);

} // close namespace panache

#endif //PANACHE_BASISSET_H
