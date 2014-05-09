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

/*!
    \defgroup MINTS libmints: Integral library
    \ingroup MINTS
*/
#include <cstdlib>

#ifdef USE_LIBINT
#include <libint/libint.h>
#define MAX_AM LIBINT_MAX_AM
#endif

#ifdef USE_LIBINT2
#include <libint2.h>
#define MAX_AM LIBINT2_MAX_AM_ERI
#endif

#include "BasisFunctionMacros.h"
#include "BasisSet.h"
#include "Molecule.h"
#include "libciomr.h"
#include "Output.h"

using namespace std;

namespace panache
{

bool BasisSet::initialized_shared_ = false;

std::vector<Vector3> BasisSet::exp_ao[MAX_AM];

// Constructs a zero AO basis set
BasisSet::BasisSet()
{
    if (initialized_shared_ == false)
        initialize_singletons();
    initialized_shared_ = true;

    // Add a dummy atom at the origin, to hold this basis function
    molecule_ = shared_ptr<Molecule>(new Molecule);
    molecule_->add_atom(0, 0.0, 0.0, 0.0);
    // Fill with data representing a single S function, at the origin, with 0 exponent
    n_uprimitive_ = 1;
    n_shells_ = 1;
    nprimitive_ = 1;
    nao_ = 1;
    nbf_ = 1;
    n_prim_per_shell_ = new int[1];
    uexponents_ = new double[1];
    ucoefficients_ = new double[1];
    uerd_coefficients_ = new double[1];
    uoriginal_coefficients_ = new double[1];
    shell_first_ao_ = new int[1];
    shell_first_basis_function_ = new int[1];
    shells_ = new GaussianShell[1];
    ao_to_shell_ = new int[1];
    function_to_shell_ = new int[1];
    function_center_ = new int[1];
    shell_center_ = new int[1];
    center_to_nshell_ = new int[1];
    center_to_shell_ = new int[1];
    xyz_ = new double[3];
    n_prim_per_shell_[0] = 1;
    uexponents_[0] = 0.0;
    ucoefficients_[0] = 1.0;
    uoriginal_coefficients_[0] = 1.0;
    shell_first_ao_[0] = 0;
    shell_first_basis_function_[0] = 0;
    ao_to_shell_[0] = 0;
    function_to_shell_[0] = 0;
    function_center_[0] = 0;
    shell_center_[0] = 0;
    center_to_nshell_[0] = 1;
    center_to_shell_[0] = 0;
    puream_ = 0;
    max_am_ = 0;
    max_nprimitive_ = 1;
    xyz_[0] = 0.0;
    xyz_[1] = 0.0;
    xyz_[2] = 0.0;
    shells_[0] = GaussianShell(0, nprimitive_, uoriginal_coefficients_, ucoefficients_, uerd_coefficients_,
                               uexponents_, ShellInfo::GaussianType(0), 0, xyz_, 0);
}

BasisSet::~BasisSet()
{
    delete_arrays();
}

void BasisSet::delete_arrays(void)
{
    delete [] n_prim_per_shell_;
    delete [] uexponents_;
    delete [] ucoefficients_;
    delete [] uerd_coefficients_;
    delete [] uoriginal_coefficients_;
    delete [] shell_first_ao_;
    delete [] shell_first_basis_function_;
    delete [] shells_;
    delete [] ao_to_shell_;
    delete [] function_to_shell_;
    delete [] function_center_;
    delete [] shell_center_;
    delete [] center_to_nshell_;
    delete [] center_to_shell_;
    delete [] xyz_;
}


void BasisSet::initialize_singletons()
{
    // Populate the exp_ao arrays
    for (int l=0; l<MAX_AM; ++l)
    {
        for (int i=0; i<=l; ++i)
        {
            int x = l-i;
            for (int j=0; j<=i; ++j)
            {
                int y = i-j;
                int z = j;

                Vector3 xyz_ao(x, y, z);
                BasisSet::exp_ao[l].push_back(xyz_ao);
            }
        }
    }
}

shared_ptr<Molecule> BasisSet::molecule() const
{
    return molecule_;
}

void BasisSet::print(FILE * out) const
{
    output::printf("  Basis Set\n");
    output::printf("    Number of shells: %d\n", nshell());
    output::printf("    Number of basis function: %d\n", nbf());
    output::printf("    Number of Cartesian functions: %d\n", nao());
    output::printf("    Spherical Harmonics?: %s\n", has_puream() ? "true" : "false");
    output::printf("    Max angular momentum: %d\n\n", max_am());
}

void BasisSet::print_summary(FILE* out) const
{
    output::printf("  -AO BASIS SET INFORMATION:\n");
    output::printf("    Total number of shells = %d\n", nshell());
    output::printf("    Number of primitives   = %d\n", nprimitive_);
    output::printf("    Number of AO           = %d\n", nao_);
    output::printf("    Number of SO           = %d\n", nbf_);
    output::printf("    Maximum AM             = %d\n", max_am_);
    output::printf("    Spherical Harmonics    = %s\n", (puream_ ? "TRUE" : "FALSE"));
    output::printf("\n");

    output::printf("  -Contraction Scheme:\n");
    output::printf("    Atom   Type   All Primitives // Shells:\n");
    output::printf("   ------ ------ --------------------------\n");

    int *nprims = new int[max_am_ + 1];
    int *nunique = new int[max_am_ + 1];
    int *nshells = new int[max_am_ + 1];
    char *amtypes = new char[max_am_ + 1];

    for (int A = 0; A < molecule_->natom(); A++)
    {

        memset((void*) nprims , '\0', (max_am_ + 1) * sizeof(int));
        memset((void*) nunique, '\0', (max_am_ + 1) * sizeof(int));
        memset((void*) nshells, '\0', (max_am_ + 1) * sizeof(int));

        output::printf("    %4d    ", A+1);
        output::printf("%2s     ", molecule_->symbol(A).c_str());

        int first_shell = center_to_shell_[A];
        int n_shell = center_to_nshell_[A];

        for (int Q = 0; Q < n_shell; Q++)
        {
            const GaussianShell& shell = shells_[Q + first_shell];
            nshells[shell.am()]++;
            nunique[shell.am()]+= shell.nprimitive();
            nprims [shell.am()]+= shell.nprimitive();
            amtypes[shell.am()] = shell.amchar();
        }

        // All Primitives
        for (int l = 0; l < max_am_ + 1; l++)
        {
            if (nprims[l] == 0)
                continue;
            output::printf("%d%c ", nprims[l], amtypes[l]);
        }
        // Shells
        output::printf("// ");
        for (int l = 0; l < max_am_ + 1; l++)
        {
            if (nshells[l] == 0)
                continue;
            output::printf("%d%c ", nshells[l], amtypes[l]);
        }
        output::printf("\n");
    }
    output::printf("\n");

    delete[] nprims;
    delete[] nunique;
    delete[] nshells;
    delete[] amtypes;
}

void BasisSet::print_detail(FILE * out) const
{
    print_summary(out);

    output::printf("  ==> AO Basis Functions <==\n");
    output::printf("\n");
    if (has_puream())
        output::printf("    spherical\n");
    else
        output::printf("    cartesian\n");
    output::printf("    ****\n");

    for (int uA = 0; uA < molecule_->nunique(); uA++)
    {
        const int A = molecule_->unique(uA);
        output::printf("   %2s %3d\n",molecule_->symbol(A).c_str(),A+1);

        int first_shell = center_to_shell_[A];
        int n_shell = center_to_nshell_[A];

        for (int Q = 0; Q < n_shell; Q++)
            shells_[Q + first_shell].print(out);

        output::printf("    ****\n");
    }
    output::printf("\n");
}

const GaussianShell& BasisSet::shell(int si) const
{
    if (si < 0 || si > nshell())
        throw RuntimeError("BasisSet::shell: requested shell is out-of-bounds.");
    return shells_[si];
}

const GaussianShell& BasisSet::shell(int center, int si) const
{
    return shell(center_to_shell_[center] + si);
}

shared_ptr<BasisSet> BasisSet::zero_ao_basis_set()
{
    // In the new implementation, we simply call the default constructor
    shared_ptr<BasisSet> new_basis(new BasisSet());
    return new_basis;
}

//shared_ptr<SOBasisSet> BasisSet::zero_so_basis_set(const shared_ptr<IntegralFactory>& factory)
//{
//    shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
//    shared_ptr<SOBasisSet> sozero(new SOBasisSet(zero, factory));
//    return sozero;
//}


BasisSet::BasisSet(SharedMolecule mol, const std::vector<std::vector<ShellInfo>> & shellmap)
    : molecule_(mol)
{
    // Singletons
    if (initialized_shared_ == false)
        initialize_singletons();
    initialized_shared_ = true;

    int natom = molecule_->natom();

    /// These will tell us where the primitives are for a given center
    std::vector<int>  primitive_start, primitive_end;

    /*
     * First, loop over the unique primitives, and store them
     */
    std::vector<double> uexps;
    std::vector<double> ucoefs;
    std::vector<double> uoriginal_coefs;
    std::vector<double> uerd_coefs;
    n_uprimitive_ = 0;
    for(auto & centerit : shellmap)
    {
        primitive_start.push_back(n_uprimitive_);

        for(auto & shell : centerit)
        {
            for(int prim = 0; prim < shell.nprimitive(); ++prim)
            {
                uexps.push_back(shell.exp(prim));
                ucoefs.push_back(shell.coef(prim));
                uoriginal_coefs.push_back(shell.original_coef(prim));
                uerd_coefs.push_back(shell.erd_coef(prim));
                n_uprimitive_++;
            }

        }

        primitive_end.push_back(n_uprimitive_);
    }


    if(primitive_end.size() != natom)
        throw RuntimeError("SIZE MISMATCH");

    /*
     * Count basis functions, shells and primitives
     */
    n_uprimitive_ = 0;
    n_shells_ = 0;
    nprimitive_ = 0;
    nao_ = 0;
    nbf_ = 0;

    for (auto & it : shellmap)
    {
        for(auto & shell : it)
        {
            int nprim = shell.nprimitive();
            n_uprimitive_ += nprim;
            nprimitive_ += nprim;
            n_shells_++;
            nao_ += shell.ncartesian();
            nbf_ += shell.nfunction();
        }
    }


    /*
     * Allocate arrays
     */
    delete_arrays();
    n_prim_per_shell_ = new int[n_shells_];
    // The unique primitives
    uexponents_ = new double[n_uprimitive_];
    ucoefficients_ = new double[n_uprimitive_];
    uoriginal_coefficients_ = new double[n_uprimitive_];
    uerd_coefficients_ = new double[n_uprimitive_];
    for(int i = 0; i < n_uprimitive_; ++i)
    {
        uexponents_[i] = uexps[i];
        ucoefficients_[i] = ucoefs[i];
        uoriginal_coefficients_[i] = uoriginal_coefs[i];
        uerd_coefficients_[i] = uerd_coefs[i];
    }

    shell_first_ao_ = new int[n_shells_];
    shell_first_basis_function_ = new int[n_shells_];
    shells_ = new GaussianShell[n_shells_];
    ao_to_shell_ = new int[nao_];
    function_to_shell_ = new int[nbf_];
    function_center_ = new int[nbf_];
    shell_center_ = new int[n_shells_];
    center_to_nshell_ = new int[natom];
    center_to_shell_ = new int[natom];
    xyz_ = new double[3*natom];

    /*
     * Now loop over all atoms, and point to the appropriate unique data
     */
    int shell_count = 0;
    int ao_count = 0;
    int bf_count = 0;
    double *xyz_ptr = xyz_;
    puream_ = false;
    max_am_ = 0;
    max_nprimitive_ = 0;
    for (int n = 0; n < natom; ++n)
    {
        const shared_ptr<CoordEntry> &atom = molecule_->atom_entry(n);
        const vector<ShellInfo>& shells = shellmap[n];
        int ustart = primitive_start[n];
        int uend = primitive_end[n];
        int nshells = shells.size();
        center_to_nshell_[n] = nshells;
        center_to_shell_[n] = shell_count;
        int atom_nprim = 0;
        for (int i = 0; i < nshells; ++i)
        {
            const ShellInfo &thisshell = shells[i];
            shell_first_ao_[shell_count] = ao_count;
            shell_first_basis_function_[shell_count] = bf_count;
            int shell_nprim = thisshell.nprimitive();
            int am = thisshell.am();
            max_nprimitive_ = shell_nprim > max_nprimitive_ ? shell_nprim : max_nprimitive_;
            max_am_ = max_am_ > am ? max_am_ : am;
            shell_center_[shell_count] = n;
            ShellInfo::GaussianType puream = thisshell.is_pure() ? ShellInfo::Pure : ShellInfo::Cartesian;
            if(puream)
                puream_ = true;
//            output::printf("atom %d basis %s shell %d nprim %d atom_nprim %d\n", n, basis.c_str(), i, shell_nprim, atom_nprim);
            shells_[shell_count] = GaussianShell(am, shell_nprim, &uoriginal_coefficients_[ustart+atom_nprim],
                                                 &ucoefficients_[ustart+atom_nprim], &uerd_coefficients_[ustart+atom_nprim], &uexponents_[ustart+atom_nprim], puream, n, xyz_ptr, bf_count);
            for(int thisbf = 0; thisbf < thisshell.nfunction(); ++thisbf)
            {
                function_to_shell_[bf_count] = shell_count;
                function_center_[bf_count++] = n;
            }
            for(int thisao = 0; thisao < thisshell.ncartesian(); ++thisao)
            {
                ao_to_shell_[ao_count++] = shell_count;
            }
            atom_nprim += shell_nprim;
            shell_count++;
        }
        Vector3 xyz = molecule_->xyz(n);
        xyz_ptr[0] = xyz[0];
        xyz_ptr[1] = xyz[1];
        xyz_ptr[2] = xyz[2];
        xyz_ptr += 3;
        if(atom_nprim != uend-ustart)
        {
            throw RuntimeError("Problem with nprimitive in basis set construction!");
        }
    }

}


void BasisSet::compute_phi(double *phi_ao, double x, double y, double z)
{
    zero_arr(phi_ao, nao());

    int ao = 0;
    for(int ns=0; ns < nshell(); ns++)
    {
        const GaussianShell& shell = shells_[ns];
        int am = shell.am();
        int nprim = shell.nprimitive();
        const double *a = shell.exps();
        const double *c = shell.coefs();

        const double *xyz = shell.center();
        double dx = x - xyz[0];
        double dy = y - xyz[1];
        double dz = z - xyz[2];
        double rr = dx*dx + dy*dy + dz*dz;

        double cexpr = 0;
        for(int np=0; np < nprim; np++)
            cexpr += c[np] * exp(-a[np] * rr);

        for(int l=0; l < INT_NCART(am); l++)
        {
            Vector3& components = exp_ao[am][l];
            phi_ao[ao+l] += pow(dx, (double) components[0]) *
                            pow(dy, (double) components[1]) *
                            pow(dz, (double) components[2]) *
                            cexpr;
        }

        ao += INT_NCART(am);
    } // nshell
}

// Free functions
int **compute_shell_map(int **atom_map, const std::shared_ptr<BasisSet> &basis)
{
    int **shell_map;

    BasisSet& gbs = *basis.get();
    Molecule& mol = *gbs.molecule().get();

    // create the character table for the point group
    CharacterTable ct = mol.point_group()->char_table();

    int natom = mol.natom();
    int ng = ct.order();

    int nshell = basis->nshell();
    shell_map = new int*[nshell];
    for (int i=0; i < nshell; i++)
        shell_map[i] = new int[ng];

    for (int i=0; i<natom; i++) {
        // hopefully, shells on equivalent centers will be numbered in the same
        // order
        for (int s=0; s < gbs.nshell_on_center(i); s++) {
            int shellnum = gbs.shell_on_center(i,s);
            for (int g=0; g < ng; g++) {
                shell_map[shellnum][g] = gbs.shell_on_center(atom_map[i][g],s);
            }
        }
    }

    return shell_map;
}

void delete_shell_map(int **shell_map, const std::shared_ptr<BasisSet> &basis)
{
    int nshell = basis->nshell();
    if (shell_map) {
        for (int i=0; i < nshell; i++)
            delete[] shell_map[i];
        delete[] shell_map;
    }
}


} // close namespace panache



