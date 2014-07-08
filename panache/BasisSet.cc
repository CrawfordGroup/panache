/*! \file
 * \brief A class to hold basis set information (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <cstdlib>
#include <cmath>

#include "panache/BasisFunctionMacros.h"
#include "panache/BasisSet.h"
#include "panache/Molecule.h"
#include "panache/Output.h"
#include "panache/BasisSetParser.h"

using namespace std;

namespace panache
{

// Constructs a zero AO basis set
BasisSet::BasisSet()
{
    iszero_ = true;

#ifdef PANACHE_USE_LIBERD
    // coef = sqrt( pi^0.75 * (2l-1)!! / (2^(2l+1.5))
    //      = sqrt( pi^0.75 / 2^1.5 )
    //      = (pi /2)^0.75
    const double coef = pow(0.5 * M_PI, 0.75);
#else
    const double coef = 1.0;
#endif

    // Add a dummy atom at the origin, to hold this basis function
    molecule_ = SharedMolecule(new Molecule);
    molecule_->add_atom(0.0, 0.0, 0.0, "");
    // Fill with data representing a single S function, at the origin, with 0 exponent
    n_uprimitive_ = 1;
    n_shells_ = 1;
    nprimitive_ = 1;
    nao_ = 1;
    nbf_ = 1;
    n_prim_per_shell_ = new int[1];
    uexponents_ = new double[1];
    ucoefficients_ = new double[1];
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
    ucoefficients_[0] = coef;
    uoriginal_coefficients_[0] = coef;
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
    shells_[0] = GaussianShell(0, nprimitive_, uoriginal_coefficients_, ucoefficients_, 
                               uexponents_, ShellInfo::GaussianType(0), 0, xyz_, 0);
}


BasisSet::BasisSet(SharedMolecule mol, const std::vector<std::vector<ShellInfo>> & shellmap)
    : molecule_(mol)
{
    construct_(shellmap);
}



BasisSet::BasisSet(const std::shared_ptr<BasisSetParser>& parser,
        const SharedMolecule& mol,
        const std::string& path) : molecule_(mol)
{
    std::vector<std::vector<ShellInfo>> basis_atom_shell;

    for(int i = 0; i < mol->natom(); i++)
    {
        std::vector<std::string> filecontents;

        filecontents = parser->load_file(path);
        basis_atom_shell.push_back(parser->parse(mol->symbol(i), filecontents));
    }

    construct_(basis_atom_shell);
}



BasisSet::~BasisSet()
{
    delete [] n_prim_per_shell_;
    delete [] uexponents_;
    delete [] ucoefficients_;
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



void BasisSet::print_summary(void) const
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
    int *nshells = new int[max_am_ + 1];
    char *amtypes = new char[max_am_ + 1];

    for (int A = 0; A < molecule_->natom(); A++)
    {

        std::fill(nprims, nprims + (max_am_ + 1), 0);
        std::fill(nshells, nshells + (max_am_ + 1), 0);

        output::printf("    %4d    ", A+1);
        output::printf("%2s     ", molecule_->symbol(A).c_str());

        int first_shell = center_to_shell_[A];
        int n_shell = center_to_nshell_[A];

        for (int Q = 0; Q < n_shell; Q++)
        {
            const GaussianShell& shell = shells_[Q + first_shell];
            nshells[shell.am()]++;
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
    delete[] nshells;
    delete[] amtypes;
}

void BasisSet::print_detail(void) const
{
    print_summary();

    output::printf("  ==> AO Basis Functions <==\n");
    output::printf("\n");
    if (has_puream())
        output::printf("    spherical\n");
    else
        output::printf("    cartesian\n");
    output::printf("    ****\n");

    for (int A = 0; A < molecule_->natom(); A++)
    {
        output::printf("   %2s %3d\n",molecule_->symbol(A).c_str(),A+1);

        int first_shell = center_to_shell_[A];
        int n_shell = center_to_nshell_[A];

        for (int Q = 0; Q < n_shell; Q++)
            shells_[Q + first_shell].print();

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



void BasisSet::construct_(const std::vector<std::vector<ShellInfo>> & shellmap)
{

    iszero_ = false;

    int natom = molecule_->natom();

    /// These will tell us where the primitives are for a given center
    std::vector<int>  primitive_start, primitive_end;

    /*
     * First, loop over the unique primitives, and store them
     */
    std::vector<double> uexps;
    std::vector<double> ucoefs;
    std::vector<double> uoriginal_coefs;
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
                n_uprimitive_++;
            }

        }

        primitive_end.push_back(n_uprimitive_);
    }


    if(static_cast<long int>(primitive_end.size()) != natom)
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
    n_prim_per_shell_ = new int[n_shells_];
    // The unique primitives
    uexponents_ = new double[n_uprimitive_];
    ucoefficients_ = new double[n_uprimitive_];
    uoriginal_coefficients_ = new double[n_uprimitive_];
    for(int i = 0; i < n_uprimitive_; ++i)
    {
        uexponents_[i] = uexps[i];
        ucoefficients_[i] = ucoefs[i];
        uoriginal_coefficients_[i] = uoriginal_coefs[i];
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
                                                 &ucoefficients_[ustart+atom_nprim], &uexponents_[ustart+atom_nprim], puream, n, xyz_ptr, bf_count);
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




} // close namespace panache



