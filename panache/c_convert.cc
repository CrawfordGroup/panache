/*! \file
 *  \brief Conversions from C-style arrays to C++ objects (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/c_convert.h"

namespace panache
{

SharedBasisSet BasisSetFromArrays(SharedMolecule molecule,
        int_t ncenters,
        int_t * nshellspercenter,
        struct C_ShellInfo * shells)
{
    // Construct the basis set info
    std::vector<std::vector<panache::ShellInfo>> shellmap;

    int counter = 0;

    for(int i = 0; i < ncenters; ++i)
    {
        std::vector<panache::ShellInfo> shellvec;

        panache::Vector3 cen = molecule->xyz(i);

        for(int j = 0; j < nshellspercenter[i]; j++)
        {
            std::vector<double> c,e;
            for(int k = 0; k < shells[counter].nprim; k++)
            {
                c.push_back(shells[counter].coef[k]);
                e.push_back(shells[counter].exp[k]);
            }
            shellvec.push_back(panache::ShellInfo(shells[counter].am, c, e,
                                                  shells[counter].ispure ? panache::ShellInfo::GaussianType::Pure : panache::ShellInfo::GaussianType::Cartesian,
                                                  i, cen)); 
            counter++;
        }

        shellmap.push_back(shellvec);

    }

    return SharedBasisSet(new panache::BasisSet(molecule, shellmap));

}

SharedMolecule MoleculeFromArrays(int_t ncenters, C_AtomCenter * atoms)
{
    SharedMolecule molecule(new panache::Molecule());

    for(auto i = 0; i < ncenters; i++)
    {
        molecule->add_atom(
            atoms[i].center[0],
            atoms[i].center[1],
            atoms[i].center[2],
            atoms[i].symbol);
    }

    return molecule;
}


} // close namespace panache

