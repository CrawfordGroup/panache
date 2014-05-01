
#include "c_convert.h"

namespace panache
{

std::shared_ptr<panache::BasisSet> BasisSetFromArrays(std::shared_ptr<panache::Molecule> molecule,
        int ncenters,
        int * nshellspercenter,
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
                                                  i, cen, 0, panache::ShellInfo::PrimitiveType::Unnormalized));
            counter++;
        }

        shellmap.push_back(shellvec);

    }

    return std::shared_ptr<panache::BasisSet>(new panache::BasisSet(molecule, shellmap));

}

std::shared_ptr<panache::Molecule> MoleculeFromArrays(int ncenters, C_AtomCenter * atoms)
{
    std::shared_ptr<panache::Molecule> molecule(new panache::Molecule());

    for(auto i = 0; i < ncenters; i++)
    {
        molecule->add_atom(
            atoms[i].Z,
            atoms[i].center[0],
            atoms[i].center[1],
            atoms[i].center[2],
            atoms[i].symbol,
            atoms[i].mass,
            0);
    }

    molecule->reset_point_group("c1");
    molecule->set_full_point_group();
    return molecule;
}


} // close namespace panache

