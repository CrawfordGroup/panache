#ifndef PANACHE_C_CONVERT_H
#define PANACHE_C_CONVERT_H

#include "BasisSet.h"
#include "Molecule.h"
#include "c_interface.h"

namespace panache {

std::shared_ptr<panache::BasisSet> BasisSetFromArrays(std::shared_ptr<panache::Molecule> molecule,
        int ncenters,
        int * nshellspercenter,
        struct C_ShellInfo * shells);

std::shared_ptr<panache::Molecule> MoleculeFromArrays(int ncenters, C_AtomCenter * atoms);


} // close namespace panache

#endif
