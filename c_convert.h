#ifndef PANACHE_C_CONVERT_H
#define PANACHE_C_CONVERT_H

#include "BasisSet.h"
#include "Molecule.h"
#include "c_interface.h"

namespace panache {

std::shared_ptr<panache::BasisSet> BasisSetFromArrays(std::shared_ptr<panache::Molecule> molecule,
        INTTYPE ncenters,
        INTTYPE * nshellspercenter,
        struct C_ShellInfo * shells);

std::shared_ptr<panache::Molecule> MoleculeFromArrays(INTTYPE ncenters, C_AtomCenter * atoms);


} // close namespace panache

#endif
