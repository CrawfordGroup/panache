/*! \file
 *  \brief Conversions from C-style arrays to C++ objects (header)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_C_CONVERT_H
#define PANACHE_C_CONVERT_H

#include "panache/BasisSet.h"
#include "panache/Molecule.h"
#include "panache/c_interface.h"

namespace panache
{

    /*!
     * \brief Converts basis set information in C-style arrays to a BasisSet object
     *
     * \note Basis set coefficients should not be normalized
     *
     * \param [in] molecule Molecule object that this basis set is for
     * \param [in] ncenters Number of basis function centers
     * \param [in] nshellspercenter  Number of shells on each center.
     *                               Expected to be of length \p ncenters.
     * \param [in] shells  Information about each shell.
     *                     Length should be the sum of \p nshellspercenter.
     *
     * \return A BasisSet object with the given basis information.
     *
     */ 
    SharedBasisSet BasisSetFromArrays(SharedMolecule molecule,
            int_t ncenters,
            int_t * nshellspercenter,
            struct C_ShellInfo * shells);



    /*!
     * \brief Converts molecule information in C-style arrays to a Molecule object
     *
     * \param [in] ncenters Number of centers in this molecule
     * \param [in] atoms Information for each center. Expected to be of length \p ncenters. 
     *
     * \return A new molecule object. 
     *
     */ 
    SharedMolecule MoleculeFromArrays(int_t ncenters, C_AtomCenter * atoms);

} // close namespace panache

#endif
