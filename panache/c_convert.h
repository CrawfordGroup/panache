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
     * \param [in] molecule Molecule object that this basis set is for
     * \param [in] ncenters Number of basis function centers
     * \param [in] nshellspercenter  Number of shells on each center.
     *                               Expected to be of length \p ncenters.
     * \param [in] shells  Information about each shell.
     *                     Length should be the sum of \p nshellspercenter.
     * \param [in] normalized True if the basis is already normalized.
     *
     * \return A BasisSet object with the given basis information.
     *
     */ 
    std::shared_ptr<panache::BasisSet> BasisSetFromArrays(std::shared_ptr<panache::Molecule> molecule,
            int_t ncenters,
            int_t * nshellspercenter,
            struct C_ShellInfo * shells, bool normalized);



    /*!
     * \brief Converts molecule information in C-style arrays to a Molecule object
     *
     * \param [in] ncenters Number of centers in this molecule
     * \param [in] atoms Information for each center. Expected to be of length \p ncenters. 
     *
     * \return A new molecule object. 
     *
     */ 
    std::shared_ptr<panache::Molecule> MoleculeFromArrays(int_t ncenters, C_AtomCenter * atoms);

} // close namespace panache

#endif
