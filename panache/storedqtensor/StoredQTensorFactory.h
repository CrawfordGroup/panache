/*! \file
 * \brief Create a three-index tensor storage object (header)
 * \ingroup storedqgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_STOREDQTENSORFACTORY_H
#define PANACHE_STOREDQTENSORFACTORY_H

#include <memory>
#include <string>

namespace panache
{

class StoredQTensor;

/*!
 * \brief Create and initialize a StoredQTensor object
 * \ingroup storedqgroup
 *
 *
 * \param [in] naux Number of auxiliary functions
 * \param [in] ndim1 Length along dimension 1
 * \param [in] ndim2 Length along dimension 2
 * \param [in] storeflags How to store (disk, memory, packed, etc)
 * \param [in] name Name of the tensor
 * \param [in] directory Directory for disk storage (if needed)
 * \return Pointer to a an object derived from StoredQTensor
 */
UniqueStoredQTensor StoredQTensorFactory(int naux, int ndim1, int ndim2, int storeflags,
        const std::string & name, const std::string & directory);


/*!
 * \brief Create a StoredQTensor object
 * \ingroup storedqgroup
 *
 *  This does not initialize the object since it isn't passed the sizes.
 *
 * \param [in] storeflags How to store (disk, memory, packed, etc)
 * \param [in] directory Directory for disk storage (if needed)
 * \return Pointer to a an object derived from StoredQTensor
 */
UniqueStoredQTensor StoredQTensorFactory(int storeflags, const std::string & directory);

} // close namespace panache

#endif
