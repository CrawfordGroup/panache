/*! \file
 * \brief Create a three-index tensor storage object (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/Flags.h"
#include "panache/Exception.h"

// All the different StoredQTensor types
#include "panache/storedqtensor/MemoryQTensor.h"
#include "panache/storedqtensor/DiskQTensor.h"

#ifdef PANACHE_CYCLOPS
#include "panache/storedqtensor/CyclopsQTensor.h"
#endif

namespace panache
{

UniqueStoredQTensor 
StoredQTensorFactory(int storeflags, const std::string & directory)
{
    if(storeflags & QSTORAGE_ONDISK)
        return UniqueStoredQTensor(new DiskQTensor(directory));

    #ifdef PANACHE_CYCLOPS
    else if(storeflags & QSTORAGE_CYCLOPS)
        return UniqueStoredQTensor(new CyclopsQTensor());
    #endif

    else
        return UniqueStoredQTensor(new MemoryQTensor());
}

UniqueStoredQTensor 
StoredQTensorFactory(int naux, int ndim1, int ndim2,
                     int storeflags, const std::string & name,
                     const std::string & directory)
{
    if(name == "")
        throw RuntimeError("NO NAME SPECIFIED");

    auto ptr = StoredQTensorFactory(storeflags, directory);

    // convert name to filename if on disk
    ptr->Init(naux, ndim1, ndim2, storeflags, name);

    return ptr;
}

} // close namespace panache
