/*! \file
 * \brief Create a three-index tensor storage object (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/Flags.h"
#include "panache/Exception.h"

// All the different StoredQTensor types
#include "panache/tensorbackend/MemoryQTensor.h"
#include "panache/tensorbackend/DiskQTensor.h"

#ifdef PANACHE_CYCLOPS
#include "panache/tensorbackend/CyclopsQTensor.h"
#endif

namespace panache
{

UniqueStoredQTensor 
StoredQTensorFactory(int storeflags)
{
    if(storeflags & QSTORAGE_ONDISK)
        return UniqueStoredQTensor(new DiskQTensor());

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

    auto ptr = StoredQTensorFactory(storeflags);

    // convert name to filename if on disk
    if(storeflags & QSTORAGE_ONDISK)
    {
        std::string filename(directory);
        filename.append("/");
        filename.append(name);
        ptr->Init(naux, ndim1, ndim2, storeflags, filename);
    }
    else
        ptr->Init(naux, ndim1, ndim2, storeflags, name);

    return ptr;
}

} // close namespace panache
