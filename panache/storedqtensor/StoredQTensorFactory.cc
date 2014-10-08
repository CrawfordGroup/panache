#include <memory>
#include <string>

#include "panache/Flags.h"
#include "panache/Exception.h"
#include "panache/storedqtensor/MemoryQTensor.h"
#include "panache/storedqtensor/DiskQTensor.h"

#ifdef PANACHE_CYCLOPS
#include "panache/storedqtensor/CyclopsQTensor.h"
#endif

namespace panache
{

std::unique_ptr<StoredQTensor> 
StoredQTensorFactory(int storeflags)
{
    if(storeflags & QSTORAGE_ONDISK)
        return std::unique_ptr<StoredQTensor>(new DiskQTensor());

    #ifdef PANACHE_CYCLOPS
    else if(storeflags & QSTORAGE_CYCLOPS)
        return std::unique_ptr<StoredQTensor>(new CyclopsQTensor());
    #endif

    else
        return std::unique_ptr<StoredQTensor>(new MemoryQTensor());
}

std::unique_ptr<StoredQTensor> 
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
