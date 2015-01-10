/*! \file
 * \brief Three-index tensor storage in memory (header)
 * \ingroup storedqgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_MEMORYQTENSOR_H
#define PANACHE_MEMORYQTENSOR_H

#include "panache/storedqtensor/LocalQTensor.h"

namespace panache
{

class DiskQTensor;

/*!
 *  \brief Class for storing a 3-index tensor in memory
 *  \ingroup storedqgroup
 */
class MemoryQTensor : public LocalQTensor
{
public:
    /*
     * \brief Construct with some basic information
     *
     * \param [in] storeflags How the tensor should be stored (packed, etc)
     * \param [in] name Some descriptive name
     */
    MemoryQTensor(int storeflags, const std::string & name);

    /*!
     * \brief Reads a tensor from disk completely into memory
     */
    MemoryQTensor(DiskQTensor * diskqt);  

protected:
    virtual void Write_(double * data, int nij, int ijstart);
    virtual void WriteByQ_(double * data, int nij, int ijstart);
    virtual void Read_(double * data, int nij, int ijstart);
    virtual void ReadByQ_(double * data, int nq, int qstart);
    virtual void Init_(void);

private:
    std::unique_ptr<double[]> data_;

};

} // close namespace panache

#endif

