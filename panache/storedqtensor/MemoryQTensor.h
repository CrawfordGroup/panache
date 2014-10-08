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

/*!
 *  \brief Class for storing a 3-index tensor in memory
 *  \ingroup storedqgroup
 */
class MemoryQTensor : public LocalQTensor
{
public:
    MemoryQTensor();

protected:
    virtual void Reset_(void);
    virtual void Write_(double * data, int nij, int ijstart);
    virtual void WriteByQ_(double * data, int nij, int ijstart);
    virtual void Read_(double * data, int nij, int ijstart);
    virtual void ReadByQ_(double * data, int nq, int qstart);
    virtual void Clear_(void);
    virtual void Init_(void);

private:
    std::unique_ptr<double[]> data_;

};

} // close namespace panache

#endif

