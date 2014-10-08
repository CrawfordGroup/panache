/*! \file
 * \brief Three-index tensor storage on disk (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_DISKQTENSOR_H
#define PANACHE_DISKQTENSOR_H

#include <fstream>
#include <memory>

#include "panache/storedqtensor/LocalQTensor.h"

namespace panache
{

/*!
 *  \brief Class for storing a 3-index tensor on disk
 */
class DiskQTensor : public LocalQTensor
{
public:
    DiskQTensor();

protected:
    virtual void Reset_(void);
    virtual void Write_(double * data, int nij, int ijstart);
    virtual void WriteByQ_(double * data, int nij, int ijstart);
    virtual void Read_(double * data, int nij, int ijstart);
    virtual void ReadByQ_(double * data, int nq, int qstart);
    virtual void Clear_(void);
    virtual void Init_(void);

private:
    std::unique_ptr<std::fstream> file_;

    void OpenFile_(void);
    void CloseFile_(void);
};

} // close namespace panache

#endif


