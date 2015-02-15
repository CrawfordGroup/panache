/*! \file
 * \brief Three-index tensor storage on disk (header)
 * \ingroup storedqgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_DISKQTENSOR_H
#define PANACHE_DISKQTENSOR_H

#include <fstream>

#include "panache/storedqtensor/LocalQTensor.h"

namespace panache
{

class MemoryQTensor;

/*!
 *  \brief Class for storing a 3-index tensor on disk
 *  \ingroup storedqgroup
 */
class DiskQTensor : public LocalQTensor
{
public:
    /*
     * \brief Construct with some basic information
     *
     * \param [in] storeflags How the tensor should be stored (packed, etc)
     * \param [in] name Some descriptive name
     * \param [in] directory Directory where to store the files
     */
    DiskQTensor(int storeflags, const std::string & name, const std::string & directory);

    DiskQTensor(MemoryQTensor * memqt);


    virtual ~DiskQTensor();

protected:
    virtual void Write_(double * data, int nij, int ijstart);
    virtual void WriteByQ_(double * data, int nq, int qstart, bool ijpacked);
    virtual void Read_(double * data, int nij, int ijstart);
    virtual void ReadByQ_(double * data, int nq, int qstart);
    virtual void Init_(void);

private:
    std::unique_ptr<std::fstream> file_;

    int f_naux_; //!< naux on the dim file
    int f_ndim1_; //!< ndim1 on the dim file
    int f_ndim2_; //!< ndim2 on the dim file
    int f_ndim12_; //!< ndim12 on the dim file
    int f_packed_; //!< ispacked on the dim file
    int f_byq_; //!< byq on the dim file

    void ReadDimFile_(void);
    void WriteDimFile_(void);
    void OpenForReadWrite_(void);
    bool OpenForRead_(bool required);
};

} // close namespace panache

#endif


