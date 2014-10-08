/*! \file
 * \brief Generic three-index tensor storage (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_STOREDQTENSOR_H
#define PANACHE_STOREDQTENSOR_H

#include <memory>
#include <string>
#include <vector>

#include "panache/Timing.h"

namespace panache
{

class FittingMetric;
class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;


/*!
 *  \brief Generic/Abstract interface for storing a 3-index tensor
 */
class StoredQTensor
{
public:
    StoredQTensor();
    virtual ~StoredQTensor();
    int StoreFlags(void) const;
    int Read(double * data, int nij, int ijstart);
    int ReadByQ(double * data, int nq, int qstart);
    void Reset(void);
    void Clear(void);
    void Init(int naux, int ndim1, int ndim2, int storeflags, const std::string & name);

    void GenDFQso(const std::shared_ptr<FittingMetric> & fit,
                  const SharedBasisSet primary,
                  const SharedBasisSet auxiliary,
                  int nthreads);

    void GenCHQso(const SharedBasisSet primary,
                  double delta,
                  int storeflags,
                  int nthreads);

    typedef std::pair<double *, int> TransformMat;

    void Transform(const std::vector<TransformMat> & left,
                   const std::vector<TransformMat> & right,
                   std::vector<StoredQTensor *> results,
                   int nthreads);

    CumulativeTime & GenTimer(void);
    CumulativeTime & GetBatchTimer(void);
    CumulativeTime & GetQBatchTimer(void);

    int naux(void) const;
    int ndim1(void) const;
    int ndim2(void) const;
    int ndim12(void) const;
    int packed(void) const;
    int byq(void) const;
    int calcindex(int i, int j) const;
    const std::string & name(void) const;

protected:
    virtual void Init_(void) = 0;
    virtual void Reset_(void) = 0;
    virtual void Clear_() = 0;
    virtual void Read_(double * data, int nij, int ijstart) = 0;
    virtual void ReadByQ_(double * data, int nq, int qstart) = 0;

    virtual void GenDFQso_(const std::shared_ptr<FittingMetric> & fit,
                           const SharedBasisSet primary,
                           const SharedBasisSet auxiliary,
                           int nthreads) = 0;

    virtual void GenCHQso_(const SharedBasisSet primary,
                           double delta,
                           int storeflags,
                           int nthreads) = 0;


    virtual void Transform_(const std::vector<TransformMat> & left,
                            const std::vector<TransformMat> & right,
                            std::vector<StoredQTensor *> results,
                            int nthreads) = 0;

    int storesize(void) const;

private:
    int naux_;   //!< Number of auxiliary functions
    int ndim1_;  //!< Length of index 1
    int ndim2_;  //!< Length of index 2
    int ndim12_; //!< Combined size of index 1 and 2 (depends on packing)
    int storeflags_; //!< How is this tensor stored
    std::string name_;
    CumulativeTime gen_timer_; //!< Timer for the generation of this tensor
    CumulativeTime getijbatch_timer_; //!< Timer for getting batch by orbital index
    CumulativeTime getqbatch_timer_; //!< Timer for getting batch by q index

    // disable copying, etc
    StoredQTensor & operator=(const StoredQTensor & rhs) = delete;
    StoredQTensor(const StoredQTensor & rhs) = delete;
    StoredQTensor(const StoredQTensor && rhs) = delete;
};

} // close namespace panache

#endif

