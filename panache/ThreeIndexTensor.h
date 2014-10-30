/*! \file
 * \brief Generic three-index tensor generation and manipulation (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_THREEINDEXTENSOR_H
#define PANACHE_THREEINDEXTENSOR_H

#include <vector>
#include <memory>
#include <string>

#include "panache/Timing.h"
#include "panache/Flags.h"
#include "panache/Iterator.h"

namespace panache
{

class FittingMetric;
class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;
class TwoBodyAOInt;
class StoredQTensor;
typedef std::unique_ptr<StoredQTensor> UniqueStoredQTensor;

namespace reorder
{
class Orderings;
class CNorm;
}


/*!
 *   \brief A generic three-index tensor
 *
 *   Specific 3-index tensors (Cholesky, Density-fitted) are
 *   generated through derived classes. The members of this class are
 *   used to store and transorm those tensors.
 *
 *   This class also stores the MO coefficient matrix and threading information.
 */
class ThreeIndexTensor
{
public:
    /*!
     * \brief Constructor
     *
     * Initializes the basis set shared pointer and other information
     *
     * \param [in] primary The primary basis set
     * \param [in] directory Full path to a directory to put scratch files
     * \param [in] qtype The type of tensor stored (Cholesky, DF)
     * \param [in] nthreads Max number of threads to use
     */ 
    ThreeIndexTensor(SharedBasisSet primary,
             const std::string & directory,
             int qtype,
             int nthreads);



    /*!
     * \brief Destructor
     *
     * Frees memory and stuff
     */
    virtual ~ThreeIndexTensor();


    // Don't need these
    ThreeIndexTensor(const ThreeIndexTensor & other) = delete;
    ThreeIndexTensor(const ThreeIndexTensor && other) = delete;
    ThreeIndexTensor & operator=(const ThreeIndexTensor & other) = delete;


    /*!
     * \brief Sets the C matrix (so-ao matrix) for use in generating Qmo, Qov, Qoo, and Qvv
     *
     * The matrix is expected be nso x nmo (MOs in the columns) in row-major order.
     * If it is nmo x nso, or the matrix is in column major order, set \p cmo_is_trans
     * to true.
     *
     * The matrix is copied by the PANACHE code, so it can be safely deleted or otherwise
     * changed after calling this function.
     *
     * \param [in] cmo Pointer to a nso x nmo matrix representing the MO coefficients
     * \param [in] nmo Number of MOs in this C matrix
     * \param [in] cmo_is_trans Set to non-zero if the matrix is the transpose (nmo x nso) or
     *                          is in column-major order.
     * \param [in] order Ordering of the basis functions
     */
    void SetCMatrix(double * cmo, int nmo, bool cmo_is_trans, int order = BSORDER_PSI4);



    /*!
     * \brief Sets the number of occupied and virtual orbitals.
     *
     * Number of virtual orbitals is taken to be the remainder after the occupied.
     * Used by Qov, etc.
     *
     * Frozen orbitals should not be counted in \p nocc.
     *
     * \note You must set the C Matrix first before calling (see SetCMatrix())
     *
     * \param [in] nocc Number of (non-frozen) occupied orbitals
     * \param [in] nfroz Number of frozen orbitals
     */
    void SetNOcc(int nocc, int nfroz = 0);




    /*!
     * \brief Sets the maximum number of (OpenMP) threads used
     *
     * Set to zero to use the value of the environment variable OMP_NUM_THREAD (or
     * set by omp_num_threads, or the default for this machine).
     *
     * \param [in] nthread Max number of threads to use
     * \return The max number of threads that will actually be used (ie if \p nthread is zero).
     */
    int SetNThread(int nthread);



    /*!
     * \brief Prints out timing information collected so far
     *
     * All times are cumulative for all operations. The output must be set
     *  first (See Output.h)
     */
    void PrintTimings(void) const;



    /*!
     * \brief Generates various 3-index tensors
     *
     * For the \p qflags and \p storeflags parameters, see Flags.h. For example, to calculate the
     * Qmo and Qov tensors on disk,
     *
     * \code{.cpp}
     * dft.GenQTensors(QGEN_QMO | QGEN_QOV, QSTORAGE_ONDISK);
     * \endcode
     *
     * Default is QSTORAGE_INMEM and not to store with QSTORAGE_BYQ
     *
     * To calculate just Qso, set \p qflags = QGEN_QSO
     *
     * \note The Qso matrix is always stored with QSTORAGE_BYQ
     * \warning Be sure to set the C-Matrix first and number of occupied orbitals first
     *          if qflags contains more than QGEN_QSO
     *
     * \param [in] qflags A combination of flags specifying which tensors to generate
     * \param [in] storeflags How to store the matrix
     */
    void GenQTensors(int qflags, int storeflags);



    /*!
     * \brief Obtain the dimensions of a tensor
     *
     * The 3-index tensor is accessed in row-major order
     * and has dimensions of (naux, ndim1, ndim2) if storage is
     * QSTORAGE_BYQ or (ndim1, ndim2, naux) otherwise.
     *
     *
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \param [out] naux Number of auxiliary index elements
     * \param [out] ndim1 First dimension
     * \param [out] ndim2 Second dimension
     * \return Total tensor size (depends on packing)
     */
    int TensorDimensions(int tensorflag, int & naux, int & ndim1, int & ndim2);



    /*!
     * \brief Obtain the batch size for GetQBatch()
     *
     * The size will be as follows
     *
     * | Tensor | Packed Size      | Unpacked Size |
     * |--------|------------------|---------------|
     * | Qso    | nso*(nso+1)/2    | nso*nso       |
     * | Qmo    | nmo*(nmo+1)/2    | nmo*nmo       |
     * | Qoo    | nocc*(nocc+1)/2  | nocc*nocc     |
     * | Qov    |                  | nocc*nvir     |
     * | Qvv    | nvir*(nvir+1)/2  | nvir*nvir     |
     *
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \return Size of batches returned by GetQBatch
     */
    int QBatchSize(int tensorflag);



    /*!
     * \brief Obtain the batch size for GetBatch()
     *
     * The size will always be naux (number of auxiliary basis functions)
     *
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \return Size of batches returned by GetBatch
     */
    int BatchSize(int tensorflag);



    /*!
     * \brief See if a particular tensor is stored packed
     *
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \return True if the tensor is in packed storage
     */
    bool IsPacked(int tensorflag);



    /*!
     * \brief Calculate a combined orbital index
     *
     * Depends on packing
     * 
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \param [in] i First orbital index
     * \param [in] j Second orbital index
     * \return ij, depending on packing
     */
    int CalcIndex(int tensorflag, int i, int j);


    /*!
     * \brief Retrieves a batch of a 3-index tensor 
     *
     * See \ref theory_page for what these tensors actually are, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*batchsize elements with
     * the index of the auxiliary basis function as the slowest index.
     *
     * The batchsize can be obtained using QBatchSize()
     *
     * Call this and process the batches, incrementing qstart by the return value,
     * until this function returns zero.
     *
     * \param [in] tensorflag Which tensor to get (see Flags.h)
     * \param [in] outbuf Memory location to store the batch of tensors
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \return The number of batches actually stored in the buffer.
     */
    int GetQBatch(int tensorflag, double * outbuf, int bufsize, int qstart);


    /*!
     * \copybrief   GetQBatch
     * \copydetails GetQBatch
     */
    int GetQBatch(int tensorflag, double * outbuf, int bufsize, QIterator qstart);


    /*!
     * \brief Retrieves a batch of a 3-index tensor 
     *
     * See \ref theory_page for what these tensors actually are, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*naux elements with the combined orbital index
     * as the slowest index.
     *
     * The batchsize can be obtained using BatchSize()
     *
     * Call this and process the batches, incrementing qstart by the return value,
     * until this function returns zero.
     *
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] tensorflag Which tensor to get (see Flags.h)
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] ijstart The starting value of the ij index
     * \return The number of batches actually stored in the buffer.
     */
    int GetBatch(int tensorflag, double * outbuf, int bufsize, int ijstart);


    /*!
     * \copybrief   GetBatch
     * \copydetails GetBatch
     */
    int GetBatch(int tensorflag, double * outbuf, int bufsize, IJIterator ijstart);



    /*!
     * \brief A class for iterating over a three-index tensor
     *
     * Iteration can happen "by Q" or "by IJ", depending on ITTYPE
     *
     * This class is not meant to be used directly, but through 
     * IteratedQTensorByQ and IteratedQTensorByIJ
     *
     * Incrementing can only be done with a prefix ++ (ie, ++it, not it++)
     *
     * Use Get() to get a pointer to the current data for the indices stored in the
     * iterator.
     **/
    template<typename ITTYPE>
    class IteratedQTensor
    {
        public:
            /// A function that takes an index, and returns the size of the batch retrieved 
            typedef std::function<int(int)> GetBatchFunc;

            /*!
             * \brief Constructor
             *
             * \param [in] tensorflag Which tensor to iterator over
             * \param [in] buf Where to place the retrieved batches
             * \param [in] bufsize The size of the buffer \p buf (in number of doubles)
             * \param [in] batchsize The size of a batch of the tensor
             * \param [in] it Base iterator
             * \param [in] gbf Function to obtain the batches
             */
            IteratedQTensor(int tensorflag, double * buf, int bufsize,
                            int batchsize, const ITTYPE & it, GetBatchFunc gbf)
                : gbf_(gbf),it_(it),buf_(buf),bufsize_(bufsize),batchsize_(batchsize)
            {
                GetBatch();
            }

            // note that we don't own the pointer
            // so no need to deep copy
            IteratedQTensor(const IteratedQTensor & rhs) = default;

            /*!
             *  \brief Is the iterator in a valid state or not
             */ 
            operator bool() const { return nbatch_ > 0 && it_; }

            /*!
             *  \brief Is the iterator in a valid state or not
             */
            bool Valid(void) const { return nbatch_ > 0 && it_; }


            // Prefix increments only
            /*!
             *  \brief Increment and obtain the next batch of the three-index tensor
             *
             *  If there is no more batches when this is called, the iterator will be
             *  placed in the 'invalid' state (and Valid() and bool() will return false).
             */
            IteratedQTensor & operator++()
            {
                ++it_;
                if(!it_)
                    return *this; // finish now

                curbatch_++;
                if(curbatch_ >= nbatch_)
                    GetBatch();
                else
                    curptr_ += batchsize_;

                return *this;
            }

            /*!
             * \brief Pointer to the data for the given indices
             */
            double * Get(void) { return curptr_; }


            /*!
             * \brief Obtain data stored in a batch
             */
            double operator[](int i) { return curptr_[i]; }

            /*!
             * \brief Obtain the size of a batch
             */
            int Batchsize(void) const { return batchsize_; }

            /*!
             * \brief Get the current (flattened) index for the current batch
             */ 
            int Index(void) const { return it_.Index(); }

        protected:

            /*!
             * \brief Get the stored iterator
             */
            const ITTYPE & Iterator(void) const { return it_; }

        private:
            GetBatchFunc gbf_;  //!< Function pointer to get batches

            ITTYPE it_;         //!< Stored iterator object
            double * curptr_;   //!< Pointer to data corresponding to the current iterator value
            int nbatch_;        //!< Number of batches currently in the buffer
            int curbatch_;      //!< The batch I'm currently on

            double * buf_;      //!< Where to store the data
            int bufsize_;       //!< Size of the buffer
            int batchsize_;     //!< Size of a single batch


            /*!
             * \brief Actually Get a batch (through gbf_)
             */
            void GetBatch(void)
            {
                nbatch_ = gbf_(it_.Index());
                curptr_ = buf_;
                curbatch_ = 0;
            }

    };


    /*!
     * \brief Iterate over a three-index tensor by Q
     */ 
    class IteratedQTensorByQ : public IteratedQTensor<QIterator>
    {
        public:
            /*!
             * \brief Constructor
             *
             * \param [in] tensorflag Which tensor to iterator over
             * \param [in] buf Where to place the retrieved batches
             * \param [in] bufsize The size of the buffer \p buf (in number of doubles)
             * \param [in] batchsize The size of a batch of the tensor
             * \param [in] it Base iterator
             * \param [in] gbf Function to obtain the batches
             */
            IteratedQTensorByQ(int tensorflag, double * buf, int bufsize,
                               int batchsize, const QIterator & it, GetBatchFunc gbf)
                    : IteratedQTensor(tensorflag, buf, bufsize, batchsize, it, gbf) { }


            /*!
             * \brief Get the current value of the auxiliary index
             */ 
            int q(void) const { return Iterator().q(); } 
    };


    /*!
     * \brief Iterate over a three-index tensor by orbital index
     */ 
    class IteratedQTensorByIJ : public IteratedQTensor<IJIterator>
    {
        public:
            /*!
             * \brief Constructor
             *
             * \param [in] tensorflag Which tensor to iterator over
             * \param [in] buf Where to place the retrieved batches
             * \param [in] bufsize The size of the buffer \p buf (in number of doubles)
             * \param [in] batchsize The size of a batch of the tensor
             * \param [in] it Base iterator
             * \param [in] gbf Function to obtain the batches
             */
            IteratedQTensorByIJ(int tensorflag, double * buf, int bufsize,
                                int batchsize, const IJIterator & it, GetBatchFunc gbf)
                    : IteratedQTensor(tensorflag, buf, bufsize, batchsize, it, gbf) { }


            /*!
             * \brief Get the current value of the first index
             */ 
            int i(void) const { return Iterator().i(); }

            /*!
             * \brief Get the current value of the second index
             */ 
            int j(void) const { return Iterator().j(); }

            /*!
             * \brief Get the current value of the combined orbital index
             */ 
            int ij(void) const { return Iterator().ij(); }
    };



    /*!
     * \brief Obtain an iterator to a specific stored three-index tensor
     *
     * Iteration will happen over the auxiliary index
     */
    IteratedQTensorByQ IterateByQ(int tensorflag, double * buf, int bufsize);


    /*!
     * \brief Obtain an iterator to a specific stored three-index tensor
     *
     * Iteration will happen over the combined orbital index
     */
    IteratedQTensorByIJ IterateByIJ(int tensorflag, double * buf, int bufsize);


protected:


    /*!
     * \brief Generate the base Qso tensor
     * \param [in] storeflags How to store (disk, memory, packed, etc)
     */
    virtual UniqueStoredQTensor GenQso(int storeflags) const = 0;

    ///@}


    /*! \name Basis set and matrix dimensions */
    ///@{  

    SharedBasisSet primary_;   //!< Primary basis set
    SharedBasisSet auxiliary_; //!< Auxiliary (DF) basis set

    int nso_;    //!< Number of SO (rows of Cmo_, number of primary basis functions)
    int nso2_;   //!< Number of SO squared

    ///@}


    /*! \name C matrix and dimensions */
    ///@{  
    std::unique_ptr<double[]> Cmo_;  //!< C matrix (nso x nmo)
    std::unique_ptr<double[]> Cmo_occ_;  //!< C matrix (occupied part, nso*nocc)
    std::unique_ptr<double[]> Cmo_vir_;  //!< C matrix (virtual part, nso*nvir)

    int nmo_;  //!< Number of MO (columns of Cmo_)
    int nmo2_; //!< Number of MO squared
    int nocc_; //!< Number of (non-frozen) occupied orbitals
    int nfroz_; //!< Number of frozen orbitals
    int nvir_; //!< Number of virtual orbitals
    int nsotri_; //!< Packed nso*nso symmetric matrix
    ///@}

    int nthreads_;  //!< Number of threads to use

    std::string directory_;  //!< Directory to use to store matrices on disk (if requested)


private:
    int qtype_;  //!< Type of tensor (DF, CH, etc. See Flags.h)


    /*!
     * \brief Splits the C matrix into occupied and virtual matrices
     *
     * Fills in Cmo_occ_ and Cmo_vir_
     *
     * nocc_, nvir_ nmo_, and nso_ must be set first!
     */
    void SplitCMat(void);


    /*!
     * \brief Reorders the whole C matrix into the specified ordering, including
     *        some renormalization if needed
     * 
     * \param [in] order The order to use
     * \param [in] cnorm Normalization factors to multiply by
     */
    void ReorderCMat(const reorder::Orderings & order, const reorder::CNorm & cnorm);




    ///@{ \name Q Tensor Storage
    UniqueStoredQTensor qso_;  //!< Qso matrix
    UniqueStoredQTensor qmo_;  //!< Qmo matrix
    UniqueStoredQTensor qoo_;  //!< Qoo matrix
    UniqueStoredQTensor qov_;  //!< Qov matrix
    UniqueStoredQTensor qvv_;  //!< Qvv matrix
    ///@}


    CumulativeTime timer_genqtensors_;   //!< Total time spent in GenQTensors()


    /*!
     * \brief Obtains a stored tensor object corresponding to a tensor flag
     *
     * \param [in] tensorflag Flag for the tensor
     * \return Stored object for the tensor
     */
    UniqueStoredQTensor & ResolveTensorFlag(int tensorflag);


    /*!
     * \brief Base function for getting batches by Q
     * 
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \param [in] qt Pointer to the StoredQTensor to retrieve
     */
    int GetQBatch_Base(double * outbuf, int bufsize, int qstart, StoredQTensor * qt);


    /*!
     * \brief Base function for getting batches by orbital index
     * 
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] ijstart Starting combined orbital index
     * \param [in] qt Pointer to the StoredQTensor to retrieve
     */
    int GetBatch_Base(double * outbuf, int bufsize, int ijstart, StoredQTensor * qt);

    

    /*!
     * \brief Prints the timer information for a given StoredQTensor
     *
     * \param [in] name Descriptive name of the tensor
     * \param [in] q Tensor whose timings to print
     * 
     */
    void PrintTimer(const char * name, const UniqueStoredQTensor & q) const;

};

} // close namespace panache

#endif //PANACHE_DFTENSOR_H

