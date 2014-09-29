/*! \file
 * \brief Density fitting tensor generation and manipulation (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_DFTENSOR_H
#define PANACHE_DFTENSOR_H

#include <vector>
#include <memory>

#include "panache/Timing.h"
#include "panache/Flags.h"
#include "panache/Iterator.h"

#ifdef PANACHE_CYCLOPS
#include <ctf.hpp>
#endif

namespace panache
{

class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;

class FittingMetric;
class TwoBodyAOInt;

namespace reorder
{
class Orderings;
class CNorm;
}

class ThreeIndexTensor
{
public:
    ThreeIndexTensor(const ThreeIndexTensor & other) = delete;
    ThreeIndexTensor(const ThreeIndexTensor && other) = delete;
    ThreeIndexTensor & operator=(const ThreeIndexTensor & other) = delete;



    /*!
     * \brief Constructor
     *
     * Initializes the basis set members, and constructs the creates the
     * FittingMetric
     *
     * \param [in] primary The primary basis set
     * \param [in] auxiliary The auxiliary (DF) basis set
     * \param [in] directory Full path to a directory to put scratch files
     * \param [in] nthreads Max number of threads to use
     */ 
    ThreeIndexTensor(SharedBasisSet primary,
             SharedBasisSet auxiliary,
             const std::string & directory,
             int nthreads);


    /*!
     * \brief Destructor
     *
     * Frees memory and stuff
     */
    ~ThreeIndexTensor();



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
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \param [out] naux Number of auxiliary basis functions
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
     * | Qso    | nso*(nso+1)/2    |               |
     * | Qmo    | nmo*(nmo+1)/2    |               |
     * | Qoo    | nocc*(nocc+1)/2  |               |
     * | Qov    |                  | nocc*nvir     |
     * | Qvv    | nvir*(nvir+1)/2  |               |
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
     * \param [in] qstart Iterator representing where to start
     * \return The number of batches actually stored in the buffer.
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
     * \param [in] qstart Iterator representing where to start
     * \return The number of batches actually stored in the buffer.
     */
    int GetBatch(int tensorflag, double * outbuf, int bufsize, IJIterator ijstart);




    template<typename ITTYPE>
    class IteratedQTensor
    {
        public:
            typedef std::function<int(int)> GetBatchFunc;

            IteratedQTensor(int tensorflag, double * buf, int bufsize,
                            int batchsize, const ITTYPE & it, GetBatchFunc gbf)
                : gbf_(gbf),it_(it),buf_(buf),bufsize_(bufsize),batchsize_(batchsize)
            {
                GetBatch();
            }

            // note that we don't own the pointer
            // so no need to deep copy
            IteratedQTensor(const IteratedQTensor & rhs) = default;

            operator bool() const { return nbatch_ > 0 && it_; }
            bool Valid(void) const { return nbatch_ > 0 && it_; }

            // Prefix increments only
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

            double * Get(void) { return curptr_; }

            double operator[](int i) { return curptr_[i]; }

            int Batchsize(void) const { return batchsize_; }

            int Index(void) const { return it_.Index(); }

        protected:
            const ITTYPE & Iterator(void) const { return it_; }

        private:
            GetBatchFunc gbf_;

            ITTYPE it_;
            double * curptr_;
            int nbatch_;   //!< Number of batches currently in the buffer
            int curbatch_; //!< The one I'm currently on

            double * buf_;
            int bufsize_;
            int batchsize_;

            void GetBatch(void)
            {
                nbatch_ = gbf_(it_.Index());
                curptr_ = buf_;
                curbatch_ = 0;
            }

    };

    class IteratedQTensorByQ : public IteratedQTensor<QIterator>
    {
        public:
            IteratedQTensorByQ(int tensorflag, double * buf, int bufsize,
                               int batchsize, const QIterator & it, GetBatchFunc gbf)
                    : IteratedQTensor(tensorflag, buf, bufsize, batchsize, it, gbf) { }


            int q(void) const { return Iterator().q(); } 
    };

    class IteratedQTensorByIJ : public IteratedQTensor<IJIterator>
    {
        public:
            IteratedQTensorByIJ(int tensorflag, double * buf, int bufsize,
                                int batchsize, const IJIterator & it, GetBatchFunc gbf)
                    : IteratedQTensor(tensorflag, buf, bufsize, batchsize, it, gbf) { }


            int i(void) const { return Iterator().i(); }
            int j(void) const { return Iterator().j(); }
            int ij(void) const { return Iterator().ij(); }
    };



    IteratedQTensorByQ IterateByQ(int tensorflag, double * buf, int bufsize);

    IteratedQTensorByIJ IterateByIJ(int tensorflag, double * buf, int bufsize);


private:

    /*! \name Basis set and matrix dimensions */
    ///@{  

    SharedBasisSet primary_;   //!< Primary basis set
    SharedBasisSet auxiliary_; //!< Auxiliary (DF) basis set

    int nso_;    //!< Number of SO (rows of Cmo_, number of primary basis functions)
    int naux_;   //!< Number of auxiliary basis functions
    int nso2_;   //!< Number of SO squared

    ///@}


    /*! \name Metric generation */
    ///@{

    std::shared_ptr<FittingMetric> fittingmetric_; //!< Fitting metric J

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
    ///@}



    /*! \brief Storage of the 3-index tensors */
    ///@{ Storage of the 3-index tensors

    class StoredQTensor
    {

    public:
        StoredQTensor(int naux, int ndim1, int ndim2, int storeflags);
        virtual ~StoredQTensor();
        int StoreFlags(void) const;
        int Read(double * data, int nij, int ijstart);
        int ReadByQ(double * data, int nq, int qstart);
        void Reset(void);
        void Clear(void);
        void Init(void);

        void GenDFQso(const std::shared_ptr<FittingMetric> & fit,
                      const SharedBasisSet primary,
                      const SharedBasisSet auxiliary,
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
        CumulativeTime gen_timer_; //!< Timer for the generation of this tensor
        CumulativeTime getijbatch_timer_; //!< Timer for getting batch by orbital index
        CumulativeTime getqbatch_timer_; //!< Timer for getting batch by q index

        // disable copying, etc
        StoredQTensor & operator=(const StoredQTensor & rhs) = delete;
        StoredQTensor(const StoredQTensor & rhs) = delete;
        StoredQTensor(const StoredQTensor && rhs) = delete;
    };


    class LocalQTensor : public StoredQTensor
    {
    public:
        LocalQTensor(int naux, int ndim1, int ndim2, int storeflags);

    protected:
        virtual void Write_(double * data, int nij, int ijstart) = 0;
        virtual void WriteByQ_(double * data, int nij, int ijstart) = 0;

        virtual void GenDFQso_(const std::shared_ptr<FittingMetric> & fit,
                             const SharedBasisSet primary,
                             const SharedBasisSet auxiliary,
                             int nthreads);

        virtual void Transform_(const std::vector<TransformMat> & left,
                                const std::vector<TransformMat> & right,
                                std::vector<StoredQTensor *> results,
                                int nthreads);
    };

    class DiskQTensor : public LocalQTensor
    {
    public:
        DiskQTensor(int naux, int ndim1, int ndim2, int storeflags, const std::string & filename);

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
        std::string filename_;

        void OpenFile_(void);
        void CloseFile_(void);
    };

    
    class MemoryQTensor : public LocalQTensor
    {
    public:
        MemoryQTensor(int naux, int ndim1, int ndim2, int storeflags);

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


    #ifdef PANACHE_CYCLOPS
    class CyclopsQTensor : public StoredQTensor
    {
    private:
        std::string name_;
        std::unique_ptr<CTF_Tensor> tensor_;

        void DecomposeIndex_(int64_t index, int & i, int & j, int & q);
        std::unique_ptr<CTF_Matrix> FillWithMatrix_(double * mat, int nrow, int ncol, int sym, const char * name);
 

    protected:
        virtual void Reset_(void);
        virtual void Read_(double * data, int nij, int ijstart);
        virtual void ReadByQ_(double * data, int nq, int qstart);
        virtual void Clear_(void);
        virtual void Init_(void);

        virtual void GenDFQso_(const std::shared_ptr<FittingMetric> & fit,
                               const SharedBasisSet primary,
                               const SharedBasisSet auxiliary,
                               int nthreads);

        virtual void Transform_(const std::vector<TransformMat> & left,
                                const std::vector<TransformMat> & right,
                                std::vector<StoredQTensor *> results,
                                int nthreads);

    public:
        CyclopsQTensor(int naux, int ndim1, int ndim2, int storeflags, const std::string & name);
    };
    #endif

    /*!
     * \brief Create a StoredQTensor object
     *
     * \param [in] naux Number of auxiliary functions
     * \param [in] ndim1 Length along dimension 1
     * \param [in] ndim2 Length along dimension 2
     * \param [in] storeflags How to store (disk, memory, packed, etc)
     * \param [in] name Name of the tensor
     * \return Pointer to a an object derived from StoredQTensor
     */
    std::unique_ptr<StoredQTensor> StoredQTensorFactory(int naux, int ndim1, int ndim2,
                                                        int storeflags,
                                                        const std::string & name);

    std::unique_ptr<StoredQTensor> qso_;  //!< Qso matrix
    std::unique_ptr<StoredQTensor> qmo_;  //!< Qmo matrix
    std::unique_ptr<StoredQTensor> qoo_;  //!< Qoo matrix
    std::unique_ptr<StoredQTensor> qov_;  //!< Qov matrix
    std::unique_ptr<StoredQTensor> qvv_;  //!< Qvv matrix

    ///@}


    CumulativeTime timer_genqtensors_;   //!< Total time spent in GenQTensors()

    int nthreads_;  //!< Number of threads to use

    std::string directory_;  //!< Directory to use to store matrices on disk (if requested)



    /*!
     * \brief Generates the basic Qso matrix
     *
     * See \ref theory_page for what Qso actually is, and memory_sec for more information
     * about memory.
     *
     * \param [in] qflags How to generate QSO (DF, Cholesky) (see Flags.h)
     * \param [in] storeflags How to store the matrix (see Flags.h)
     */
    void GenQso(int qflags, int storeflags);


    /*!
     * \brief Obtains a stored tensor object corresponding to a tensor flag
     *
     * \param [in] tensorflag Flag for the tensor
     * \return Stored object for the tensor
     */
    std::unique_ptr<StoredQTensor> & ResolveTensorFlag(int tensorflag);


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
    void PrintTimer(const char * name, const std::unique_ptr<StoredQTensor> & q) const;

};

}

#endif //PANACHE_DFTENSOR_H

