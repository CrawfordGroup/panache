/*! \file
 * \brief Density fitting tensor generation and manipulation (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_DFTENSOR2_H
#define PANACHE_DFTENSOR2_H

#include <fstream>

#include "panache/Timing.h"

#include <future>


namespace panache
{

class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;

class FittingMetric;

namespace reorder
{
class Orderings;
}

class DFTensor2
{

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
    int nocc_; //!< Number of occupied orbitals
    int nvir_; //!< Number of virtual orbitals

    /*!
     * \brief Splits the C matrix into occupied and virtual matrices
     *
     * Fills in Cmo_occ_ and Cmo_vir_
     *
     * nocc_, nvir_ nmo_, and nso_ must be set first!
     */
    void SplitCMat(void);

    /*!
     * \brief Reorders the whole C matrix into the specified ordering
     * 
     * \param [in] order The order to use
     */
    void ReorderCMat(reorder::Orderings & order);
    ///@}



    class enum QStorage
    {
        INMEM,
        ONDISK,
        ONFLY
    };



    struct StoredQTensor
    {
        int ndim1_;
        int ndim2_;
        QStorage storetype_;

        // for ONDISK
        string filename_;
        std::unique_ptr<std::fstream> file_;
        int curij_;

        // for INMEM
        unique_ptr<double *> data_;

        void OpenFile(void)
        {
            if(file_is_open())
                return;
   

            if(filename.length() == 0)
                throw RuntimeError("Error - no file specified!");
 
            file_ = std::unique_ptr<std::fstream>(new std::fstream(filename_.c_str()),
                                                  std::fstream::in | std::fstream::out |
                                                  std::fstream::binary | std::fstream::trunc );
            if(!file_->is_open())
                throw RuntimeError(filename_);

            file_->exceptions(std::fstream::failbit | std::fstream::badbit | std::fstream::eofbit);
            curij_ = 0;
        }

        void CloseFile(void)
        {
            if(file_ && file_->is_open())
            {
                file_->close();
                matfile_.reset();
            }
        }


        void ResetFile(void)
        {
            if(file_)
            {
                file_->seekg(0);
                file_->seekp(0);
                curij_ = 0;
            }
        }

        void Reset(void)
        {
            curij_ = 0;
            if(storetype_ == QStorage::ONDISK)
                ResetFile();
        }
    };


    StoredQTensor qso_;
    StoredQTensor qmo_;
    StoredQTensor qoo_;
    StoredQTensor qsv_;
    StoredQTensor qvv_;

    ///@}








    /*! \name Memory, Buffers and placeholders for transformations */
    ///@{

    std::unique_ptr<double[]> qc_;  //!< Holds C(T) Q or QC intermediate
    std::unique_ptr<double[]> q_;   //!< Holds a batch of Q (packed storage)

    ///@}



    /*! \name Batch retreival and processing */
    ///@{

    int curq_;  //!< next Q to be (batch) read


    /*!
     * \brief Retrieves batches and stores then in DFTensor2::q_
     *
     * The number of Q returned may be less than \p ntoget (ie there are no more
     * Q left)
     *
     * \param [in] ntoget Number of Q in the batch (max to get)
     * \return Actual number of Q stored in DFTensor2::q_
     */
    int GetBatch_Base(int ntoget);

    /*!
     * \brief Obtains batches with a tranformation applied
     *
     * Transformation is the standard \f$ U_l^\dagger Q U_r \f$
     * where \f$ U_l \f$ and \f$ U_r \f$ are left and right transformation
     * matrices.
     *
     * Results are stored in outbuffer_
     *
     * \param [in] left   Left matrix
     * \param [in] lncols Number of columns in \p left
     * \param [in] right  Right matrix
     * \param [in] rncols Number of columns in \p right
     * \param [in] timer  A timer for this functionality
     * \param [in] timername  Some descriptive name for the timer
     * \param [in] nthreads Number of threads to use for the transformation
     * \return Number of Q stored in outbuffer_
     */
    int GetBatch_transform(double * left, int lncols, 
                           double * right, int rncols,
                           Timer & timer, const char * timername,
                           int nthreads);

    ///@}


    /*! \name Timing */
    ///@{
    Timer timer_genqso;        //!< Total time spent in GenQso()
    Timer timer_getbatch_qso;  //!< Total time spent in GetBatch_Qso()
    Timer timer_getbatch_qmo;  //!< Total time spent in GetBatch_Qmo()
    Timer timer_getbatch_qoo;  //!< Total time spent in GetBatch_Qoo()
    Timer timer_getbatch_qvv;  //!< Total time spent in GetBatch_Qvv()
    Timer timer_getbatch_qov;  //!< Total time spent in GetBatch_Qov()

    ///@}



    /*! \name Threading and Prefetching */
    ///@{

    int nthreads_;  //!< Number of threads to use

    std::unique_ptr<double[]> q2_;  //!< Buffer for prefetching from disk
    std::future<int> fill_future_;  //!< Asynchronous future object

    ///@}




public:
    enum class BSOrder
    {
        Psi4,
        GAMESS
    };


    DFTensor2 & operator=(const DFTensor2 & other) = delete;
    DFTensor2(const DFTensor2 & other) = delete;

    /*!
     * \brief Constructor
     *
     * \param [in] primary The primary basis set
     * \param [in] auxiliary The auxiliary (DF) basis set
     * \param [in] filename Full path to a file for disk storage (may not be used)
     * \param [in] nthreads Max number of threads to use
     */ 
    DFTensor2(SharedBasisSet primary,
             SharedBasisSet auxiliary,
             const std::string & filename,
             int nthreads);

    /*!
     * \brief Destructor
     *
     * Frees memory and stuff
     */
    ~DFTensor2();


    /*!
     * \brief Queries information about the expected matrix dimensions
     *
     * Useful for determining buffer sizes or determining if Qso should be placed in memory.
     * The size of Qso (unpacked) will be naux * nso2, and batches of Qso will be read in
     * multiples of nso2 (for GetBatch_Qso()).
     *
     * \param [out] naux Number of auxiliary basis functions
     * \param [out] nso2 Number of primary basis functions squared (nso*nso)
     * \return Total size of an unpacked Qso tensor (ie naux * nso2)
     */
    int QsoDimensions(int & naux, int & nso2);



    /*!
     * \brief Generates the basic Qso matrix
     *
     * See \ref theory_page for what Qso actually is, and memory_sec for more information
     * about memory.
     *
     * \param [in] inmem If nonzero, store the Qso matrix in memoryof auxiliary basis functions
     */
    void GenQso(bool inmem);



    /*!
     * \brief Sets the C matrix (so-ao matrix) for use in generating Qmo, Qov, Qoo, and Qvv
     *
     * The matrix is expected be nso x nmo (MOs in the columns) in row-major order.
     * If it is nmo x nso, or the matrix is in column major order, set \p cmo_is_trans.
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
    void SetCMatrix(double * cmo, int nmo, bool cmo_is_trans, BSOrder order = BSOrder::Psi4);



    /*!
     * \brief Sets the number of occupied and virtual orbitals.
     *
     * Number of virtual orbitals is taken to be the remainder after the occupied.
     * Used by Qov, etc.
     *
     * \note You must set the C Matrix first before calling (see SetCMatrix())
     *
     * \param [in] nocc Number of occupied orbitals
     */
    void SetNOcc(int nocc);




    /*!
     * \brief Sets the buffer used for storing batches of Qso, Qmo, Qov, Qoo, Qvv
     *
     * Batches are read in multiples of:
     *
     * - nso2 (GetBatch_Qso(), see QsoDimensions()),
     * - nmo*nmo (GetBatch_Qmo())
     * - nocc*nvir (GetBatch_Qov())
     * - nocc*nocc (GetBatch_Qoo())
     * - nvir*nvir (GetBatch_Qvv())
     *
     * How many can fit in the buffer is determined automatically
     * from the matsize parameter. Any 'left over' buffer space is not used.
     *
     * \param [in] buf A pointer to memory for a buffer (of \p bufsize size)
     * \param [in] size Number of elements in \p buffer (not number of bytes)
     */
    void SetOutputBuffer(double * buf, long int size);





    /*!
     * \brief Retrieves a batch of Qso
     *
     * The batches are stored in the matrix set by SetOutputBuffer().
     * See \ref theory_page for what Qso actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nso2 elements (see QsoDimensions()).
     *
     * Call this and process the batches until this function returns zero.
     *
     * \return The number of batches actually stored in the buffer.
     */
    int GetBatch_Qso(void);




    /*!
     * \brief Retrieves a batch of Qmo
     *
     * The batches are stored in the matrix set by SetOutputBuffer().
     * See \ref theory_page for what Qmo actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nmo*nmo elements.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \return The number of batches actually stored in the buffer.
     */
    int GetBatch_Qmo(void);



    /*!
     * \brief Retrieves a batch of Qov (occupied-virtual)
     *
     * The batches are stored in the matrix set by SetOutputBuffer().
     * See \ref theory_page for what Qov actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nocc*nvir elements.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \return The number of batches actually stored in the buffer.
     */
    int GetBatch_Qov(void);


    /*!
     * \brief Retrieves a batch of Qoo (occupied-occupied)
     *
     * The batches are stored in the matrix set by SetOutputBuffer().
     * See \ref theory_page for what Qov actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nocc*nocc elements.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \return The number of batches actually stored in the buffer.
     */
    int GetBatch_Qoo(void);


    /*!
     * \brief Retrieves a batch of Qvv (occupied-virtual)
     *
     * The batches are stored in the matrix set by SetOutputBuffer().
     * See \ref theory_page for what Qvv actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nvir*nvir elements.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \return The number of batches actually stored in the buffer.
     */
    int GetBatch_Qvv(void);



    /*!
     * \brief Resets information about batches
     *
     * Once called, the various GetBatch functions will start from the beginning of Qso
     */
    void ResetBatches(void);




    /*!
     * \brief Sets the maximum number of (OpenMP) threads used
     *
     * Set to zero to use the maximum number of threads for this machine.
     *
     * \param [in] nthread Max number of threads to use
     * \return The max number of threads that will actually be used (ie if \p nthread is zero).
     */
    int SetNThread(int nthread);

/*
    int CalculateERI(double * qso, int qsosize, int shell1, int shell2, int shell3, int shell4, double * outbuffer, int buffersize);

    int CalculateERIMulti(double * qso, int qsosize,
                          int shell1, int nshell1,
                          int shell2, int nshell2,
                          int shell3, int nshell3,
                          int shell4, int nshell4,
                          double * outbuffer, int buffersize);

    void ReorderQ(double * qso, int qsosize, const reorder::Orderings & order);
    void ReorderQ_GAMESS(double * qso, int qsosize);
*/
};

}

#endif //PANACHE_DFTENSOR_H
