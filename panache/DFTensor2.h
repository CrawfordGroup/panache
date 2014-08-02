/*! \file
 * \brief Density fitting tensor generation and manipulation (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_DFTENSOR2_H
#define PANACHE_DFTENSOR2_H

#include <fstream>

#include "panache/Timing.h"

#include "Exception.h"

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
    int nfroz_; //!< Number of frozen orbitals
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



    enum class QStorage
    {
        INMEM,
        ONDISK,
        ONFLY
    };



    class StoredQTensor
    {
    private:
        int naux_;
        int ndim1_;
        int ndim2_;
        int ndim12_;
        bool packed_;
        QStorage storetype_;


    protected:
        virtual void Init_(void) = 0;
        virtual void Reset_(void) = 0;
        virtual void Write_(double * data, int ij) = 0;
        virtual void Read_(double * data, int ij) = 0;
        virtual void ReadByQ_(double * data, int nq, int qstart) = 0;
        virtual void Clear_() = 0;

        int naux(void) const { return naux_; }
        int ndim1(void) const { return ndim1_; }
        int ndim2(void) const { return ndim2_; }
        int ndim12(void) const { return ndim12_; }
        int storesize(void) const { return ndim12_*naux_; }
        int packed(void) const { return packed_; }
        int calcindex(int i, int j)
        {
            if(!packed_)
                return (i*ndim2_+j);
            else if(i >= j)
                return ((i*(i+1))>>1) + j;
            else
                return ((j*(j+1))>>1) + i;
        }

    public:
        StoredQTensor(int naux, int ndim1, int ndim2, bool packed, QStorage storetype)
        {
            naux_ = naux;
            ndim1_ = ndim1;
            ndim2_ = ndim2;
            storetype_ = storetype;
            packed_ = packed;

            if(packed && ndim1 != ndim2)
                throw RuntimeError("non square packed matrices?");

            ndim12_ = (packed ? (ndim1_ * (ndim2_+1))/2 : ndim1_*ndim2_);
        }

        virtual ~StoredQTensor()
        {
        }

        QStorage StoreType(void) const
        {
            return storetype_;
        }

        virtual void Write(double * data, int i, int j)
        {
            Write_(data, calcindex(i,j));
        }

        virtual void Read(double * data, int i, int j)
        {
            Read_(data, calcindex(i,j));
        }

        virtual int ReadByQ(double * data, int nq, int qstart)
        {
            if(qstart + nq >= naux_)
                nq = naux_-qstart;

            ReadByQ_(data, nq, qstart);
            return nq;
        }

        virtual void Reset(void)
        {
            Reset_();
        }

        virtual void Clear(void)
        {
            Clear_();
        }

        virtual void Init(void)
        {
            Init_();
        }

    };


/*
    class DiskQTensor : public StoredQTensor
    {
    private:
        string filename_;
        std::unique_ptr<std::fstream> file_;

    protected:
        virtual void Reset_(void)
        {
            if(file_)
            {
                file_->seekg(0);
                file_->seekp(0);
            }
        }

        virtual void Write_(double * data, size_t nq, int ij)
        {
               // \todo write to file 
        }

        virtual void Read_(double * data, size_t nq, int ij)
        {
               // \todo write to file 
        }

    public:
        DiskQTensor(int naux, int ndim1, int ndim2, bool packed, const string & filename)
            : StoredQTensor(naux, ndim1, ndim2, packed, QStorage::ONDISK)
        {
            filename_ = filename;
        }

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


    };
*/

    
    class MemoryQTensor : public StoredQTensor
    {
    private:
        std::unique_ptr<double[]> data_;

    protected:
        virtual void Reset_(void)
        {
            // nothing needed
        }

        virtual void Write_(double * data, int ij)
        {
            double * start = data_.get() + ij * naux();
            std::copy(data, data+naux(), start);
        }

        virtual void Read_(double * data, int ij)
        {
            double * start = data_.get() + ij * naux();
            std::copy(start, start+naux(), data);
        }

        virtual void ReadByQ_(double * data, int nq, int qstart)
        {
            // ugly. Can't really use std::copy
            if(packed())
            {
                for(int q0 = 0, q = qstart; q0 < nq; q++, q0++)
                for(int m = 0; m < ndim1(); m++)
                for(int n = 0; n < ndim2(); n++)
                    data[q0*ndim1()*ndim2() + m*ndim2() + n]
                        = data[q0*ndim1()*ndim2() + n*ndim1() + m]
                        = data_[calcindex(m,n)*naux() + q];
            }
            else
            {
                for(int q0 = 0, q = qstart; q0 < nq; q++, q0++)
                for(int m = 0; m < ndim1(); m++)
                for(int n = 0; n < ndim2(); n++)
                    data[q0 * ndim12() + calcindex(m,n)]
                        = data_[calcindex(m,n)*naux() + q];
            }
        }

        virtual void Clear_(void)
        {
            data_.reset();
        }

        virtual void Init_(void)
        {
            if(!data_)
                data_ = std::unique_ptr<double []>(new double[storesize()]);
        }

    public:
        MemoryQTensor(int naux, int ndim1, int ndim2, bool packed)
            : StoredQTensor(naux, ndim1, ndim2, packed, QStorage::INMEM)
        {
        }

    };

    std::unique_ptr<StoredQTensor> qso_;
    std::unique_ptr<StoredQTensor> qmo_;
    std::unique_ptr<StoredQTensor> qoo_;
    std::unique_ptr<StoredQTensor> qsv_;
    std::unique_ptr<StoredQTensor> qvv_;

    ///@}








    /*! \name Memory, Buffers and placeholders for transformations */
    ///@{

    std::unique_ptr<double[]> qc_;  //!< Holds C(T) Q or QC intermediate
    std::unique_ptr<double[]> q_;   //!< Holds a batch of Qso (packed storage)

    ///@}



    /*! \name Timing */
    ///@{
    Timer timer_genqso;        //!< Total time spent in GenQso()
    /*
    Timer timer_getbatch_qso;  //!< Total time spent in GetBatch_Qso()
    Timer timer_getbatch_qmo;  //!< Total time spent in GetBatch_Qmo()
    Timer timer_getbatch_qoo;  //!< Total time spent in GetBatch_Qoo()
    Timer timer_getbatch_qvv;  //!< Total time spent in GetBatch_Qvv()
    Timer timer_getbatch_qov;  //!< Total time spent in GetBatch_Qov()
    */
    ///@}



    /*! \name Threading and Prefetching */
    ///@{

    int nthreads_;  //!< Number of threads to use

    ///@}



    /*!
     * \brief Generates the basic Qso matrix
     *
     * See \ref theory_page for what Qso actually is, and memory_sec for more information
     * about memory.
     *
     * \param [in] inmem If nonzero, store the Qso matrix in memoryof auxiliary basis functions
     */
    void GenQso(QStorage storetype);


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
     * \param [in] directory Full path to a directory to put scratch files
     * \param [in] nthreads Max number of threads to use
     */ 
    DFTensor2(SharedBasisSet primary,
             SharedBasisSet auxiliary,
             const std::string & directory,
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
     * \param [in] nfroz Number of frozen orbitals
     */
    void SetNOcc(int nocc, int nfroz = 0);




    /*!
     * \brief Sets the maximum number of (OpenMP) threads used
     *
     * Set to zero to use the maximum number of threads for this machine.
     *
     * \param [in] nthread Max number of threads to use
     * \return The max number of threads that will actually be used (ie if \p nthread is zero).
     */
    int SetNThread(int nthread);


    void GenQTensors(int qflags);
    int GetBatch_Qso(double * outbuf, int bufsize, int qstart);
};

}

#endif //PANACHE_DFTENSOR_H

