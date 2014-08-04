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
public:
    enum class BSOrder
    {
        Psi4,
        GAMESS
    };

    enum class QStorage
    {
        INMEM,
        ONDISK,
        ONFLY
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


    int GetQBatch_Qso(double * outbuf, int bufsize, int qstart);
    int GetQBatch_Qmo(double * outbuf, int bufsize, int qstart);
    int GetQBatch_Qoo(double * outbuf, int bufsize, int qstart);
    int GetQBatch_Qov(double * outbuf, int bufsize, int qstart);
    int GetQBatch_Qvv(double * outbuf, int bufsize, int qstart);


    int GetBatch_Qso(double * outbuf, int bufsize, int ijstart);
    int GetBatch_Qmo(double * outbuf, int bufsize, int ijstart);
    int GetBatch_Qoo(double * outbuf, int bufsize, int ijstart);
    int GetBatch_Qov(double * outbuf, int bufsize, int ijstart);
    int GetBatch_Qvv(double * outbuf, int bufsize, int ijstart);


    void GenQTensors(int qflags, QStorage storetype);



    void PrintTimings(void) const;


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
     * \brief Reorders the whole C matrix into the specified ordering
     * 
     * \param [in] order The order to use
     */
    void ReorderCMat(reorder::Orderings & order);
    ///@}





    class StoredQTensor
    {
    private:
        int naux_;
        int ndim1_;
        int ndim2_;
        int ndim12_;
        bool packed_;
        bool byq_;
        QStorage storetype_;
        Timer gen_timer_;
        Timer getijbatch_timer_;
        Timer getqbatch_timer_;

    protected:
        virtual void Init_(void) = 0;
        virtual void Reset_(void) = 0;
        virtual void Write_(double * data, int nij, int ijstart) = 0;
        virtual void WriteByQ_(double * data, int nq, int qstart) = 0;
        virtual void Read_(double * data, int nij, int ijstart) = 0;
        virtual void ReadByQ_(double * data, int nq, int qstart) = 0;
        virtual void Clear_() = 0;

        int storesize(void) const;
        int byq(void) const;

    public:
        StoredQTensor(int naux, int ndim1, int ndim2, bool packed, bool byq, QStorage storetype);
        virtual ~StoredQTensor();
        QStorage StoreType(void) const;
        void Write(double * data, int nij, int ijstart);
        void WriteByQ(double * data, int nq, int qstart);
        int Read(double * data, int nij, int ijstart);
        int ReadByQ(double * data, int nq, int qstart);
        void Reset(void);
        void Clear(void);
        void Init(void);

        Timer & GenTimer(void);
        Timer & GetBatchTimer(void);
        Timer & GetQBatchTimer(void);

        int naux(void) const;
        int ndim1(void) const;
        int ndim2(void) const;
        int ndim12(void) const;
        int packed(void) const;
        int calcindex(int i, int j) const;
    };


    class DiskQTensor : public StoredQTensor
    {
    private:
        std::unique_ptr<std::fstream> file_;
        std::string filename_;

        void OpenFile_(void);
        void CloseFile_(void);

    protected:
        virtual void Reset_(void);
        virtual void Write_(double * data, int nij, int ijstart);
        virtual void WriteByQ_(double * data, int nq, int qstart);
        virtual void Read_(double * data, int nij, int ijstart);
        virtual void ReadByQ_(double * data, int nq, int qstart);
        virtual void Clear_(void);
        virtual void Init_(void);

    public:
        DiskQTensor(int naux, int ndim1, int ndim2, bool packed, bool byq, const std::string & filename);

    };

    
    class MemoryQTensor : public StoredQTensor
    {
    private:
        std::unique_ptr<double[]> data_;

    protected:
        virtual void Reset_(void);
        virtual void Write_(double * data, int nij, int ijstart);
        virtual void WriteByQ_(double * data, int nq, int qstart);
        virtual void Read_(double * data, int nij, int ijstart);
        virtual void ReadByQ_(double * data, int nq, int qstart);
        virtual void Clear_(void);
        virtual void Init_(void);

    public:
        MemoryQTensor(int naux, int ndim1, int ndim2, bool packed, bool byq);

    };


    std::unique_ptr<StoredQTensor> StoredQTensorFactory(int naux, int ndim1, int ndim2,
                                                        bool packed, bool byq, QStorage storetype, const std::string & name);
    std::unique_ptr<StoredQTensor> qso_;
    std::unique_ptr<StoredQTensor> qmo_;
    std::unique_ptr<StoredQTensor> qoo_;
    std::unique_ptr<StoredQTensor> qov_;
    std::unique_ptr<StoredQTensor> qvv_;

    ///@}








    ///@}


    void TransformQTensor(double * left, int lncols,
                          double * right, int rncols,
                          std::unique_ptr<StoredQTensor> & qout,
                          int q,
                          double * qso, double * qc, double * cqc);

    /*! \name Timing */
    ///@{
    Timer timer_genqtensors_;   //!< Total time spent in GenQTensors()
    ///@}



    /*! \name Threading and Prefetching */
    ///@{

    int nthreads_;  //!< Number of threads to use

    ///@}

    std::string directory_; 


    /*!
     * \brief Generates the basic Qso matrix
     *
     * See \ref theory_page for what Qso actually is, and memory_sec for more information
     * about memory.
     *
     * \param [in] inmem If nonzero, store the Qso matrix in memoryof auxiliary basis functions
     */
    void GenQso(QStorage storetype);


    int GetQBatch_Base(double * outbuf, int bufsize, int qstart, StoredQTensor * qt);
    int GetBatch_Base(double * outbuf, int bufsize, int ijstart, StoredQTensor * qt);

    // helper
    void PrintTimer(const char * name, const std::unique_ptr<StoredQTensor> & q) const;

};

}

#endif //PANACHE_DFTENSOR_H

