/*! \file
 *  \brief C interface to the PANACHE library (header)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef PANACHE_C_INTERFACE_H
#define PANACHE_C_INTERFACE_H

#include "panache/int_t.h"
#include "panache/Flags.h"

extern "C" {

    /*!
     * \brief Holds information about a single basis set shell
     */
    struct C_ShellInfo
    {
        int_t nprim;    //!< Number of primitives in this shell
        int_t am;       //!< Angular momentum of the shell
        int_t ispure;   //!< Non-zero if this shell is a pure (spherical) shell
        double * exp;   //!< Exponents of the primitives of this shell (expected to be of length nprim)
        double * coef;  //!< Coefficients of the primitives of this shell (expected to be of length nprim)
    };



    /*!
     * \brief Information about a basis function center
     */
    struct C_AtomCenter
    {
        const char * symbol;   //!< Atomic symbol (used in printing, debugging, and if the auxiliary basis set is from a file)
        double center[3];      //!< x, y, and z Coordinates (in atomic units)
    };




    /*!
     * \brief Initializes a new density-fitting calculation
     *
     * Sets up the basis set information and returns a handle that
     * is used to identify this particular calculation.
     *
     * Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
     *
     * \note Basis set coefficients should NOT be normalized
     *
     * \param [in] ncenters    The number of basis function centers
     * \param [in] atoms       Information about the centers. This is expected to be of length \p ncenters.
     * \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
     *                                       Expected to be of length \p ncenters.
     * \param [in] primary_shells  Information about each shell in the primary basis.
     *                             Length should be the sum of \p primary_nshellspercenter.
     * \param [in] aux_nshellspercenter  Number of shells on each center for the auxiliary (density fitting) basis.
     *                                   Expected to be of length \p ncenters.
     * \param [in] aux_shells  Information about each shell in the auxiliary (density fitting) basis.
     *                         Length should be the sum of \p aux_nshellspercenter.
     * \param [in] directory A full path to a directory to be used for storing matrices to disk.
     *                       Not referenced if the disk is not used. Should not be set to "NULL", but
     *                       may be set to an empty string if disk is not to be used.
     *                       If used, any existing files will be overwritten.
     * \param [in] nthreads Max number of threads to use
     *
     * \return A handle representing this particular density-fitting calculation.
     */
    int_t panache_dfinit(int_t ncenters,
                       C_AtomCenter * atoms,
                       int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                       int_t * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
                       const char * directory, int_t nthreads);




    /*!
     * \brief Initializes a new density-fitting calculation using an auxiliary basis set file
     *
     * Sets up the basis set information and returns a handle that
     * is used to identify this particular calculation.
     *
     * Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
     *
     * \note Basis set coefficients should NOT be normalized
     *
     * \param [in] ncenters    The number of basis function centers
     * \param [in] atoms       Information about the centers. This is expected to be of length \p ncenters.
     * \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
     *                                       Expected to be of length \p ncenters.
     * \param [in] primary_shells  Information about each shell in the primary basis.
     *                             Length should be the sum of \p primary_nshellspercenter.
     * \param [in] auxfilename A full path to a file containing the auxiliary basis set (in Gaussian94 format)
     * \param [in] directory A full path to a file to be used for storing matrices to disk.
     *                       Not referenced if the disk is not used. Should not be set to "NULL", but
     *                       may be set to an empty string if disk is not to be used.
     *                       If used, any existing files will be overwritten.
     * \param [in] nthreads Max number of threads to use
     *
     * \return A handle representing this particular density-fitting calculation.
     */
    int_t panache_dfinit2(int_t ncenters,
                        C_AtomCenter * atoms,
                        int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                        const char * auxfilename, const char * directory, int_t nthreads);


    /*!
     * \brief Initializes a new cholesky calculation using an auxiliary basis set file
     *
     * Sets up the basis set information and returns a handle that
     * is used to identify this particular calculation.
     *
     * Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
     *
     * \note Basis set coefficients should NOT be normalized
     *
     * \param [in] ncenters    The number of basis function centers
     * \param [in] atoms       Information about the centers. This is expected to be of length \p ncenters.
     * \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
     *                                       Expected to be of length \p ncenters.
     * \param [in] primary_shells  Information about each shell in the primary basis.
     *                             Length should be the sum of \p primary_nshellspercenter.
     * \param [in] delta     Maximum error in the Cholesky procedure                            
     * \param [in] directory A full path to a file to be used for storing matrices to disk.
     *                       Not referenced if the disk is not used. Should not be set to "NULL", but
     *                       may be set to an empty string if disk is not to be used.
     *                       If used, any existing files will be overwritten.
     * \param [in] nthreads Max number of threads to use
     *
     * \return A handle representing this particular density-fitting calculation.
     */
    int_t panache_chinit(int_t ncenters,
                        C_AtomCenter * atoms,
                        int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                        double delta, const char * directory, int_t nthreads);



    /*!
     * \brief Clean up a particular density-fitting calculation and free memory
     *
     * You should not attempt to use the handle afterwards
     *
     * \param [in] df_handle A handle (returned from an init function) for the DF
     *                       calculation to be cleaned up
     */
    void panache_cleanup(int_t df_handle);



    /*!
     * \brief Cleans up all density fitting calculations
     *
     * All handles are invalid after this point
     */
    void panache_cleanup_all(void);



    /*!
     * \brief Sets the text output of PANACHE to the specified file pointer
     *
     * No automatic handling, closing, etc, is done by PANACHE.
     *
     * \param [in] out Output file pointer to use.
     *
     */
    void panache_output(FILE * out);


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
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] cmo Pointer to a nso x nmo matrix representing the MO coefficients
     * \param [in] nmo Number of MOs in this C matrix
     * \param [in] cmo_is_trans Set to non-zero if the matrix is the transpose (nmo x nso) or
     *                          is in column-major order.
     * \param [in] bsorder Ordering of the basis function in the c-matrix. See Flags.h
     */
    void panache_setcmatrix(int_t df_handle, double * cmo, int_t nmo, int_t cmo_is_trans,
                            int_t bsorder = BSORDER_PSI4);


    /*!
     * \brief Sets the number of occupied and virtual orbitals.
     *
     * Number of virtual orbitals is taken to be the remainder after the occupied.
     * Used by Qov, etc.
     *
     * Frozen orbitals should not be counted in \p nocc.
     *
     * \note You must set the C Matrix first before calling (see panache_setcmatrix())
     *
     * \param [in] df_handle A handle (returned from an init function) for the DF calculation 
     * \param [in] nocc Number of (non-frozen) occupied orbitals
     * \param [in] nfroz Number of frozen occupied orbitals
     */
     void panache_setnocc(int_t df_handle, int_t nocc, int_t nfroz = 0);


    /*!
     * \brief Sets the maximum number of (OpenMP) threads used
     *
     * Set to zero to use the value of the environment variable OMP_NUM_THREAD (or
     * set by omp_num_threads, or the default for this machine).
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] nthread Max number of threads to use
     * \return The max number of threads that will actually be used (ie if \p nthread is zero).
     */ 
    int_t panache_setnthread(int_t df_handle, int_t nthread);



    /*!
     * \brief Prints out timing information collected so far
     *
     * All times are cumulative for all operations. The output must be set
     *  first (See Output.h)
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     */
    void panache_printtimings(int_t df_handle);



    /*!
     * \brief Generates various 3-index tensors
     *
     * For the \p qflags and \p storeflags parameters, see Flags.h. For example, to calculate the
     * Qmo and Qov tensors on disk,
     *
     * \code{.c}
     * panache_genqtensors(df_handle, QGEN_QMO | QGEN_QOV, QSTORAGE_ONDISK);
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
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] qflags A combination of flags specifying which tensors to generate
     * \param [in] storeflags How to store the matrix
     */
    void panache_genqtensors(int_t df_handle, int_t qflags, int_t storeflags);

    /*!
     * \brief Obtain the batch size for panache_getqbatch()
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
     * \param [in] df_handle A handle (returned from an init function) for the DF calculation 
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \return Size of batches returned by panache_getqbatch
     */
    int_t panache_qbatchsize(int_t df_handle, int_t tensorflag);


    /*!
     * \brief Obtain the batch size for panache_getbatch()
     *
     * The size will always be naux (number of auxiliary basis functions)
     *
     * \param [in] df_handle A handle (returned from an init function) for the DF calculation 
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \return Size of batches returned by panache_getbatch()
     */
    int_t panache_batchsize(int_t df_handle, int_t tensorflag);



    /*!
     * \brief See if a particular tensor is stored packed
     *
     * \param [in] df_handle A handle (returned from an init function) for the DF calculation 
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \return Nonzero if the tensor is in packed storage
     */
    int_t panache_ispacked(int_t df_handle, int_t tensorflag);

    /*!
     * \brief Calculate a combined orbital index
     *
     * Depends on packing
     * 
     * \param [in] df_handle A handle (returned from an init function) for the DF calculation 
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \param [in] i First orbital index
     * \param [in] j Second orbital index
     * \return ij, depending on packing
     */
    int_t panache_calcindex(int_t df_handle, int_t tensorflag, int_t i, int_t j);
    


    /*!
     * \brief Obtain the dimensions of a tensor
     *
     * \param [in] df_handle A handle (returned from an init function) for the DF calculation 
     * \param [in] tensorflag Which tensor to query (see Flags.h)
     * \param [out] naux Number of auxiliary basis functions
     * \param [out] ndim1 First dimension for a particular q
     * \param [out] ndim2 Second dimension for a particular q
     * \return Total tensor size (depends on packing)
     */
    int_t panache_tensordimensions(int_t df_handle, int_t tensorflag,
                                   int_t * naux, int_t * ndim1, int_t * ndim2);






    /*!
     * \brief Retrieves a batch of a 3-index tensor 
     *
     * See \ref theory_page for what these tensors actually are, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*qbatchsize elements with
     * the index of the auxiliary basis function as the slowest index.
     *
     * The batchsize can be obtained using QBatchSize()
     *
     * Call this and process the batches, incrementing qstart by the return value,
     * until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] tensorflag Which tensor to get (see Flags.h)
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \return The number of batches actually stored in the buffer.
     */
    int_t panache_getqbatch(int_t df_handle, int_t tensorflag, double * outbuf, int_t bufsize, int_t qstart);



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
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] tensorflag Which tensor to get (see Flags.h)
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] ijstart The starting value of q
     * \return The number of batches actually stored in the buffer.
     */
    int_t panache_getbatch(int_t df_handle, int_t tensorflag, double * outbuf, int_t bufsize, int_t ijstart);


} // end extern "C"


#endif

