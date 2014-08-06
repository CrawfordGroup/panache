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
     *
     * \todo Rename?
     */
    struct C_AtomCenter
    {
        const char * symbol;   //!< Atomic symbol (used in printing, debugging, and if the auxiliary basis set is from a file)
        double center[3];      //!< x, y, and z Coordinates (in atomic units)
    };




    /*!
     * \brief Initializes a new density-fitting calculation
     *
     * Sets up the basis set information and calculates the metric. It returns a handle that
     * is used to identify this particular calculation.
     *
     * Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
     *
     * \param [in] ncenters    The number of basis function centers
     * \param [in] atoms       Information about the centers. This is expected to be of length \p ncenters.
     * \param [in] normalized  Are these basis functions normalized or not. Nonzero = No normalization needed.
     * \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
     *                                       Expected to be of length \p ncenters.
     * \param [in] primary_shells  Information about each shell in the primary basis.
     *                             Length should be the sum of \p primary_nshellspercenter.
     * \param [in] aux_nshellspercenter  Number of shells on each center for the auxiliary (density fitting) basis.
     *                                   Expected to be of length \p ncenters.
     * \param [in] aux_shells  Information about each shell in the auxiliary (density fitting) basis.
     *                         Length should be the sum of \p aux_nshellspercenter.
     * \param [in] filename A full path to a file to be used if storing matrices to disk.
     *                      Not referenced if the disk is not used. Should not be set to "NULL", but
     *                      may be set to an empty string if disk is not to be used.
     * \param [in] nthreads Max number of threads to use
     *
     * \return A handle representing this particular density-fitting calculation.
     */
    int_t panache_init(int_t ncenters,
                       C_AtomCenter * atoms, int_t normalized,
                       int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                       int_t * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
                       const char * filename, int_t nthreads);




    /*!
     * \brief Initializes a new density-fitting calculation using an auxiliary basis set file
     *
     * Sets up the basis set information and calculates the metric. It returns a handle that
     * is used to identify this particular calculation.
     *
     * Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
     *
     * \param [in] ncenters    The number of basis function centers
     * \param [in] atoms       Information about the centers. This is expected to be of length \p ncenters.
     * \param [in] normalized  Are these basis functions normalized or not. Nonzero = No normalization needed.
     * \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
     *                                       Expected to be of length \p ncenters.
     * \param [in] primary_shells  Information about each shell in the primary basis.
     *                             Length should be the sum of \p primary_nshellspercenter.
     * \param [in] auxfilename A full path to a file containing the auxiliary basis set (in Gaussian94 format)
     * \param [in] filename A full path to a file to be used if storing matrices to disk.
     *                      Not referenced if the disk is not used. Should not be set to "NULL", but
     *                      may be set to an empty string if disk is not to be used. If used, any existing
     *                      file will be overwritten.
     * \param [in] nthreads Max number of threads to use
     *
     * \return A handle representing this particular density-fitting calculation.
     */
    int_t panache_init2(int_t ncenters,
                        C_AtomCenter * atoms, int_t normalized,
                        int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                        const char * auxfilename, const char * filename, int_t nthreads);





    /*!
     * \brief Sets the C matrix (so-ao matrix) for use in generating Qmo, etc 
     *
     * The matrix is expected be nso x nmo (MOs in the columns) in row-major order.
     * If it is nmo x nso, or the matrix is in column major order, set \p cmo_is_trans.
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
                            int bsorder = BSORDER_PSI4);



    /*!
     * \brief Queries information about the expected matrix dimensions
     *
     * Useful for determining buffer sizes or determining if Qso should be placed in memory.
     * The size of Qso (unpacked) will be naux * nso2, and batches of Qso will be read in
     * multiples of nso2 (for panache_getqbatch_qso()).
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [out] naux Number of auxiliary basis functions
     * \param [out] nso2 Number of primary basis functions squared (nso*nso)
     * \return Total size of an unpacked Qso tensor (ie naux * nso2)
     */
    int_t panache_qsodimensions(int_t df_handle, int_t * naux, int_t * nso2);





    /*!
     * \brief Sets the number of OpenMP threads used in DF routines
     *
     * Set to zero to use the maximum number of threads for this machine.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] nthread Max number of threads to use
     * \return The max number of threads that will actually be used (ie if \p nthread is zero).
     */ 
    int_t panache_setnthread(int_t df_handle, int_t nthread);



    /*!
     * \brief Generates various 3-index tensors
     *
     * For the \p qflags and \p storetype parameters, see Flags.h. For example, to calculate the
     * Qmo and Qov tensors on disk,
     *
     * \code{.c}
     * panache_genqtensors(df_handle, QFLAGS_QMO | QFLAGS_QOV, QSTORAGE_ONDISK);
     * \endcode
     *
     * \note Be sure to set the C-Matrix first!
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] qflags A combination of flags specifying which tensors to generate
     * \param [in] storetype How to store the matrix
     *
     */
    void panache_genqtensors(int_t df_handle, int_t qflags, int_t storetype);



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
     * \brief Sets the number of occupied and virtual orbitals.
     *
     * Number of virtual orbitals is taken to be the remainder after the occupied.
     * Used by Qov, etc.
     *
     * The number of frozen orbitals should be counted in the \p nocc parameter as well.
     *
     * \note You must set the C Matrix first before calling (see SetCMatrix())
     *
     * \param [in] df_handle A handle (returned from an init function) for the DF calculation 
     * \param [in] nocc Number of occupied orbitals. This includes frozen orbitals
     * \param [in] nfroz Number of frozen occupied orbitals
     */
     void panache_setnocc(int_t df_handle, int_t nocc, int_t nfroz);


    /*! \name Retrieving by auxiliary basis index */
    ///@{
    /*!
     * \brief Retrieves a batch of Qso
     *
     * See \ref theory_page for what Qso actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nso2 elements with
     * the index of the auxiliary basis function as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getqbatch_qso(int_t df_handle, double * outbuf, int bufsize, int qstart);



    /*!
     * \brief Retrieves a batch of Qmo
     *
     * See \ref theory_page for what Qmo actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nmo2 elements with
     * the index of the auxiliary basis function as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getqbatch_qmo(int_t df_handle, double * outbuf, int bufsize, int qstart);


    /*!
     * \brief Retrieves a batch of Qoo
     *
     * See \ref theory_page for what Qoo actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nocc*nocc elements with
     * the index of the auxiliary basis function as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getqbatch_qoo(int_t df_handle, double * outbuf, int bufsize, int qstart);

    /*!
     * \brief Retrieves a batch of Qov
     *
     * See \ref theory_page for what Qov actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nocc*nvir elements with
     * the index of the auxiliary basis function as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getqbatch_qov(int_t df_handle, double * outbuf, int bufsize, int qstart);
    ///@}



    /*! \name Retrieving by orbital index */
    ///@{
    /*!
     * \brief Retrieves a batch of Qvv
     *
     * See \ref theory_page for what Qvv actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nvir*nvir elements with
     * the index of the auxiliary basis function as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getqbatch_qvv(int_t df_handle, double * outbuf, int bufsize, int qstart);


    /*!
     * \brief Retrieves a batch of Qso
     *
     * See \ref theory_page for what Qso actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*naux elements with
     * the combined orbital index as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] ijstart The starting value of the ij index
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getbatch_qso(int_t df_handle, double * outbuf, int bufsize, int ijstart);



    /*!
     * \brief Retrieves a batch of Qmo
     *
     * See \ref theory_page for what Qmo actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*naux elements with
     * the combined orbital index as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] ijstart The starting value of the ij index
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getbatch_qmo(int_t df_handle, double * outbuf, int bufsize, int ijstart);


    /*!
     * \brief Retrieves a batch of Qoo
     *
     * See \ref theory_page for what Qoo actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*naux elements with
     * the combined orbital index as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] ijstart The starting value of the ij index
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getbatch_qoo(int_t df_handle, double * outbuf, int bufsize, int ijstart);


    /*!
     * \brief Retrieves a batch of Qov
     *
     * See \ref theory_page for what Qov actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*naux elements with
     * the combined orbital index as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] ijstart The starting value of the ij index
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getbatch_qov(int_t df_handle, double * outbuf, int bufsize, int ijstart);


    /*!
     * \brief Retrieves a batch of Qvv
     *
     * See \ref theory_page for what Qvv actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*naux elements with
     * the combined orbital index as the slowest index.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] ijstart The starting value of the ij index
     * \return The number of batches returned in the \p outbuf buffer.
     */
    int_t panache_getbatch_qvv(int_t df_handle, double * outbuf, int bufsize, int ijstart);
    ///@}

} // end extern "C"


#endif

