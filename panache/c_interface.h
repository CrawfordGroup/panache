/*! \file
 *  \brief C interface to the PANACHE library
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef PANACHE_C_INTERFACE_H
#define PANACHE_C_INTERFACE_H

#include "int_t.h"


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
     * \param [in] atoms       Information about the centers. This is expected to be of length ncenters.
     * \param [in] normalized  Are these basis functions normalized or not. Nonzero = No normalization needed.
     * \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
     *                                       Expected to be of length ncenters.
     * \param [in] primary_shells  Information about each shell in the primary basis.
     *                             Length should be the sum of primary_nshellspercenter.
     * \param [in] aux_nshellspercenter  Number of shells on each center for the auxiliary (density fitting) basis.
     *                                   Expected to be of length ncenters.
     * \param [in] aux_shells  Information about each shell in the auxiliary (density fitting) basis.
     *                         Length should be the sum of aux_nshellspercenter.
     * \param [in] filename A full path to a file to be used if storing matrices to disk.
     *                      Not referenced if the disk is not used. Should not be set to "NULL", but
     *                      may be set to an empty string if disk is not to be used.
     *
     * \return A handle representing this particular density-fitting calculation.
     */
    int_t panache_init(int_t ncenters,
                       C_AtomCenter * atoms, int_t normalized,  
                       int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                       int_t * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
                       const char * filename);




    /*!
     * \brief Initializes a new density-fitting calculation using an auxiliary basis set file
     *
     * Sets up the basis set information and calculates the metric. It returns a handle that
     * is used to identify this particular calculation.
     *
     * Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
     *
     * \param [in] ncenters    The number of basis function centers
     * \param [in] atoms       Information about the centers. This is expected to be of length ncenters.
     * \param [in] normalized  Are these basis functions normalized or not. Nonzero = No normalization needed.
     * \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
     *                                       Expected to be of length ncenters.
     * \param [in] primary_shells  Information about each shell in the primary basis.
     *                             Length should be the sum of primary_nshellspercenter.
     * \param [in] auxfilename A full path to a file containing the auxiliary basis set (in Gaussian94 format)
     * \param [in] filename A full path to a file to be used if storing matrices to disk.
     *                      Not referenced if the disk is not used. Should not be set to "NULL", but
     *                      may be set to an empty string if disk is not to be used.
     *
     * \return A handle representing this particular density-fitting calculation.
     */
    int_t panache_init2(int_t ncenters,
                        C_AtomCenter * atoms, int_t normalized,
                        int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                        const char * auxfilename, const char * filename);





    /*!
     * \brief Sets the C matrix (so-ao matrix) for use in generating Qmo
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
     */
    void panache_setcmatrix(int_t df_handle, double * cmo, int_t nmo, int_t cmo_is_trans);




    /*! 
     * \brief Queries information about the expected matrix dimensions
     *
     * Useful for determining buffer sizes or determining if Qso should be placed in memory.
     * The size of Qso (unpacked) will be naux * nso2, and batches of Qso will be read in
     * multiples of nso2 (for panache_getbatch_qso()).
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [out] naux Number of auxiliary basis functions
     * \param [out] nso2 Number of primary basis functions squared (nso*nso)
     * \return Total size of an unpacked Qso tensor (ie naux * nso2)
     */
    int_t panache_qsodimensions(int_t df_handle, int_t * naux, int_t * nso2);





    /*!
     * \brief Generates the basic Qso matrix
     *
     * See \ref theory_page for what Qso actually is, and memory_sec for more information
     * about memory. 
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] inmem If nonzero, store the Qso matrix in memoryof auxiliary basis functions
     */
    void panache_genqso(int_t df_handle, int_t inmem);



    /*!
     * \brief Sets the buffer used for storing batches of Qso or Qmo
     *
     * Batches are read in multiples of either nso2 (panache_getbatch_qso(), see panache_qsodimensions()) or
     * nmo*nmo (panache_getbatch_qmo()). How many can fit in the buffer is determined automatically
     * from the matsize parameter. Any 'left over' buffer space is not used.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] buffer A pointer to memory for a buffer (of \p bufsize size)
     * \param [in] bufsize Number of elements in \p buffer (not number of bytes)
     */ 
    void panache_setoutputbuffer(int_t df_handle, double * matout, int_t matsize);




    /*!
     * \brief Retrieves a batch of Qso
     *
     * The batches are stored in the matrix set by panache_setoutputbuffer().
     * See \ref theory_page for what Qso actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nso2 elements (see panache_qsodimensions()).
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \return The number of batches actually stored in the buffer.
     */
    int_t panache_getbatch_qso(int_t df_handle);



    /*!
     * \brief Retrieves a batch of Qmo
     *
     * The batches are stored in the matrix set by panache_setoutputbuffer().
     * See \ref theory_page for what Qmo actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nmo*nmo elements.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \return The number of batches actually stored in the buffer.
     */
    int_t panache_getbatch_qmo(int_t dfhandle);


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




/*
    int_t panache_CalculateERI(int_t df_handle, double * qso, int_t qsosize, int_t shell1, int_t shell2, int_t shell3, int_t shell4, double * outbuffer, int_t buffersize);

    int_t C_CalculateERIMulti(int_t df_handle,
                            double * qso, int_t qsosize,
                            int_t shell1, int_t nshell1,
                            int_t shell2, int_t nshell2,
                            int_t shell3, int_t nshell3,
                            int_t shell4, int_t nshell4,
                            double * outbuffer, int_t buffersize);

    void C_ReorderQ_GAMESS(int_t df_handle, double * qso, int_t qsosize);
*/
}

    
#endif
