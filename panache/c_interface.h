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

    void panache_setcmatrix(int_t df_handle, double * cmo, int_t nmo, int_t cmo_is_trans); 
    void panache_genqso(int_t df_handle, int_t inmem);

    void panache_setoutputbuffer(int_t df_handle, double * matout, int_t matsize);
    int_t panache_getbatch_qso(int_t df_handle);
    int_t panache_getbatch_qmo(int_t dfhandle);

    void panache_cleanup(int_t df_handle);
    void panache_cleanup_all(void);

    int_t panache_qsodimensions(int_t df_handle, int_t * naux, int_t * nso2);

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
