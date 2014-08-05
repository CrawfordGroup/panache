/*! \file
 *  \brief Fortran interface to the PANACHE library
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */


#include <iostream>
#include <cstring> // strncpy

#include "panache/c_interface.h"
#include "panache/Output.h"


extern "C" {


    /*!
     * \brief Initializes a new density-fitting calculation
     *
     * Sets up the basis set information and calculates the metric. It returns a handle that
     * is used to identify this particular calculation.
     *
     * Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
     *
     * \param [in] ncenters    The number of basis function centers
     * \param [in] xyz         Coordinates of the basis function centers. In order:
     *                         (x1, y1, z1, x2, y2, z2, ..., xN, yN, zN)
     * \param [in] symbols     Atomic symbols for each center, as a set of \p ncenters strings of length \p symbollen
     * \param [in] symbollen   Length of each symbol in \p symbols
     * \param [in] normalized  Are these basis functions normalized or not. Nonzero = No normalization needed.
     *
     * \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
     *                                       Expected to be of length ncenters.
     * \param [in] primary_am  Angular momentum of each shell (s = 0, p = 1, etc) in the primary basis. 
     *                         Length should be the sum of primary_nshellspercenter.
     * \param [in] primary_is_pure  Whether each shell is pure/spherical or not (primary basis).
     *                              Length should be the sum of primary_nshellspercenter.
     * \param [in] primary_nprimpershell  Number of primitives in each shell of the primary basis. 
     *                                    Length should be the sum of primary_nshellspercenter.
     * \param [in] primary_exp  All exponents for all shells of the primary basis. 
     *                          Length should be the sum of primary_nprimpershell, with grouping
     *                          by shell.
     * \param [in] primary_coef All basis function coefficients for all shells of the primary basis. 
     *                          Length should be the sum of primary_nprimpershell, with grouping
     *                          by shell.
     *
     * \param [in] aux_nshellspercenter  Number of shells on each center for the auxiliary basis.
     *                                   Expected to be of length ncenters.
     * \param [in] aux_am  Angular momentum of each shell (s = 0, p = 1, etc) in the auxiliary basis. 
     *                     Length should be the sum of aux_nshellspercenter.
     * \param [in] aux_is_pure  Whether each shell is pure/spherical or not (auxiliary basis).
     *                          Length should be the sum of aux_nshellspercenter.
     * \param [in] aux_nprimpershell  Number of primitives in each shell of the auxiliary basis. 
     *                                Length should be the sum of aux_nshellspercenter.
     * \param [in] aux_exp  All exponents for all shells of the auxiliary basis. 
     *                      Length should be the sum of aux_nprimpershell, with grouping
     *                      by shell.
     * \param [in] aux_coef All basis function coefficients for all shells of the auxiliary basis. 
     *                      Length should be the sum of aux_nprimpershell, with grouping
     *                      by shell.
     * \param [in] filename A full path to a file to be used if storing matrices to disk.
     *                      Not referenced if the disk is not used. Should not be set to "NULL", but
     *                      may be set to an empty string if disk is not to be used.
     * \param [in] filenamelen Actual length of \p filename
     * \param [in] nthreads Number of threads to use
     *
     * \param [out] dfhandle A handle representing this particular density-fitting calculation.
     */
    void panachef_init_(int_t * ncenters,
                       double * xyz, char * symbols, int_t * symbollen, int_t * normalized, 
                       int_t * primary_nshellspercenter, int_t * primary_am, int_t * primary_is_pure,
                       int_t * primary_nprimpershell, double * primary_exp, double * primary_coef,
                       int_t * aux_nshellspercenter, int_t * aux_am, int_t * aux_is_pure,
                       int_t * aux_nprimpershell, double * aux_exp, double * aux_coef,
                       const char * filename, int_t * filenamelen, int_t * nthreads, int_t * dfhandle)
    {
        // Make molecule struct
        C_AtomCenter * atoms = new C_AtomCenter[*ncenters];

        char ** csymbols = new char*[*ncenters];
        for(int_t i = 0; i < *ncenters; i++)
        {
            csymbols[i] = new char[*symbollen+1];
            for(int_t j = 0; j < ((*symbollen)+1); j++)
                csymbols[i][j] = '\0';
        }


        // symbols should be basically a long, space-delimited string
        // length is ncenters * 4, but is NOT null delimited.
        // fortran inserts spaces instead.

        for(int_t i = 0; i < *ncenters; i++)
        {
            for(int_t j = 0; j < *symbollen; j++)
            {
                char c = symbols[i * (*symbollen) + j];
                if(c != ' ')
                    csymbols[i][j] = c;
            }

            atoms[i].symbol = csymbols[i]; 

            atoms[i].center[0] = xyz[i*(*ncenters)];
            atoms[i].center[1] = xyz[i*(*ncenters)+1];
            atoms[i].center[2] = xyz[i*(*ncenters)+2];
        }

        int_t p_nshell = 0;
        int_t a_nshell = 0;
        for(int_t i = 0; i < *ncenters; i++)
        {
            p_nshell += primary_nshellspercenter[i];
            a_nshell += aux_nshellspercenter[i];
        }

        C_ShellInfo * primary_shells = new C_ShellInfo[p_nshell];
        C_ShellInfo * aux_shells = new C_ShellInfo[a_nshell];

        int_t prim_count = 0;
        for(int_t i = 0; i < p_nshell; i++)
        {
            primary_shells[i].nprim   = primary_nprimpershell[i];
            primary_shells[i].am      = primary_am[i];
            primary_shells[i].ispure  = primary_is_pure[i];
            primary_shells[i].exp     = new double[primary_shells[i].nprim];
            primary_shells[i].coef    = new double[primary_shells[i].nprim];
            for(int_t j = 0; j < primary_shells[i].nprim; j++)
            {
                primary_shells[i].exp[j] = primary_exp[prim_count];
                primary_shells[i].coef[j] = primary_coef[prim_count++];
            }
        }

        prim_count = 0;
        for(int_t i = 0; i < a_nshell; i++)
        {
            aux_shells[i].nprim   = aux_nprimpershell[i];
            aux_shells[i].am      = aux_am[i];
            aux_shells[i].ispure  = aux_is_pure[i];
            aux_shells[i].exp     = new double[aux_shells[i].nprim];
            aux_shells[i].coef    = new double[aux_shells[i].nprim];
            for(int_t j = 0; j < aux_shells[i].nprim; j++)
            {
                aux_shells[i].exp[j] = aux_exp[prim_count];
                aux_shells[i].coef[j] = aux_coef[prim_count++];
            }
        }

        *dfhandle = panache_init(*ncenters, atoms, *normalized,
                                 primary_nshellspercenter, primary_shells,
                                 aux_nshellspercenter, aux_shells, filename, *nthreads);

        // Free memory
        for(int_t i = 0; i < p_nshell; i++)
        {
            delete [] primary_shells[i].exp;
            delete [] primary_shells[i].coef;
        }

        for(int_t i = 0; i < a_nshell; i++)
        {
            delete [] aux_shells[i].exp;
            delete [] aux_shells[i].coef;
        }

        delete [] primary_shells;
        delete [] aux_shells;

        delete [] atoms;

        for(int_t i = 0; i < *ncenters; i++)
            delete [] csymbols[i];
        delete [] csymbols;


    }


    /*!
     * \brief Initializes a new density-fitting calculation using an auxiliary basis set file
     *
     * Sets up the basis set information and calculates the metric. It returns a handle that
     * is used to identify this particular calculation.
     *
     * Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
     *
     * \param [in] ncenters    The number of basis function centers
     * \param [in] xyz         Coordinates of the basis function centers. In order:
     *                         (x1, y1, z1, x2, y2, z2, ..., xN, yN, zN)
     * \param [in] symbols     Atomic symbols for each center, as a set of \p ncenters strings of length \p symbollen
     * \param [in] symbollen   Length of each symbol in \p symbols
     * \param [in] normalized  Are these basis functions normalized or not. Nonzero = No normalization needed.
     *
     * \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
     *                                       Expected to be of length ncenters.
     * \param [in] primary_am  Angular momentum of each shell (s = 0, p = 1, etc) in the primary basis. 
     *                         Length should be the sum of primary_nshellspercenter.
     * \param [in] primary_is_pure  Whether each shell is pure/spherical or not (primary basis).
     *                              Length should be the sum of primary_nshellspercenter.
     * \param [in] primary_nprimpershell  Number of primitives in each shell of the primary basis. 
     *                                    Length should be the sum of primary_nshellspercenter.
     * \param [in] primary_exp  All exponents for all shells of the primary basis. 
     *                          Length should be the sum of primary_nprimpershell, with grouping
     *                          by shell.
     * \param [in] primary_coef All basis function coefficients for all shells of the primary basis. 
     *                          Length should be the sum of primary_nprimpershell, with grouping
     *                          by shell.
     * \param [in] auxfilename A full path to a file containing the auxiliary basis set (in Gaussian94 format)
     * \param [in] auxfilenamelen Actual length of \p auxfilename
     *
     * \param [in] matfilename A full path to a file to be used if storing matrices to disk.
     *                      Not referenced if the disk is not used. Should not be set to "NULL", but
     *                      may be set to an empty string if disk is not to be used.
     * \param [in] matfilenamelen Actual length of \p filename
     * \param [in] nthreads Number of threads to use
     *
     * \param [out] dfhandle A handle representing this particular density-fitting calculation.
     */
    void panachef_init2_(int_t * ncenters,
                        double * xyz, char * symbols, int_t * symbollen, int_t * normalized,
                        int_t * primary_nshellspercenter, int_t * primary_am, int_t * primary_is_pure,
                        int_t * primary_nprimpershell, double * primary_exp, double * primary_coef,
                        const char * auxfilename, int_t * auxfilenamelen, const char * matfilename,
                        int_t * matfilenamelen, int_t * nthreads, int_t * dfhandle)
    {
        // DUMP
        /*
        std::cout << "ncenters: " << *ncenters << "\n";
        std::cout << "XYZ: \n";
        int cc = 0;
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < (*ncenters); j++)
                std::cout << xyz[cc++] << "   ";
            std::cout << "\n";
        }
        std::cout << "symbol len: " << *symbollen << "\n";
        std::cout << "Symbols: \n";
        cc = 0;
        for(int i = 0; i < *ncenters; i++)
        {
            std::cout << "\"";
            for(int j = 0; j < *symbollen; j++)
                std::cout << symbols[cc++];
            std::cout << "\"\n";
        }
        */

        // Make molecule struct
        C_AtomCenter * atoms = new C_AtomCenter[*ncenters];

        char ** csymbols = new char*[*ncenters];
        for(int_t i = 0; i < *ncenters; i++)
        {
            csymbols[i] = new char[*symbollen+1];
            for(int_t j = 0; j < ((*symbollen)+1); j++)
                csymbols[i][j] = '\0';
        }

        // symbols should be basically a long, space-delimited string
        // length is ncenters * 4, but is NOT null delimited.
        // fortran inserts spaces instead.
        // So the extra element from above guarentee room for a null
        // character

        for(int_t i = 0; i < *ncenters; i++)
        {
            for(int_t j = 0; j < *symbollen; j++)
            {
                char c = symbols[i * (*symbollen) + j];
                if(c != ' ')
                    csymbols[i][j] = c;
            }

            atoms[i].symbol = csymbols[i]; 

            atoms[i].center[0] = xyz[i*(*ncenters)];
            atoms[i].center[1] = xyz[i*(*ncenters)+1];
            atoms[i].center[2] = xyz[i*(*ncenters)+2];
            //std::cout << "Adding atom " << i << "(" << atoms[i].symbol << ") with center:\n";
            //std::cout << atoms[i].center[0] << "  " << atoms[i].center[1] << "  " << atoms[i].center[2] << "\n";
        }

        int_t p_nshell = 0;
        for(int_t i = 0; i < *ncenters; i++)
        {
            p_nshell += primary_nshellspercenter[i];
            //std::cout << "NSHELLPERCENTER " << i << " = " << primary_nshellspercenter[i] << "\n";
        }

        C_ShellInfo * primary_shells = new C_ShellInfo[p_nshell];

        int_t prim_count = 0;
        /*
        std::cout << "P_NSHELL: " << p_nshell << "\n";
        std::cout << "NPRIMPERSHELL: \n";
        for(int i = 0; i < p_nshell; i++)
            std::cout << primary_nprimpershell[i] << "  ";
        std::cout << "\n";
        */

        for(int_t i = 0; i < p_nshell; i++)
        {
            //std::cout << "shell " << i << " nprim " << primary_nprimpershell[i]
            //          << " am " << primary_am[i] << " pure? " << primary_is_pure[i] << "\n";
            primary_shells[i].nprim   = primary_nprimpershell[i];
            primary_shells[i].am      = primary_am[i];
            primary_shells[i].ispure  = primary_is_pure[i];
            primary_shells[i].exp     = new double[primary_shells[i].nprim];
            primary_shells[i].coef    = new double[primary_shells[i].nprim];
            for(int_t j = 0; j < primary_shells[i].nprim; j++)
            {
                primary_shells[i].exp[j] = primary_exp[prim_count];
                primary_shells[i].coef[j] = primary_coef[prim_count++];
            }
        }

        /*
        std::cout << "PRIMCOUNT: " << prim_count << "\n";
        std::cout << "exp, coef\n";
        for(int i = 0; i < prim_count; i++)
            std::cout << primary_exp[i] << "  " << primary_coef[i] << "\n";
        std::cout << "\n";
        */

        // create a real, null-terminated c string
        char * cfname = new char[*auxfilenamelen+1];
        char * mfname = new char[*matfilenamelen+1];
        strncpy(cfname, auxfilename, *auxfilenamelen);
        strncpy(mfname, matfilename, *matfilenamelen);
        cfname[*auxfilenamelen] = '\0';
        mfname[*matfilenamelen] = '\0';

        *dfhandle = panache_init2(*ncenters, atoms, *normalized,
                                  primary_nshellspercenter, primary_shells,
                                  cfname, mfname, *nthreads);

        // Free memory
        delete [] cfname;
        delete [] mfname;

        for(int_t i = 0; i < p_nshell; i++)
        {
            delete [] primary_shells[i].exp;
            delete [] primary_shells[i].coef;
        }

        delete [] primary_shells;

        delete [] atoms;

        for(int_t i = 0; i < *ncenters; i++)
            delete [] csymbols[i];
        delete [] csymbols;


    }


    
    /*!
     * \brief Sets the C matrix (so-ao matrix) for use in generating Qmo, etc
     *
     * The matrix is expected be nso x nmo (MOs in the columns) in column-major order.
     * If it is nmo x nso, or the matrix is in row major order, set \p cmo_is_trans.
     *
     * \note This is different from panache_setcmatrix(), as fortran column-major order
     * matrices are handled automatically.
     *
     * The matrix is copied by the PANACHE code, so it can be safely deleted or otherwise
     * changed after calling this function.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] cmo Pointer to a nso x nmo matrix representing the MO coefficients
     * \param [in] nmo Number of MOs in this C matrix
     * \param [in] cmo_is_trans Set to non-zero if the matrix is the transpose (nmo x nso) or
     *                          is in row-major order.
     * \param [in] bsorder Ordering of the C matrix passed in (see BSORDER_* in Flags.h)
     */
    void panachef_setcmatrix_(int_t * df_handle, double * cmo, int_t * nmo, int_t * cmo_is_trans,
                              int_t * bsorder)
    {
        if(*cmo_is_trans)
          panache_setcmatrix(*df_handle, cmo, *nmo, 0, *bsorder);
        else
          panache_setcmatrix(*df_handle, cmo, *nmo, 1, *bsorder);
    }



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
     * \param [out] matsize size of an unpacked Qso tensor (ie naux * nso2)
     */
    void panachef_qsodimensions_(int_t * df_handle, int_t * naux, int_t * nso2, int_t * matsize)
    {
        *matsize = panache_qsodimensions(*df_handle, naux, nso2);
    }



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
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \param [out] nq The number of batches actually stored in the \p outbuf buffer.
     */
    void panachef_getqbatch_qso_(int_t * df_handle, double * outbuf, int_t * bufsize,
                                 int_t * qstart, int_t * nq)
    {
        *nq = panache_getqbatch_qso(*df_handle, outbuf, *bufsize, *qstart);
    }



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
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \param [out] nq The number of batches actually stored in the \p outbuf buffer.
     */
    void panachef_getqbatch_qmo_(int_t * df_handle, double * outbuf, int_t * bufsize,
                                 int_t * qstart, int_t * nq)
    {
        *nq = panache_getqbatch_qmo(*df_handle, outbuf, *bufsize, *qstart);
    }

    /*!
     * \brief Retrieves a batch of Qoo
     *
     * The batches are stored in the matrix set by panache_setoutputbuffer().
     * See \ref theory_page for what Qoo actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nmo*nmo elements.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \param [out] nq The number of batches actually stored in the \p outbuf buffer.
     */
    void panachef_getqbatch_qoo_(int_t * df_handle, double * outbuf, int_t * bufsize,
                                 int_t * qstart, int_t * nq)
    {
        *nq = panache_getqbatch_qoo(*df_handle, outbuf, *bufsize, *qstart);
    }

    /*!
     * \brief Retrieves a batch of Qov
     *
     * The batches are stored in the matrix set by panache_setoutputbuffer().
     * See \ref theory_page for what Qov actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nmo*nmo elements.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \param [out] nq The number of batches actually stored in the \p outbuf buffer.
     */
    void panachef_getqbatch_qov_(int_t * df_handle, double * outbuf, int_t * bufsize,
                                 int_t * qstart, int_t * nq)
    {
        *nq = panache_getqbatch_qov(*df_handle, outbuf, *bufsize, *qstart);
    }

    /*!
     * \brief Retrieves a batch of Qvv
     *
     * The batches are stored in the matrix set by panache_setoutputbuffer().
     * See \ref theory_page for what Qvv actually is, and memory_sec for more information
     * about memory.
     *
     * This function returns the number of batches it has stored in the buffer. The buffer
     * will contain (number of batches)*nmo*nmo elements.
     *
     * Call this and process the batches until this function returns zero.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] outbuf Memory location to store the tensor
     * \param [in] bufsize The size of \p outbuf (in number of doubles)
     * \param [in] qstart The starting value of q
     * \param [out] nq The number of batches actually stored in the \p outbuf buffer.
     */
    void panachef_getqbatch_qvv_(int_t * df_handle, double * outbuf, int_t * bufsize,
                                 int_t * qstart, int_t * nq)
    {
        *nq = panache_getqbatch_qvv(*df_handle, outbuf, *bufsize, *qstart);
    }


    /*!
     * \brief Clean up a particular density-fitting calculation and free memory
     *
     * You should not attempt to use the handle afterwards
     *
     * \param [in] df_handle A handle (returned from an init function) for the DF
     *                       calculation to be cleaned up
     */
    void panachef_cleanup_(int_t * df_handle)
    {
        panache_cleanup(*df_handle);
    }




    /*!
     * \brief Cleans up all density fitting calculations
     *
     * All handles are invalid after this point
     */
    void panachef_cleanup_all_(void)
    {
        panache_cleanup_all();
    }


    /*!
     * \brief Sets the number of OpenMP threads used in DF routines
     *
     * Set to zero to use the maximum number of threads for this machine.
     *
     * \param [in] df_handle A handle (returned from an init function) for this DF calculation
     * \param [in] nthread Max number of threads to use
     * \param [out] actual The max number of threads that will actually be used (ie if \p nthread is zero).
     */ 
    void panachef_setnthread_(int_t * df_handle, int_t * nthread, int_t * actual)
    {
        *actual = panache_setnthread(*df_handle, *nthread);
    }


    /*!
     * \brief Sets the text output of PANACHE to stout
     */ 
    void panachef_stdout_(void)
    {
        panache::output::SetOutput(&std::cout);
    }



    /*!
     * \brief Sets the number of occupied and virtual orbitals.
     *
     * Number of virtual orbitals is taken to be the remainder after the occupied.
     * Used by Qov, etc.
     *
     * \note You must set the C Matrix first before calling (see SetCMatrix())
     *
     * \param [in] df_handle A handle (returned from an init function) for the DF calculation 
     * \param [in] nocc Number of occupied orbitals
     */
     void panachef_setnocc_(int_t * df_handle, int_t * nocc, int_t * nfroz)
     {
        panache_setnocc(*df_handle, *nocc, *nfroz);
     }


}

