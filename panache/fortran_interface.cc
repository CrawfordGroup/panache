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
     * \brief Initializes a new cholesky calculation
     *
     * Sets up the basis set information and returns a handle that
     * is used to identify this particular calculation.
     *
     * Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
     *
     * \note Basis set coefficients should NOT be normalized
     *
     * \param [in] ncenters    The number of basis function centers
     * \param [in] xyz         Coordinates of the basis function centers. In order:
     *                         (x1, y1, z1, x2, y2, z2, ..., xN, yN, zN)
     * \param [in] symbols     Atomic symbols for each center, as a set of \p ncenters strings of length \p symbollen
     * \param [in] symbollen   Length of each symbol in \p symbols
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
     * \param [in] delta     Maximum error in the Cholesky procedure
     * \param [in] directory A full path to a file to be used if storing matrices to disk.
     *                       Not referenced if the disk is not used. Should not be set to "NULL", but
     *                       may be set to an empty string if disk is not to be used.
     * \param [in] directorylen Actual length of \p filename
     * \param [in] bsorder Basis function ordering flag
     * \param [in] nthreads Number of threads to use
     *
     * \param [out] dfhandle A handle representing this particular density-fitting calculation.
     */
    void panachef_chinit_(int_t * ncenters,
                        double * xyz, char * symbols, int_t * symbollen,
                        int_t * primary_nshellspercenter, int_t * primary_am, int_t * primary_is_pure,
                        int_t * primary_nprimpershell, double * primary_exp, double * primary_coef,
                        double * delta, const char * directory,
                        int_t * directorylen, int_t * bsorder, int_t * nthreads, int_t * dfhandle)
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

            strncpy(atoms[i].symbol, csymbols[i], 4); 
            atoms[i].symbol[4] = '\0';

            atoms[i].center[0] = xyz[3*i];
            atoms[i].center[1] = xyz[3*i+1];
            atoms[i].center[2] = xyz[3*i+2];
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
        char * mdir = new char[*directorylen+1];
        strncpy(mdir, directory, *directorylen);
        mdir[*directorylen] = '\0';

        *dfhandle = panache_chinit(*ncenters, atoms,
                                  primary_nshellspercenter, primary_shells,
                                  *delta, mdir, *bsorder, *nthreads);

        // Free memory
        delete [] mdir;

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
}

