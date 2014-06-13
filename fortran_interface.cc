#include <iostream>
#include <cstring> // strncpy

#include "c_interface.h"
#include "Output.h"

const std::string Z_to_sym [] = {
//0-10
"XXX ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"
};



extern "C" {


    void panachef_init_(int_t * ncenters,
                       double * xyz, char * symbols, int_t * symbollen, int_t * normalized, 
                       int_t * primary_nshellspercenter, int_t * primary_am, int_t * primary_is_pure,
                       int_t * primary_nprimpershell, double * primary_exp, double * primary_coef,
                       int_t * aux_nshellspercenter, int_t * aux_am, int_t * aux_is_pure,
                       int_t * aux_nprimpershell, double * aux_exp, double * aux_coef,
                       const char * filename, int_t * filenamelen, int_t * dfhandle)
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

            atoms[i].center[0] = xyz[i];
            atoms[i].center[1] = xyz[i+(*ncenters)];
            atoms[i].center[2] = xyz[i+2*(*ncenters)];
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
                                 aux_nshellspercenter, aux_shells, filename);

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


    void panachef_init2_(int_t * ncenters,
                        double * xyz, char * symbols, int_t * symbollen, int_t * normalized,
                        int_t * primary_nshellspercenter, int_t * primary_am, int_t * primary_is_pure,
                        int_t * primary_nprimpershell, double * primary_exp, double * primary_coef,
                        const char * auxfilename, int_t * auxfilenamelen, const char * matfilename,
                        int_t * matfilenamelen, int_t * dfhandle)
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

            atoms[i].center[0] = xyz[i];
            atoms[i].center[1] = xyz[i+(*ncenters)];
            atoms[i].center[2] = xyz[i+2*(*ncenters)];
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
        cfname[*matfilenamelen] = '\0';

        *dfhandle = panache_init2(*ncenters, atoms, *normalized,
                                  primary_nshellspercenter, primary_shells,
                                  cfname, mfname);

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

    void panachef_qsodimensions_(int_t * df_handle, int_t * naux, int_t * nso2, int_t * matsize)
    {
        *matsize = panache_qsodimensions(*df_handle, naux, nso2);
    }


    void panachef_genqso_(int_t * df_handle, int_t * inmem)
    {
        panache_genqso(*df_handle, *inmem);
    }

    void panachef_setcmatrix_(int_t * df_handle, double * cmo, int_t * nmo, int_t * cmo_is_trans)
    {
        panache_setcmatrix(*df_handle, cmo, *nmo, *cmo_is_trans);
    }

    void panachef_getbatch_qso_(int_t * df_handle, double * matout, int_t * matsize, int_t * nq)
    {
        *nq = panache_getbatch_qso(*df_handle, matout, *matsize);
    }

    void panachef_getbatch_qmo_(int_t * df_handle, double * matout, int_t * matsize, int_t * nq)
    {
        *nq = panache_getbatch_qmo(*df_handle, matout, *matsize);
    }

    void panachef_cleanup_(int_t * df_handle)
    {
        panache_cleanup(*df_handle);
    }

    void panachef_cleanup_all_(void)
    {
        panache_cleanup_all();
    }

/*
    void panachef_eri_(int_t * df_handle, double * qso, int_t * qsosize,
               int_t * shell1, int_t * shell2, int_t * shell3, int_t * shell4,
               double * outbuffer, int_t * buffersize, int_t * ncalc)
    {
        *ncalc = C_CalculateERI(*df_handle, qso, *qsosize, *shell1, *shell2, *shell3, *shell4,
                       outbuffer, *buffersize);
    }

    void panachef_eri_multi_(int_t * df_handle,
                           double * qso, int_t * qsosize,
                           int_t * shell1, int_t * nshell1,
                           int_t * shell2, int_t * nshell2,
                           int_t * shell3, int_t * nshell3,
                           int_t * shell4, int_t * nshell4,
                           double * outbuffer, int_t * buffersize, int_t * ncalc)
    {
        *ncalc = C_CalculateERIMulti(*df_handle, qso, *qsosize,
                                     *shell1, *nshell1, *shell2, *nshell2,
                                     *shell3, *nshell3, *shell4, *nshell4,
                                     outbuffer, *buffersize);
    }

    void panachef_reorderq_gamess_(int_t * df_handle,
                                  double * qso, int_t * qsosize)
    {
        C_ReorderQ_GAMESS(*df_handle, qso, *qsosize);
    }
*/

    void panachef_stdout_(void)
    {
        panache::output::SetOutput(&std::cout);
    }

}

