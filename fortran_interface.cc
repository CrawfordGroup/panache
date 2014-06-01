#include <iostream>
#include <cstring> // strncpy

#include "c_interface.h"
#include "Output.h"

const std::string Z_to_sym [] = {
//0-10
"XXX ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"
};



extern "C" {


    void fortran_init_(INTTYPE * ncenters,
                       DBLTYPE * xyz, char * symbols,
                       INTTYPE * primary_nshellspercenter, INTTYPE * primary_am, INTTYPE * primary_is_pure,
                       INTTYPE * primary_nprimpershell, DBLTYPE * primary_exp, DBLTYPE * primary_coef,
                       INTTYPE * aux_nshellspercenter, INTTYPE * aux_am, INTTYPE * aux_is_pure,
                       INTTYPE * aux_nprimpershell, DBLTYPE * aux_exp, DBLTYPE * aux_coef, INTTYPE * dfhandle)
    {
        // Make molecule struct
        C_AtomCenter * atoms = new C_AtomCenter[*ncenters];

        char ** csymbols = new char*[*ncenters];
        for(INTTYPE i = 0; i < *ncenters; i++)
        {
            csymbols[i] = new char[5];
            csymbols[i][0] = csymbols[i][1] = csymbols[i][2] 
                           = csymbols[i][3] = csymbols[i][4] = '\0';
        }

        // symbols should be basically a long, space-delimited string
        // length is ncenters * 4, but is NOT null delimited.
        // fortran inserts spaces instead.

        for(INTTYPE i = 0; i < *ncenters; i++)
        {
            for(INTTYPE j = 0; j < 4; j++)
            {
                char c = symbols[i * 4 + j];
                if(c != ' ')
                    csymbols[i][j] = c;
            }

            atoms[i].symbol = csymbols[i]; 

            atoms[i].center[0] = xyz[i];
            atoms[i].center[1] = xyz[i+(*ncenters)];
            atoms[i].center[2] = xyz[i+2*(*ncenters)];
        }

        INTTYPE p_nshell = 0;
        INTTYPE a_nshell = 0;
        for(INTTYPE i = 0; i < *ncenters; i++)
        {
            p_nshell += primary_nshellspercenter[i];
            a_nshell += aux_nshellspercenter[i];
        }

        C_ShellInfo * primary_shells = new C_ShellInfo[p_nshell];
        C_ShellInfo * aux_shells = new C_ShellInfo[a_nshell];

        INTTYPE prim_count = 0;
        for(INTTYPE i = 0; i < p_nshell; i++)
        {
            primary_shells[i].nprim   = primary_nprimpershell[i];
            primary_shells[i].am      = primary_am[i];
            primary_shells[i].ispure  = primary_is_pure[i];
            primary_shells[i].exp     = new DBLTYPE[primary_shells[i].nprim];
            primary_shells[i].coef    = new DBLTYPE[primary_shells[i].nprim];
            for(INTTYPE j = 0; j < primary_shells[i].nprim; j++)
            {
                primary_shells[i].exp[j] = primary_exp[prim_count];
                primary_shells[i].coef[j] = primary_coef[prim_count++];
            }
        }

        prim_count = 0;
        for(INTTYPE i = 0; i < a_nshell; i++)
        {
            aux_shells[i].nprim   = aux_nprimpershell[i];
            aux_shells[i].am      = aux_am[i];
            aux_shells[i].ispure  = aux_is_pure[i];
            aux_shells[i].exp     = new DBLTYPE[aux_shells[i].nprim];
            aux_shells[i].coef    = new DBLTYPE[aux_shells[i].nprim];
            for(INTTYPE j = 0; j < aux_shells[i].nprim; j++)
            {
                aux_shells[i].exp[j] = aux_exp[prim_count];
                aux_shells[i].coef[j] = aux_coef[prim_count++];
            }
        }

        *dfhandle = C_init(*ncenters, atoms,
                           primary_nshellspercenter, primary_shells,
                           aux_nshellspercenter, aux_shells);

        // Free memory
        for(INTTYPE i = 0; i < p_nshell; i++)
        {
            delete [] primary_shells[i].exp;
            delete [] primary_shells[i].coef;
        }

        for(INTTYPE i = 0; i < a_nshell; i++)
        {
            delete [] aux_shells[i].exp;
            delete [] aux_shells[i].coef;
        }

        delete [] primary_shells;
        delete [] aux_shells;

        delete [] atoms;

        for(INTTYPE i = 0; i < *ncenters; i++)
            delete [] csymbols[i];
        delete [] csymbols;


    }


    void fortran_init2_(INTTYPE * ncenters,
                        DBLTYPE * xyz, char * symbols, INTTYPE * symbollen,
                        INTTYPE * primary_nshellspercenter, INTTYPE * primary_am, INTTYPE * primary_is_pure,
                        INTTYPE * primary_nprimpershell, DBLTYPE * primary_exp, DBLTYPE * primary_coef,
                        const char * auxfilename, INTTYPE * auxfilenamelen, INTTYPE * dfhandle)
    {
        // Make molecule struct
        C_AtomCenter * atoms = new C_AtomCenter[*ncenters];

        char ** csymbols = new char*[*ncenters];
        for(INTTYPE i = 0; i < *ncenters; i++)
        {
            csymbols[i] = new char[*symbollen+1];
            for(INTTYPE j = 0; j < ((*symbollen)+1); j++)
                csymbols[i][j] = '\0';
        }

        // symbols should be basically a long, space-delimited string
        // length is ncenters * 4, but is NOT null delimited.
        // fortran inserts spaces instead.
        // So the extra element from above guarentee room for a null
        // character

        for(INTTYPE i = 0; i < *ncenters; i++)
        {
            for(INTTYPE j = 0; j < *symbollen; j++)
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

        INTTYPE p_nshell = 0;
        for(INTTYPE i = 0; i < *ncenters; i++)
        {
            p_nshell += primary_nshellspercenter[i];
            std::cout << "NSHELLPERCENTER " << i << " = " << primary_nshellspercenter[i] << "\n";
        }

        C_ShellInfo * primary_shells = new C_ShellInfo[p_nshell];

        INTTYPE prim_count = 0;
        for(INTTYPE i = 0; i < p_nshell; i++)
        {
            primary_shells[i].nprim   = primary_nprimpershell[i];
            primary_shells[i].am      = primary_am[i];
            primary_shells[i].ispure  = primary_is_pure[i];
            primary_shells[i].exp     = new DBLTYPE[primary_shells[i].nprim];
            primary_shells[i].coef    = new DBLTYPE[primary_shells[i].nprim];
            for(INTTYPE j = 0; j < primary_shells[i].nprim; j++)
            {
                primary_shells[i].exp[j] = primary_exp[prim_count];
                primary_shells[i].coef[j] = primary_coef[prim_count++];
            }
        }

        std::cout << "SIZEOF INT: " << sizeof(INTTYPE) << "\n";
        std::cout << "NCENTER: " << *ncenters << "\n";
        std::cout << "PRIMCOUNT: " << prim_count << "\n";
        std::cout << "P_NSHELL: " << p_nshell << "\n";

        // create a real, null-terminated c string
        char * cfname = new char[*auxfilenamelen+1];
        strncpy(cfname, auxfilename, *auxfilenamelen);
        cfname[*auxfilenamelen] = '\0';

        *dfhandle = C_init2(*ncenters, atoms,
                            primary_nshellspercenter, primary_shells,
                            cfname);


        // Free memory
        delete [] cfname;
        for(INTTYPE i = 0; i < p_nshell; i++)
        {
            delete [] primary_shells[i].exp;
            delete [] primary_shells[i].coef;
        }

        delete [] primary_shells;

        delete [] atoms;

        for(INTTYPE i = 0; i < *ncenters; i++)
            delete [] csymbols[i];
        delete [] csymbols;


    }

    void fortran_init3_(INTTYPE * ncenters,
                        DBLTYPE * xyz, DBLTYPE * Z,
                        INTTYPE * primary_nshellspercenter, INTTYPE * primary_am, INTTYPE * primary_is_pure,
                        INTTYPE * primary_nprimpershell, DBLTYPE * primary_exp, DBLTYPE * primary_coef,
                        const char * auxfilename, INTTYPE * auxfilenamelen, INTTYPE * dfhandle)
    {
        std::cout << "HEREHERE\n";
        // Make molecule struct
        C_AtomCenter * atoms = new C_AtomCenter[*ncenters];

        for(INTTYPE i = 0; i < *ncenters; i++)
        {
            atoms[i].symbol = Z_to_sym[static_cast<INTTYPE>(Z[i])].c_str(); 

            atoms[i].center[0] = xyz[i];
            atoms[i].center[1] = xyz[i+(*ncenters)];
            atoms[i].center[2] = xyz[i+2*(*ncenters)];
        }

        INTTYPE p_nshell = 0;
        for(INTTYPE i = 0; i < *ncenters; i++)
        {
            p_nshell += primary_nshellspercenter[i];
        }

        C_ShellInfo * primary_shells = new C_ShellInfo[p_nshell];

        INTTYPE prim_count = 0;
        for(INTTYPE i = 0; i < p_nshell; i++)
        {
            primary_shells[i].nprim   = primary_nprimpershell[i];
            primary_shells[i].am      = primary_am[i];
            primary_shells[i].ispure  = primary_is_pure[i];
            primary_shells[i].exp     = new DBLTYPE[primary_shells[i].nprim];
            primary_shells[i].coef    = new DBLTYPE[primary_shells[i].nprim];
            for(INTTYPE j = 0; j < primary_shells[i].nprim; j++)
            {
                primary_shells[i].exp[j] = primary_exp[prim_count];
                primary_shells[i].coef[j] = primary_coef[prim_count++];
            }
        }


        // create a real, null-terminated c string
        char * cfname = new char[*auxfilenamelen+1];
        strncpy(cfname, auxfilename, *auxfilenamelen);
        cfname[*auxfilenamelen] = '\0';

        *dfhandle = C_init2(*ncenters, atoms,
                            primary_nshellspercenter, primary_shells,
                            cfname);


        // Free memory
        delete [] cfname;
        for(INTTYPE i = 0; i < p_nshell; i++)
        {
            delete [] primary_shells[i].exp;
            delete [] primary_shells[i].coef;
        }

        delete [] primary_shells;

        delete [] atoms;
        
    }

    void fortran_tensordimensions_(INTTYPE * df_handle, INTTYPE * d1, INTTYPE * d2, INTTYPE * d3, INTTYPE * matsize)
    {
        *matsize = C_TensorDimensions(*df_handle, d1, d2, d3);
    }


    void fortran_qao_(INTTYPE * df_handle, DBLTYPE * matout, INTTYPE * matsize)
    {
        C_QAO(*df_handle, matout, *matsize);
    }

    void fortran_cleanup_(INTTYPE * df_handle)
    {
        C_cleanup(*df_handle);
    }

    void fortran_cleanup_all_(void)
    {
        C_cleanup_all();
    }

    void fortran_eri_(INTTYPE * df_handle, DBLTYPE * qso, INTTYPE * qsosize,
               INTTYPE * shell1, INTTYPE * shell2, INTTYPE * shell3, INTTYPE * shell4,
               DBLTYPE * outbuffer, INTTYPE * buffersize, INTTYPE * ncalc)
    {
        *ncalc = C_CalculateERI(*df_handle, qso, *qsosize, *shell1, *shell2, *shell3, *shell4,
                       outbuffer, *buffersize);
    }

    void fortran_eri_multi_(INTTYPE * df_handle,
                           DBLTYPE * qso, INTTYPE * qsosize,
                           INTTYPE * shell1, INTTYPE * nshell1,
                           INTTYPE * shell2, INTTYPE * nshell2,
                           INTTYPE * shell3, INTTYPE * nshell3,
                           INTTYPE * shell4, INTTYPE * nshell4,
                           DBLTYPE * outbuffer, INTTYPE * buffersize, INTTYPE * ncalc)
    {
        *ncalc = C_CalculateERIMulti(*df_handle, qso, *qsosize,
                                     *shell1, *nshell1, *shell2, *nshell2,
                                     *shell3, *nshell3, *shell4, *nshell4,
                                     outbuffer, *buffersize);
    }

    void fortran_reorderq_gamess_(INTTYPE * df_handle,
                                  DBLTYPE * qso, INTTYPE * qsosize)
    {
        C_ReorderQ_GAMESS(*df_handle, qso, *qsosize);
    }

    void fortran_stdout_(void)
    {
        panache::output::SetOutput(&std::cout);
    }

}

