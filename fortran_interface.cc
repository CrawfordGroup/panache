#include <iostream>

#include "c_interface.h"



extern "C" {


    void fortran_qao_(int * ncenters,
                      double * xyz, double * Z, double * masses,
                      int * primary_nshellspercenter, int * primary_am, int * primary_is_pure,
                      int * primary_nprimpershell, double * primary_exp, double * primary_coef,
                      int * aux_nshellspercenter, int * aux_am, int * aux_is_pure,
                      int * aux_nprimpershell, double * aux_exp, double * aux_coef,
                      double * matout, int * matsize)
    {
        std::cout << "In fortran_qao: " << *ncenters << "\n";

        // Make molecule struct
        C_AtomCenter * atoms = new C_AtomCenter[*ncenters];
        for(int i = 0; i < *ncenters; i++)
        {
            atoms[i].Z = Z[i];
            atoms[i].symbol = "XX";
            atoms[i].mass = masses[i];
            atoms[i].center[0] = xyz[i];
            atoms[i].center[1] = xyz[i+(*ncenters)];
            atoms[i].center[2] = xyz[i+2*(*ncenters)];
        }

        int p_nshell = 0;
        int a_nshell = 0;
        for(int i = 0; i < *ncenters; i++)
        {
            p_nshell += primary_nshellspercenter[i];
            a_nshell += aux_nshellspercenter[i];
        }

        std::cout << "NShells: " << p_nshell << " " << a_nshell << "\n";

        C_ShellInfo * primary_shells = new C_ShellInfo[p_nshell];
        C_ShellInfo * aux_shells = new C_ShellInfo[a_nshell];

        int prim_count = 0;
        for(int i = 0; i < p_nshell; i++)
        {
            primary_shells[i].nprim   = primary_nprimpershell[i];
            primary_shells[i].am      = primary_am[i];
            primary_shells[i].ispure  = primary_is_pure[i];
            primary_shells[i].exp     = new double[primary_shells[i].nprim];
            primary_shells[i].coef    = new double[primary_shells[i].nprim];
            for(int j = 0; j < primary_shells[i].nprim; j++)
            {
                primary_shells[i].exp[j] = primary_exp[prim_count];
                primary_shells[i].coef[j] = primary_coef[prim_count++];
            }
        }

        prim_count = 0;
        for(int i = 0; i < a_nshell; i++)
        {
            aux_shells[i].nprim   = aux_nprimpershell[i];
            aux_shells[i].am      = aux_am[i];
            aux_shells[i].ispure  = aux_is_pure[i];
            aux_shells[i].exp     = new double[aux_shells[i].nprim];
            aux_shells[i].coef    = new double[aux_shells[i].nprim];
            for(int j = 0; j < aux_shells[i].nprim; j++)
            {
                aux_shells[i].exp[j] = aux_exp[prim_count];
                aux_shells[i].coef[j] = aux_coef[prim_count++];
            }
        }

        C_QAO(*ncenters, atoms,
              primary_nshellspercenter, primary_shells,
              aux_nshellspercenter, aux_shells, matout, *matsize);

        std::cout << "element 3359: " << matout[3359] << "\n";
        std::cout << "element 36640: " << matout[36640] << "\n";
        std::cout << "element 4074: " << matout[4074] << "\n";



        // Free memory
        for(int i = 0; i < p_nshell; i++)
        {
            delete [] primary_shells[i].exp;
            delete [] primary_shells[i].coef;
        }

        for(int i = 0; i < a_nshell; i++)
        {
            delete [] aux_shells[i].exp;
            delete [] aux_shells[i].coef;
        }

        delete [] primary_shells;
        delete [] aux_shells;

    }


}

