#ifndef PANACHE_C_INTERFACE_H
#define PANACHE_C_INTERFACE_H

extern "C" {

    struct C_ShellInfo
    {
        int nprim;
        int am;
        int ispure;
        double * exp; // of length nprim
        double * coef; // of length nprim
    };

    struct C_AtomCenter
    {
        double Z;
        const char * symbol;
        double mass;
        double center[3];
    };

    double * C_QAO(int ncenters,
                 C_AtomCenter * atoms,
                 int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                 int * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
                 int * nrow_out, int * ncol_out);

    void free_matrix(double * mat);

}

#endif
