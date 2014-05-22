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



    int C_init(int ncenters,
               C_AtomCenter * atoms,
               int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               int * aux_nshellspercenter, struct C_ShellInfo * aux_shells);

    void C_QAO(int df_handle, double * matout, int matsize);

    void C_cleanup(int df_handle);
    void C_cleanup_all(void);

    int C_TensorDimensions(int df_handle, int * d1, int * d2, int * d3);

    int C_ERI(int df_handle, double * qso, int qsosize, int shell1, int shell2, int shell3, int shell4, double * outbuffer, int buffersize);

}

#endif
