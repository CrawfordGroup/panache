#ifndef PANACHE_C_INTERFACE_H
#define PANACHE_C_INTERFACE_H

#ifndef INTTYPE
    #define INTTYPE int64_t
#endif

#ifndef DBLTYPE
    #define DBLTYPE double
#endif


extern "C" {

    struct C_ShellInfo
    {
        INTTYPE nprim;
        INTTYPE am;
        INTTYPE ispure;
        double * exp; // of length nprim
        double * coef; // of length nprim
    };

    struct C_AtomCenter
    {
        const char * symbol;
        double center[3];
    };



    INTTYPE C_init(INTTYPE ncenters,
               C_AtomCenter * atoms,
               INTTYPE * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               INTTYPE * aux_nshellspercenter, struct C_ShellInfo * aux_shells);

    INTTYPE C_init2(INTTYPE ncenters,
                C_AtomCenter * atoms,
                INTTYPE * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                const char * auxfilename);

    void C_QAO(INTTYPE df_handle, double * matout, INTTYPE matsize);

    void C_cleanup(INTTYPE df_handle);
    void C_cleanup_all(void);

    INTTYPE C_TensorDimensions(INTTYPE df_handle, INTTYPE * d1, INTTYPE * d2, INTTYPE * d3);

    INTTYPE C_CalculateERI(INTTYPE df_handle, double * qso, INTTYPE qsosize, INTTYPE shell1, INTTYPE shell2, INTTYPE shell3, INTTYPE shell4, double * outbuffer, INTTYPE buffersize);

    INTTYPE C_CalculateERIMulti(INTTYPE df_handle,
                            double * qso, INTTYPE qsosize,
                            INTTYPE shell1, INTTYPE nshell1,
                            INTTYPE shell2, INTTYPE nshell2,
                            INTTYPE shell3, INTTYPE nshell3,
                            INTTYPE shell4, INTTYPE nshell4,
                            double * outbuffer, INTTYPE buffersize);

    void C_ReorderQ_GAMESS(INTTYPE df_handle, double * qso, INTTYPE qsosize);
}

    
#endif
