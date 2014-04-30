extern "C" {

    struct C_ShellInfo
    {
        int nprim;
        int am;
        int ispure;
        double * exp; // of length nprim
        double * coef; // of length nprim
    };


    double * C_QAO(int ncenters, int nao, int nmo,
                 double * Z, char const ** symbols,
                 double * masses, double ** centers,
                 int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                 int * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
                 int * nrow_out, int * ncol_out);

    void free_matrix(double * mat);

}

