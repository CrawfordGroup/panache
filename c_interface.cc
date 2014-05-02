#include "c_interface.h"
#include "c_convert.h"
#include "DFTensor.h"



extern "C" {



    double * C_QAO(int ncenters,
                 C_AtomCenter * atoms,
                 int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                 int * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
                 int * nrow_out, int * ncol_out)
    {
        // Molecule
        std::shared_ptr<panache::Molecule> molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        auto primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                                        primary_nshellspercenter, primary_shells);

        auto auxBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                                    aux_nshellspercenter, aux_shells);

        panache::DFTensor dft(primaryBasis, auxBasis);
        auto mat = dft.Qso();

        // convert the resulting matrix
        *nrow_out = mat->rowspi()[0];
        *ncol_out = mat->colspi()[0]; 
        //return mat->give_up();

        double * ret = (double *)malloc(mat->rowspi()[0]*mat->colspi()[0] * sizeof(double));
        memcpy(ret, &(mat->pointer(0)[0][0]), mat->rowspi()[0]*mat->colspi()[0]*sizeof(double));
        return ret;

    }




    void free_matrix(double * mat)
    {
        free(mat);
    }

}

