#include "c_interface.h"
#include "c_convert.h"
#include "DFTensor.h"



extern "C" {



    void C_QAO(int ncenters,
               C_AtomCenter * atoms,
               int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               int * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
               double * matout, int matsize)
    {
        // Molecule
        std::shared_ptr<panache::Molecule> molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        auto primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                                        primary_nshellspercenter, primary_shells);

        auto auxBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                                    aux_nshellspercenter, aux_shells);

        panache::DFTensor dft(primaryBasis, auxBasis);

        // matsize is checked in here
        dft.Qso(matout, matsize);
    }

}

