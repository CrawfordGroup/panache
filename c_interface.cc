#include <map>

#include "c_interface.h"
#include "c_convert.h"
#include "DFTensor.h"
#include "Exception.h"

using panache::RuntimeError;

namespace {

int tensor_index_ = 0;
std::map<int, panache::DFTensor *> dftensors_;

} // close anonymous namespace


extern "C" {


    int C_init(int ncenters,
               C_AtomCenter * atoms,
               int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               int * aux_nshellspercenter, struct C_ShellInfo * aux_shells)
    {
        // Molecule
        std::shared_ptr<panache::Molecule> molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        auto primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                                        primary_nshellspercenter, primary_shells);

        auto auxBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                                    aux_nshellspercenter, aux_shells);

        panache::DFTensor * dft = new panache::DFTensor(primaryBasis, auxBasis);
        dftensors_[tensor_index_] = dft;

        return tensor_index_++;
    }

    void C_cleanup(int df_handle)
    {
        if(dftensors_.count(df_handle) > 0)
        {
            delete dftensors_[df_handle];
            dftensors_.erase(df_handle);
        }
        else
            throw RuntimeError("Error - cannot erase DFTensor object with that handle!");
            
    }

    void C_cleanup_all(void)
    {
        for(auto it : dftensors_)
            delete it.second;

        dftensors_.clear();
    }

    void C_QAO(int df_handle, double * matout, int matsize)
    {
        // matsize is checked in here
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->Qso(matout, matsize);
    }

}

