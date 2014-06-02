#include <map>

#include "c_interface.h"
#include "c_convert.h"
#include "DFTensor.h"

#include "Exception.h"
#include "BasisSetParser.h"

using panache::RuntimeError;

namespace
{

int tensor_index_ = 0;
std::map<int, panache::DFTensor *> dftensors_;

} // close anonymous namespace


extern "C" {


    int_t C_init(int_t ncenters,
               C_AtomCenter * atoms,
               int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               int_t * aux_nshellspercenter, struct C_ShellInfo * aux_shells)
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


    int_t C_init2(int_t ncenters,
                    C_AtomCenter * atoms,
                    int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                    const char * auxfilename)
    {
        // Molecule
        std::shared_ptr<panache::Molecule> molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        auto primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                            primary_nshellspercenter, primary_shells);

        // Gaussian input file parser for the auxiliary basis
        std::shared_ptr<panache::Gaussian94BasisSetParser> parser(new panache::Gaussian94BasisSetParser);
        auto auxBasis = panache::BasisSet::construct(parser, molecule, auxfilename);

        panache::DFTensor * dft = new panache::DFTensor(primaryBasis, auxBasis);
        dftensors_[tensor_index_] = dft;

        return tensor_index_++;
    }




    void C_cleanup(int_t df_handle)
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

    void C_QAO(int_t df_handle, double * matout, int_t matsize)
    {
        // matsize is checked in here
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->Qso(matout, matsize);
    }

    int_t C_TensorDimensions(int_t df_handle, int_t * d1, int_t * d2, int_t * d3)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        //DFTensor class takes d1-d3 by reference
        // use regular ints first
        int t1, t2, t3;
        dftensors_[df_handle]->TensorDimensions(t1, t2, t3);

        *d1 = t1;
        *d2 = t2;
        *d3 = t3;
        return t1 * t2 * t3;
    }


    int_t C_CalculateERI(int_t df_handle, double * qso, int_t qsosize, int_t shell1, int_t shell2, int_t shell3, int_t shell4, double * outbuffer, int_t buffersize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->CalculateERI(qso, qsosize, shell1, shell2, shell3, shell4, outbuffer, buffersize);
    }


    int_t C_CalculateERIMulti(int_t df_handle,
                            double * qso, int_t qsosize,
                            int_t shell1, int_t nshell1,
                            int_t shell2, int_t nshell2,
                            int_t shell3, int_t nshell3,
                            int_t shell4, int_t nshell4,
                            double * outbuffer, int_t buffersize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        return dftensors_[df_handle]->CalculateERIMulti(qso, qsosize, shell1, nshell1, shell2, nshell2,
                                                        shell3, nshell3, shell4, nshell4, outbuffer, buffersize);
    }

    void C_ReorderQ_GAMESS(int_t df_handle, double * qso, int_t qsosize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->ReorderQ_GAMESS(qso, qsosize);
    }
}

