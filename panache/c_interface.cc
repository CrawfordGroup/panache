/*! \file
 *  \brief C interface to the PANACHE library (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <map>

#include "panache/c_interface.h"
#include "panache/c_convert.h"
#include "panache/DFTensor2.h"
#include "panache/Output.h"
#include "panache/Exception.h"
#include "panache/BasisSetParser.h"

using panache::RuntimeError;
using panache::BasisSet;
using panache::Gaussian94BasisSetParser;
using panache::SharedBasisSet;
using panache::DFTensor2;
using panache::Molecule;
using panache::SharedMolecule;

namespace
{

int tensor_index_ = 0;
std::map<int, DFTensor2 *> dftensors_;

} // close anonymous namespace


extern "C" {


    int_t panache_init(int_t ncenters,
               C_AtomCenter * atoms, int_t normalized,
               int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               int_t * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
               const char * filename, int_t nthreads )
    {
        // Molecule
        SharedMolecule molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        auto primaryBasis = BasisSetFromArrays(molecule, ncenters,
                            primary_nshellspercenter, primary_shells, normalized);

        auto auxBasis = BasisSetFromArrays(molecule, ncenters,
                        aux_nshellspercenter, aux_shells, normalized);

        DFTensor2 * dft = new DFTensor2(primaryBasis, auxBasis, filename, nthreads);
        dftensors_[tensor_index_] = dft;

        return tensor_index_++;
    }



    int_t panache_init2(int_t ncenters,
                    C_AtomCenter * atoms, int_t normalized,
                    int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                    const char * auxfilename, const char * matfilename, int_t nthreads)
    {
        // Molecule
        SharedMolecule molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        SharedBasisSet primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                      primary_nshellspercenter, primary_shells, normalized);

        // Gaussian input file parser for the auxiliary basis
        std::shared_ptr<Gaussian94BasisSetParser> parser(new Gaussian94BasisSetParser);
        SharedBasisSet auxBasis(new BasisSet(parser, molecule, auxfilename));

        DFTensor2 * dft = new DFTensor2(primaryBasis, auxBasis, matfilename, nthreads);
        dftensors_[tensor_index_] = dft;

        return tensor_index_++;
    }




    void panache_cleanup(int_t df_handle)
    {
        if(dftensors_.count(df_handle) > 0)
        {
            delete dftensors_[df_handle];
            dftensors_.erase(df_handle);
        }
        else
            throw RuntimeError("Error - cannot erase DFTensor2 object with that handle!");

    }



    void panache_cleanup_all(void)
    {
        for(auto it : dftensors_)
            delete it.second;

        dftensors_.clear();
    }

    void panache_setcmatrix(int_t df_handle, double * cmo, int_t nmo, int_t cmo_is_trans, int_t bsorder)
    {
        // matsize is checked in here
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");

        dftensors_[df_handle]->SetCMatrix(cmo, nmo, cmo_is_trans, bsorder);
    }

    void panache_genqtensors(int_t df_handle, int_t qflags, int_t storetype)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");

        return dftensors_[df_handle]->GenQTensors(qflags, storetype); 
    }


    int_t panache_getqbatch_qso(int_t df_handle, double * outbuf, int bufsize, int qstart)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");

        return dftensors_[df_handle]->GetQBatch_Qso(outbuf, bufsize, qstart);
    }

    int_t panache_getqbatch_qmo(int_t df_handle, double * outbuf, int bufsize, int qstart)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");

        return dftensors_[df_handle]->GetQBatch_Qmo(outbuf, bufsize, qstart);
    }

    int_t panache_getqbatch_qoo(int_t df_handle, double * outbuf, int bufsize, int qstart)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");

        return dftensors_[df_handle]->GetQBatch_Qoo(outbuf, bufsize, qstart);
    }

    int_t panache_getqbatch_qov(int_t df_handle, double * outbuf, int bufsize, int qstart)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");

        return dftensors_[df_handle]->GetQBatch_Qov(outbuf, bufsize, qstart);
    }

    int_t panache_getqbatch_qvv(int_t df_handle, double * outbuf, int bufsize, int qstart)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");

        return dftensors_[df_handle]->GetQBatch_Qvv(outbuf, bufsize, qstart);
    }


    int_t panache_qsodimensions(int_t df_handle, int_t * naux, int_t * nso2)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");

        //DFTensor2 class takes d1-d3 by reference
        // use regular ints first
        int t1, t2;
        int tot = dftensors_[df_handle]->QsoDimensions(t1, t2);

        *naux = t1;
        *nso2 = t2;
        return tot;
    }

    void panache_output(FILE * out)
    {
       panache::output::SetOutput(out); 
    }

    int_t panache_setnthread(int_t df_handle, int_t nthread)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");
        return dftensors_[df_handle]->SetNThread(nthread); 
    }

    void panache_setnocc(int_t df_handle, int_t nocc, int_t nfroz)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor2 object with that handle!");
        dftensors_[df_handle]->SetNOcc(nocc, nfroz); 
    }

}

