/*! \file
 *  \brief C interface to the PANACHE library (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <map>

#include "c_interface.h"
#include "c_convert.h"
#include "DFTensor.h"
#include "Output.h"
#include "Exception.h"
#include "BasisSetParser.h"

using panache::RuntimeError;
using panache::BasisSet;
using panache::Gaussian94BasisSetParser;
using panache::SharedBasisSet;
using panache::DFTensor;
using panache::Molecule;
using panache::SharedMolecule;

namespace
{

int tensor_index_ = 0;
std::map<int, DFTensor *> dftensors_;

} // close anonymous namespace


extern "C" {


    int_t panache_init(int_t ncenters,
               C_AtomCenter * atoms, int_t normalized,
               int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               int_t * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
               const char * filename)
    {
        // Molecule
        SharedMolecule molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        auto primaryBasis = BasisSetFromArrays(molecule, ncenters,
                            primary_nshellspercenter, primary_shells, normalized);

        auto auxBasis = BasisSetFromArrays(molecule, ncenters,
                        aux_nshellspercenter, aux_shells, normalized);

        DFTensor * dft = new DFTensor(primaryBasis, auxBasis, filename);
        dftensors_[tensor_index_] = dft;

        return tensor_index_++;
    }



    int_t panache_init2(int_t ncenters,
                    C_AtomCenter * atoms, int_t normalized,
                    int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                    const char * auxfilename, const char * matfilename)
    {
        // Molecule
        SharedMolecule molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        SharedBasisSet primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                      primary_nshellspercenter, primary_shells, normalized);

        // Gaussian input file parser for the auxiliary basis
        std::shared_ptr<Gaussian94BasisSetParser> parser(new Gaussian94BasisSetParser);
        SharedBasisSet auxBasis(new BasisSet(parser, molecule, auxfilename));

        DFTensor * dft = new DFTensor(primaryBasis, auxBasis, matfilename);
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
            throw RuntimeError("Error - cannot erase DFTensor object with that handle!");

    }



    void panache_cleanup_all(void)
    {
        for(auto it : dftensors_)
            delete it.second;

        dftensors_.clear();
    }

    void panache_setcmatrix(int_t df_handle, double * cmo, int_t nmo, int_t cmo_is_trans)
    {
        // matsize is checked in here
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->SetCMatrix(cmo, nmo, cmo_is_trans);
    }

    void panache_setoutputbuffer(int_t df_handle, double * buffer, int_t bufsize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->SetOutputBuffer(buffer, bufsize);
    
    }

    void panache_genqso(int_t df_handle, int_t inmem)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->GenQso(inmem);
    }


    int_t panache_getbatch_qso(int_t df_handle)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        return dftensors_[df_handle]->GetBatch_Qso();
    }

    int_t panache_getbatch_qmo(int_t df_handle)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        return dftensors_[df_handle]->GetBatch_Qmo();
    }


    int_t panache_qsodimensions(int_t df_handle, int_t * naux, int_t * nso2)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        //DFTensor class takes d1-d3 by reference
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
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");
        return dftensors_[df_handle]->SetNThread(nthread); 
    }


/*
    int_t panache_CalculateERI(int_t df_handle, double * qso, int_t qsosize, int_t shell1, int_t shell2, int_t shell3, int_t shell4, double * outbuffer, int_t buffersize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->CalculateERI(qso, qsosize, shell1, shell2, shell3, shell4, outbuffer, buffersize);
    }


    int_t panache_CalculateERIMulti(int_t df_handle,
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

    void panache_ReorderQ_GAMESS(int_t df_handle, double * qso, int_t qsosize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->ReorderQ_GAMESS(qso, qsosize);
    }
*/
}

