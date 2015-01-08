/*! \file
 *  \brief C interface to the PANACHE library (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <map>
#include <sstream>
#include <iostream> // for std::cout

// included from c_convert.h
//#include "panache/c_interface.h"
#include "panache/c_convert.h"
#include "panache/DFTensor.h"
#include "panache/CHTensor.h"
#include "panache/Output.h"
#include "panache/Exception.h"
#include "panache/BasisSetParser.h"

using panache::RuntimeError;
using panache::BasisSet;
using panache::SharedBasisSetParser;
using panache::Gaussian94BasisSetParser;
using panache::SharedBasisSet;
using panache::ThreeIndexTensor;
using panache::DFTensor;
using panache::CHTensor;
using panache::Molecule;
using panache::SharedMolecule;

namespace
{

int tensor_index_ = 0;
std::map<int, ThreeIndexTensor *> xtensors_;


void CheckHandle(int_t df_handle, const char * func)
{
    if(xtensors_.count(df_handle) == 0)
    {
        std::stringstream ss;
        ss << "Function: " << func << ": Error - cannot find ThreeIndexTensor object with that handle!";
        throw RuntimeError(ss.str());
    }
}


} // close anonymous namespace


extern "C" {


    int_t panache_dfinit(int_t ncenters,
               C_AtomCenter * atoms,
               int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               int_t * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
               const char * directory, int_t metricflag, int_t bsorder, int_t nthreads )
    {
        // Molecule
        SharedMolecule molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        auto primaryBasis = BasisSetFromArrays(molecule, ncenters,
                            primary_nshellspercenter, primary_shells);

        auto auxBasis = BasisSetFromArrays(molecule, ncenters,
                        aux_nshellspercenter, aux_shells);

        ThreeIndexTensor * dft = new DFTensor(primaryBasis, auxBasis, directory, metricflag, bsorder, nthreads);
        xtensors_[tensor_index_] = dft;

        return tensor_index_++;
    }



    int_t panache_dfinit2(int_t ncenters,
                    C_AtomCenter * atoms,
                    int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                    const char * auxfilename, const char * directory, int_t metricflag, int_t bsorder, int_t nthreads)
    {
        // Molecule
        SharedMolecule molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        SharedBasisSet primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                      primary_nshellspercenter, primary_shells);

        // Gaussian input file parser for the auxiliary basis
        SharedBasisSetParser parser(new Gaussian94BasisSetParser);
        SharedBasisSet auxBasis(new BasisSet(parser, molecule, auxfilename));

        ThreeIndexTensor * dft = new DFTensor(primaryBasis, auxBasis, directory, metricflag, bsorder, nthreads);
        xtensors_[tensor_index_] = dft;

        return tensor_index_++;
    }


    int_t panache_chinit(int_t ncenters,
                    C_AtomCenter * atoms,
                    int_t * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                    double delta, const char * directory, int_t bsorder, int_t nthreads)
    {
        // Molecule
        SharedMolecule molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        SharedBasisSet primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                                      primary_nshellspercenter, primary_shells);

        ThreeIndexTensor * dft = new CHTensor(primaryBasis, delta, directory, bsorder, nthreads);
        xtensors_[tensor_index_] = dft;

        return tensor_index_++;
    }


    int_t panache_qbatchsize(int_t df_handle, int_t tensorflag)
    {
        CheckHandle(df_handle, __FUNCTION__);
        return xtensors_[df_handle]->QBatchSize(tensorflag); 
    }

    int_t panache_batchsize(int_t df_handle, int_t tensorflag)
    {
        CheckHandle(df_handle, __FUNCTION__);
        return xtensors_[df_handle]->BatchSize(tensorflag); 
    }

    int_t panache_ispacked(int_t df_handle, int_t tensorflag)
    {
        CheckHandle(df_handle, __FUNCTION__);
        return xtensors_[df_handle]->IsPacked(tensorflag); 
    }

    int_t panache_calcindex(int_t df_handle, int_t tensorflag, int_t i, int_t j)
    {
        CheckHandle(df_handle, __FUNCTION__);
        return xtensors_[df_handle]->CalcIndex(tensorflag, i, j); 
    }

    int_t panache_tensordimensions(int_t df_handle, int_t tensorflag,
                                   int_t * naux, int_t * ndim1, int_t * ndim2)
    {
        CheckHandle(df_handle, __FUNCTION__);
        int nauxtmp, ndim1tmp, ndim2tmp;
        int ret = xtensors_[df_handle]->TensorDimensions(tensorflag, nauxtmp, ndim1tmp, ndim2tmp);
        *naux = nauxtmp;
        *ndim1 = ndim1tmp;
        *ndim2 = ndim2tmp;
        return ret;
    }

    void panache_cleanup(int_t df_handle)
    {
        if(xtensors_.count(df_handle) > 0)
        {
            delete xtensors_[df_handle];
            xtensors_.erase(df_handle);
        }
    }



    void panache_cleanup_all(void)
    {
        for(auto it : xtensors_)
            delete it.second;

        xtensors_.clear();
    }

    void panache_setcmatrix(int_t df_handle, double * cmo, int_t nmo, int_t cmo_is_trans)
    {
        CheckHandle(df_handle, __FUNCTION__);
        xtensors_[df_handle]->SetCMatrix(cmo, nmo, cmo_is_trans);
    }

    void panache_genqtensors(int_t df_handle, int_t qflags, int_t storeflags)
    {
        CheckHandle(df_handle, __FUNCTION__);
        xtensors_[df_handle]->GenQTensors(qflags, storeflags); 
    }

    void panache_delete(int_t df_handle, int_t qflags)
    {
        CheckHandle(df_handle, __FUNCTION__);
        xtensors_[df_handle]->Delete(qflags); 
    }

    void panache_output(FILE * out)
    {
       panache::output::SetOutput(out); 
    }

    void panache_stdout(void)
    {
       panache::output::SetOutput(&std::cout); 
    }

    int_t panache_setnthread(int_t df_handle, int_t nthread)
    {
        CheckHandle(df_handle, __FUNCTION__);
        return xtensors_[df_handle]->SetNThread(nthread); 
    }


    void panache_setnocc(int_t df_handle, int_t nocc, int_t nfroz)
    {
        CheckHandle(df_handle, __FUNCTION__);
        xtensors_[df_handle]->SetNOcc(nocc, nfroz); 
    }

    void panache_printtimings(int_t df_handle)
    {
        CheckHandle(df_handle, __FUNCTION__);
        xtensors_[df_handle]->PrintTimings(); 
    }

    int_t panache_getqbatch(int_t df_handle, int_t tensorflag, double * outbuf, int_t bufsize, int_t qstart)
    {
        CheckHandle(df_handle, __FUNCTION__);
        return xtensors_[df_handle]->GetQBatch(tensorflag, outbuf, bufsize, qstart);
    }


    int_t panache_getbatch(int_t df_handle, int_t tensorflag, double * outbuf, int_t bufsize, int_t ijstart)
    {
        CheckHandle(df_handle, __FUNCTION__);
        return xtensors_[df_handle]->GetBatch(tensorflag, outbuf, bufsize, ijstart);
    }
}

