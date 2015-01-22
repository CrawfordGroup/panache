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


void CheckHandle(int handle, const char * func)
{
    if(xtensors_.count(handle) == 0)
    {
        std::stringstream ss;
        ss << "Function: " << func << ": Error - cannot find ThreeIndexTensor object with that handle!";
        throw RuntimeError(ss.str());
    }
}


} // close anonymous namespace


extern "C" {


    int panache_dfinit(int ncenters,
               C_AtomCenter * atoms,
               int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               int * aux_nshellspercenter, struct C_ShellInfo * aux_shells,
               const char * directory, int metricflag, int bsorder, int nthreads )
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



    int panache_dfinit2(int ncenters,
                    C_AtomCenter * atoms,
                    int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                    const char * auxfilename, const char * directory, int metricflag, int bsorder, int nthreads)
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


    int panache_chinit(int ncenters,
                    C_AtomCenter * atoms,
                    int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                    double delta, const char * directory, int bsorder, int nthreads)
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


    int panache_qbatchsize(int handle, int tensorflag)
    {
        CheckHandle(handle, __FUNCTION__);
        return xtensors_[handle]->QBatchSize(tensorflag); 
    }

    int panache_batchsize(int handle, int tensorflag)
    {
        CheckHandle(handle, __FUNCTION__);
        return xtensors_[handle]->BatchSize(tensorflag); 
    }

    int panache_ispacked(int handle, int tensorflag)
    {
        CheckHandle(handle, __FUNCTION__);
        return xtensors_[handle]->IsPacked(tensorflag); 
    }

    int panache_calcindex(int handle, int tensorflag, int i, int j)
    {
        CheckHandle(handle, __FUNCTION__);
        return xtensors_[handle]->CalcIndex(tensorflag, i, j); 
    }

    int panache_tensordimensions(int handle, int tensorflag,
                                   int * naux, int * ndim1, int * ndim2)
    {
        CheckHandle(handle, __FUNCTION__);
        int ret = xtensors_[handle]->TensorDimensions(tensorflag, *naux, *ndim1, *ndim2);
        return ret;
    }

    void panache_cleanup(int handle)
    {
        if(xtensors_.count(handle) > 0)
        {
            delete xtensors_[handle];
            xtensors_.erase(handle);
        }
    }



    void panache_cleanup_all(void)
    {
        for(auto it : xtensors_)
            delete it.second;

        xtensors_.clear();
    }

    void panache_setcmatrix(int handle, double * cmo, int nmo, int cmo_is_trans)
    {
        CheckHandle(handle, __FUNCTION__);
        xtensors_[handle]->SetCMatrix(cmo, nmo, cmo_is_trans);
    }

    void panache_genqtensors(int handle, int qflags, int storeflags)
    {
        CheckHandle(handle, __FUNCTION__);
        xtensors_[handle]->GenQTensors(qflags, storeflags); 
    }

    void panache_delete(int handle, int qflags)
    {
        CheckHandle(handle, __FUNCTION__);
        xtensors_[handle]->Delete(qflags); 
    }

    void panache_output(FILE * out)
    {
       panache::output::SetOutput(out); 
    }

    void panache_stdout(void)
    {
       panache::output::SetOutput(&std::cout); 
    }

    int panache_setnthread(int handle, int nthread)
    {
        CheckHandle(handle, __FUNCTION__);
        return xtensors_[handle]->SetNThread(nthread); 
    }


    void panache_setnocc(int handle, int nocc, int nfroz)
    {
        CheckHandle(handle, __FUNCTION__);
        xtensors_[handle]->SetNOcc(nocc, nfroz); 
    }

    void panache_printtimings(int handle)
    {
        CheckHandle(handle, __FUNCTION__);
        xtensors_[handle]->PrintTimings(); 
    }

    int panache_getqbatch(int handle, int tensorflag, double * outbuf, int bufsize, int qstart)
    {
        CheckHandle(handle, __FUNCTION__);
        return xtensors_[handle]->GetQBatch(tensorflag, outbuf, bufsize, qstart);
    }


    int panache_getbatch(int handle, int tensorflag, double * outbuf, int bufsize, int ijstart)
    {
        CheckHandle(handle, __FUNCTION__);
        return xtensors_[handle]->GetBatch(tensorflag, outbuf, bufsize, ijstart);
    }
}

