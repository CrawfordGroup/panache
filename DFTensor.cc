/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */
#include <sstream>
#include <fstream>
#include <algorithm>
#include <utility>
#include <cctype>

#include "DFTensor.h"
#include "FittingMetric.h"
#include "Molecule.h"
#include "BasisSet.h"
#include "qt.h"

// for reordering
#include "Orderings.h"
#include "MemorySwapper.h"

#include "ERI.h"

#include "ERDERI.h"
#include "Output.h"

namespace panache
{

DFTensor::DFTensor(std::shared_ptr<BasisSet> primary,
                   std::shared_ptr<BasisSet> auxiliary)
    : primary_(primary), auxiliary_(auxiliary)
{
    common_init();
}

DFTensor::~DFTensor()
{
}

void DFTensor::common_init()
{
    //! \todo remote print?
    print_ = 0;
    debug_ = 0;

    print_header();
    output::printf("\n===========================\nINSIDE BPPSI::DFTENSOR::COMMON_INIT\n===========================\n");

    molecule_ = primary_->molecule();

    build_metric();
}

void DFTensor::print_header()
{
    output::printf("  ==> DF Tensor (by Rob Parrish) <==\n\n");

    output::printf(" => Primary Basis Set <= \n\n");
    primary_->print_detail();

    output::printf(" => Auxiliary Basis Set <= \n\n");
    auxiliary_->print_detail();
}

void DFTensor::build_metric()
{
    std::shared_ptr<FittingMetric> met(new FittingMetric(auxiliary_, true));
    met->form_eig_inverse();
    metric_ = met->get_metric();

    if (debug_)
    {
        //metric_->print();
    }
}

int DFTensor::TensorDimensions(int & d1, int & d2, int & d3)
{
    d1 = d2 = primary_->nbf();
    d3 = auxiliary_->nbf();
    return d1 * d2 * d3;
}


void DFTensor::Qso(double * A, size_t length)
{
    //! \todo I think this should be nbf, but I'm not positive
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();

    if(length != naux * nso * nso)
        throw RuntimeError("Incorrect length array given to Qso");

    SharedMatrix B(new Matrix("Bso", naux, nso * nso));
    double** Bp = B->pointer();
    double** Jp = metric_->pointer();

    output::printf("\n===========================\nINSIDE BPPSI::DFTENSOR::QSO\n===========================\n");


    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    std::shared_ptr<TwoBodyAOInt> eri = GetERI(auxiliary_,zero,primary_,primary_);

    const double* buffer = eri->buffer();

    for (int P = 0; P < auxiliary_->nshell(); P++)
    {
        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        for (int M = 0; M < primary_->nshell(); M++)
        {
            int nm = primary_->shell(M).nfunction();
            int mstart = primary_->shell(M).function_index();
            for (int N = 0; N < primary_->nshell(); N++)
            {
                int nn = primary_->shell(N).nfunction();
                int nstart = primary_->shell(N).function_index();

                eri->compute_shell(P,0,M,N);

                for (int p = 0, index = 0; p < np; p++)
                {
                    for (int m = 0; m < nm; m++)
                    {
                        for (int n = 0; n < nn; n++, index++)
                        {
                            Bp[p + pstart][(m + mstart) * nso + (n + nstart)] = buffer[index];
                        }
                    }
                }
            }
        }
    }

//    C_DGEMM('N','N',naux, nso * nso, naux, 1.0, Jp[0], naux, Bp[0], nso * nso, 0.0,
//            A, nso * nso);

    // Do the transpose of above - this makes ERI generation "faster"
    C_DGEMM('T','T',nso*nso, naux, naux, 1.0, Bp[0], nso*nso, Jp[0], naux, 0.0,
            A, naux);

    if (debug_)
    {
        //metric_->print();
        //B->print();
        //A->print();
    }
}

int DFTensor::CalculateERI(double * qso, int qsosize, int shell1, int shell2, int shell3, int shell4, double * outbuffer, int buffersize)
{
    //! \todo do something with qsosize

    int nfa = primary_->shell(shell1).nfunction();
    int astart = primary_->shell(shell1).function_index();

    int nfb = primary_->shell(shell2).nfunction();
    int bstart = primary_->shell(shell2).function_index();

    int nfc = primary_->shell(shell3).nfunction();
    int cstart = primary_->shell(shell3).function_index();

    int nfd = primary_->shell(shell4).nfunction();
    int dstart = primary_->shell(shell4).function_index();

    int nint = nfa * nfb * nfc * nfd;

    int nbf = primary_->nbf();
    int naux = auxiliary_->nbf();


    if(nint > buffersize)
        throw RuntimeError("Error - ERI buffer not large enough!");

    int bufindex = 0;

    //!\todo replace with DGEMM?
    for(int a = 0; a < nfa; a++)
        for(int b = 0; b < nfb; b++)
            for(int c = 0; c < nfc; c++)
                for(int d = 0; d < nfd; d++)
                {
                    outbuffer[bufindex] = C_DDOT(naux,
                                                 qso + (astart+a)*nbf*naux+(bstart+b)*naux, 1,
                                                 qso + (cstart+c)*nbf*naux+(dstart+d)*naux, 1);
                    bufindex++;
                }

    return nint;
}


int DFTensor::CalculateERIMulti(double * qso, int qsosize,
                                int shell1, int nshell1,
                                int shell2, int nshell2,
                                int shell3, int nshell3,
                                int shell4, int nshell4,
                                double * outbuffer, int buffersize)
{
    //! \todo do something with qsosize
    int nint = 0;

    int nbf = primary_->nbf();
    int naux = auxiliary_->nbf();

    int bufindex = 0;

    for(int i = 0; i < nshell1; i++)
    {
        int nfa = primary_->shell(shell1+i).nfunction();
        int astart = primary_->shell(shell1+i).function_index();

        for(int a = 0; a < nfa; a++)
        {
            for(int j = 0; j < nshell2; j++)
            {
                int nfb = primary_->shell(shell2+j).nfunction();
                int bstart = primary_->shell(shell2+j).function_index();

                for(int b = 0; b < nfb; b++)
                {
                    for(int k = 0; k < nshell3; k++)
                    {
                        int nfc = primary_->shell(shell3+k).nfunction();
                        int cstart = primary_->shell(shell3+k).function_index();

                        for(int c = 0; c < nfc; c++)
                        {
                            for(int l = 0; l < nshell4; l++)
                            {
                                int nfd = primary_->shell(shell4+l).nfunction();
                                int dstart = primary_->shell(shell4+l).function_index();

                                nint += nfd;

                                if(nint > buffersize)
                                    throw RuntimeError("Error - ERI buffer not large enough!");

                                for(int d = 0; d < nfd; d++)
                                {
                                    outbuffer[bufindex] = C_DDOT(naux,
                                                                 qso + (astart+a)*nbf*naux+(bstart+b)*naux, 1,
                                                                 qso + (cstart+c)*nbf*naux+(dstart+d)*naux, 1);
                                    bufindex++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return nint;
}




// note - passing by value for the vector
static void Reorder(std::vector<unsigned short> order, std::vector<double *> pointers,
                    reorder::MemorySwapper & sf)
{
    size_t size = order.size();

    // original order is 1 2 3 4 5 6....
    std::vector<unsigned short> currentorder(size);

    for(int i = 0; i < size; i++)
        currentorder[i] = i+1;

    for(int i = 0; i < size; i++)
    {
        // find the index in the current order
        size_t cindex = 0;
        bool found = false;

        for(int j = 0; j < size; j++)
        {
            if(currentorder[j] == order[i])
            {
                found = true;
                cindex = j;
                break;
            }
        }
        if(!found)
            throw RuntimeError("Error in reordering - index not found?");


        // we shouldn't swap anything that was previously put in place...
        if(cindex < i)
            throw RuntimeError("Error in reordering - going to swap something I shouldn't");

        //swap
        if(cindex != i)
        {
            sf.swap(pointers[i], pointers[cindex]);
            std::swap(currentorder[i], currentorder[cindex]);
        }
    }

    // double check
    for(int i = 0; i < size; i++)
    {
        if(currentorder[i] != order[i])
            throw RuntimeError("Reordering failed!");
    }
}




void DFTensor::ReorderQ(double * qso, int qsosize, const reorder::Orderings & order)
{
    using namespace reorder;
    using std::placeholders::_1;
    using std::placeholders::_2;

    int nso = primary_->nbf();
    int nq = auxiliary_->nbf();

    if((nq * nso * nso) != qsosize)
        throw RuntimeError("Incompatible Qso matrix in ReorderQ");

    // Dimensions on Q:
    // nso * nso * nq


    TotalMemorySwapper sf1(nq);  // swaps rows
    LimitedMemorySwapper sf2(nso*nq, 1e7); // swaps 'tables' (limited to ~80MB extra memory)

    std::vector<PointerMap> vpm;

    //go through what would need to be changed in the primary basis
    for(int i = 0; i < primary_->nshell(); i++)
    {
        const GaussianShell & s = primary_->shell(i);
        if(order.NeedsReordering(s.am()))
            vpm.push_back(PointerMap(s.function_index(), order.GetOrder(s.am())));
    }


    std::vector<double *> pointers(primary_->max_function_per_shell());


    // for each i
    for(size_t i = 0; i < nso; i++)
    {
        // Swap rows
        for(auto & it : vpm)
        {
            size_t ntoswap = it.order.size();
            size_t start = i*nso*nq;

            for(size_t n = 0; n < ntoswap; n++)
                pointers[n] = qso + start + (it.start+n)*nq;

            Reorder(it.order, pointers, sf1);

        }

    }


    // swap 'tables'
    for(auto & it : vpm)
    {
        size_t ntoswap = it.order.size();

        for(size_t n = 0; n < ntoswap; n++)
            pointers[n] = qso + (it.start+n)*nso*nq;

        Reorder(it.order, pointers, sf2);

    }

}

void DFTensor::ReorderQ_GAMESS(double * qso, int qsosize)
{
    reorder::GAMESS_Ordering go;
    ReorderQ(qso, qsosize, go);
}


}







