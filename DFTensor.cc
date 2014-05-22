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

#ifdef USE_LIBINT
#include "ERI.h"
#endif

#ifdef USE_LIBINT2
#include "ERI2.h"
#endif

#ifdef USE_SLOWERI
#include "SlowERI.h"
#endif

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
    //primary_->print_by_level(outfile,print_);
    //fflush(outfile);

    output::printf(" => Auxiliary Basis Set <= \n\n");
    auxiliary_->print_detail();
    //auxiliary_->print_by_level(outfile,print_);
    //fflush(outfile);
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
    d1 = auxiliary_->nbf();
    d2 = d3 = primary_->nbf();
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

#ifdef USE_LIBINT
    ERI eri(auxiliary_,zero,primary_,primary_);
#endif

#ifdef USE_LIBINT2
    ERI2 eri(auxiliary_,zero,primary_,primary_);
#endif

#ifdef USE_SLOWERI
    SlowERI eri(auxiliary_,zero,primary_,primary_);
#endif

    //ERDERI eri(auxiliary_,zero,primary_,primary_);

    const double* buffer = eri.buffer();

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

                eri.compute_shell(P,0,M,N);

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

    C_DGEMM('N','N',naux, nso * nso, naux, 1.0, Jp[0], naux, Bp[0], nso * nso, 0.0,
            A, nso * nso);

    if (debug_)
    {
        //metric_->print();
        //B->print();
        //A->print();
    }
}

int DFTensor::ERI(double * qso, int qsosize, int shell1, int shell2, int shell3, int shell4, double * outbuffer, int buffersize)
{
    //! \todo do something with qsosize

    int nfp = primary_->shell(shell1).nfunction();
    int pstart = primary_->shell(shell1).function_index();

    int nfq = primary_->shell(shell2).nfunction();
    int qstart = primary_->shell(shell2).function_index();

    int nfr = primary_->shell(shell3).nfunction();
    int rstart = primary_->shell(shell3).function_index();

    int nfs = primary_->shell(shell4).nfunction();
    int sstart = primary_->shell(shell4).function_index();

    int nint = nfp * nfq * nfr * nfs;

    int nbf = primary_->nbf();
    int nbf2 = nbf*nbf;


    if(nint > buffersize)
        throw RuntimeError("Error - ERI buffer not large enough!");

    for(int i = 0; i < nint; i++)
        outbuffer[i] = 0.0;


    int naux = auxiliary_->nbf();

    int bufindex = 0;

    //!\todo replace with BLAS call
    for(int a = 0; a < nfp; a++)
    for(int b = 0; b < nfq; b++)
    for(int c = 0; c < nfr; c++)
    for(int d = 0; d < nfs; d++)
    {
        for(int Q = 0; Q < naux; Q++)
        {
            outbuffer[bufindex] += qso[Q*nbf2+(pstart+a)*nbf+(qstart+b)]
                                 * qso[Q*nbf2+(rstart+c)*nbf+(sstart+d)];
        }
        bufindex++;
    }

    return nint;
}

}

