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

#include <algorithm>
#include <utility>
#include <cstring> // memset

#include "ERI.h"

#include "FittingMetric.h"
#include "BasisSet.h"
#include "Lapack.h"

#ifdef _OPENMP
#include <omp.h>
#endif


namespace panache
{

FittingMetric::FittingMetric(std::shared_ptr<BasisSet> aux) :
    aux_(aux), naux_(aux->nbf()),is_poisson_(false), is_inverted_(false), omega_(0.0),
    metric_(new double[naux_*naux_])
{
}

FittingMetric::FittingMetric(std::shared_ptr<BasisSet> aux, double omega) :
    aux_(aux), naux_(aux->nbf()),is_poisson_(false), is_inverted_(false), omega_(omega)
{
}

FittingMetric::~FittingMetric()
{
    delete [] metric_;
}

void FittingMetric::form_fitting_metric()
{
    is_inverted_ = false;
    algorithm_ = "NONE";

    // Build the full DF/Poisson matrix in the AO basis first
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    // Only thread if not already in parallel (handy for local fitting)
    int nthread = 1;

    #ifdef _OPENMP
    nthread = omp_get_max_threads();
    #endif

    // == (A|B) Block == //
    const double **Jbuffer = new const double*[nthread];
    std::shared_ptr<TwoBodyAOInt> *Jint = new std::shared_ptr<TwoBodyAOInt>[nthread];

    for (int Q = 0; Q<nthread; Q++)
    {
        if (omega_ > 0.0)
        {
            Jint[Q] = GetErfERI(omega_, aux_, zero, aux_, zero);
        }
        else
        {
            Jint[Q] = GetERI(aux_, zero, aux_, zero);
        }
        Jbuffer[Q] = Jint[Q]->buffer();
    }

    for (int MU=0; MU < aux_->nshell(); ++MU)
    {
        int nummu = aux_->shell(MU).nfunction();

        int thread = 0;
        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif

        for (int NU=0; NU <= MU; ++NU)
        {
            int numnu = aux_->shell(NU).nfunction();

            Jint[thread]->compute_shell(MU, 0, NU, 0);

            int index = 0;
            for (int mu=0; mu < nummu; ++mu)
            {
                int omu = aux_->shell(MU).function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index)
                {
                    int onu = aux_->shell(NU).function_index() + nu;

                    //! \todo packed storage?
                    metric_[omu*naux_+onu] = metric_[onu*naux_+omu] = Jbuffer[thread][index];
                }
            }
        }
    }
    delete[] Jbuffer;
    delete[] Jint;

    pivots_.resize(naux_);
    rev_pivots_.resize(naux_);
    for (int Q = 0; Q < naux_; Q++)
        pivots_[Q] = rev_pivots_[Q] = Q;
}

void FittingMetric::form_eig_inverse(double tol)
{
    is_inverted_ = true;
    algorithm_ = "EIG";

    form_fitting_metric();

    //metric_->print();

    int n = naux_;
    int n2 = n*n;

    // Copy J to W
    double * W = new double[n2];
    C_DCOPY(n*n,metric_,1,W,1);

    double* eigval = new double[n];
    int lwork = n * 3;
    double* work = new double[lwork];
    C_DSYEV('v','u',n,W,n,eigval,work,lwork);
    delete[] work;

    double * Jcopy = new double[n2];

    C_DCOPY(n*n,W,1,Jcopy,1);

    // Now form Jp^{-1/2} = U(T)*j'^{-1/2}*U,
    // where j'^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of J'
    double max_J = eigval[n-1];

    int nsig = 0;
    for (int ind=0; ind<n; ind++)
    {
        if (eigval[ind] / max_J < tol || eigval[ind] <= 0.0)
            eigval[ind] = 0.0;
        else
        {
            nsig++;
            eigval[ind] = 1.0 / sqrt(eigval[ind]);
        }
        // scale one set of eigenvectors by the diagonal elements j^{-1/2}
        C_DSCAL(n, eigval[ind], W + ind*n, 1);
    }
    delete[] eigval;

    C_DGEMM('T','N',n,n,n,1.0,Jcopy,n,W,n,0.0,metric_,n);

    delete [] Jcopy;
    delete [] W;
}

/*
void FittingMetric::form_cholesky_inverse()
{
    is_inverted_ = true;
    algorithm_ = "CHOLESKY";

    form_fitting_metric();

    pivot();
    for (int h = 0; h < metric_->nirrep(); h++)
    {

        if (metric_->colspi()[h] == 0) continue;

        // Cholesky Decomposition
        double** J = metric_->pointer(h);
        C_DPOTRF('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);
        for (int A = 0; A < metric_->colspi()[h]; A++)
            for (int B = 0; B < A; B++)
                J[A][B] = 0.0;
    }
    metric_->set_name("SO Basis Fitting Inverse (Cholesky)");
}
void FittingMetric::form_QR_inverse(double tol)
{
    is_inverted_ = true;
    algorithm_ = "QR";

    form_fitting_metric();

    pivot();
    for (int h = 0; h < metric_->nirrep(); h++)
    {

        if (metric_->colspi()[h] == 0) continue;

//        metric_->print();

        double** J = metric_->pointer(h);
        int n = metric_->colspi()[h];

        // Copy the J matrix to R (actually R')
        SharedMatrix R(new Matrix("R", n, n));
        double** Rp = R->pointer();
        C_DCOPY(n*(unsigned long int)n, J[0], 1, Rp[0], 1);

        // QR Decomposition
        double* tau = new double[n];

        // First, find out how much workspace to provide
        double work_size;
        C_DGEQRF(n,n,Rp[0],n,tau,&work_size, -1);

        // Now, do the QR decomposition
        int lwork = (int)work_size;
        double *work = new double[lwork];
        C_DGEQRF(n,n,Rp[0],n,tau,work,lwork);
        delete[] work;

        // Copy Jcopy to Q (actually Q')
        SharedMatrix Q(new Matrix("Q", n, n));
        double** Qp = Q->pointer();
        C_DCOPY(n*(unsigned long int)n, Rp[0], 1, Qp[0], 1);

        // Put R in the upper triangle where it belongs
        for (int i = 1; i < n; i++)
            for (int j = 0; j < i; j++)
            {
                Rp[j][i] = 0.0;
            }

        // First, find out how much workspace to provide
        C_DORGQR(n,n,n,Qp[0],n,tau,&work_size,-1);

        // Now, form Q
        lwork = (int)work_size;
        work = new double[lwork];
        C_DORGQR(n,n,n,Qp[0],n,tau,work,lwork);
        delete[] work;

        //Q->print();
        //R->print();

        // Find the number of significant basis functions
        int nsig = 0;
        double R_max = fabs(Rp[0][0]);
        for (int A = 0; A < n; A++)
        {
            if ((fabs(Rp[A][A]) / R_max) < tol)
                break;
            nsig++;
        }

        // Transform into the reduced basis
        // Just use R's memory, don't need it anymore
        C_DGEMM('N','N',nsig,n,n,1.0,Qp[0],n,J[0],n,0.0,Rp[0],n);
        C_DGEMM('N','T',nsig,nsig,n,1.0,Rp[0],n,Qp[0],n,0.0,J[0],nsig);

        // Find the Cholesky factor in the reduced basis
        C_DPOTRF('L',nsig,J[0],nsig);

        // Backsolve the triangular factor against the change of basis matrix
        C_DTRSM('L','U','N','N',nsig,n,1.0,J[0],nsig,Qp[0],n);

        // Zero out the metric
        memset(static_cast<void*>(J[0]), '\0', n*(unsigned long int)n);

        // Copy the top bit in
        C_DCOPY(n*(unsigned long int)nsig, Qp[0], 1, J[0], 1);

        delete[] tau;
    }
    metric_->set_name("SO Basis Fitting Inverse (QR)");
}

void FittingMetric::form_full_eig_inverse(double tol)
{
    is_inverted_ = true;
    algorithm_ = "EIG";

    form_fitting_metric();

    //metric_->print();

    for (int h = 0; h < metric_->nirrep(); h++)
    {

        if (metric_->colspi()[h] == 0) continue;

        double** J = metric_->pointer(h);
        int n = metric_->colspi()[h];

        // Copy J to W
        SharedMatrix W(new Matrix("W", n, n));
        double** Wp = W->pointer();
        C_DCOPY(n*(unsigned long int)n,J[0],1,Wp[0],1);

        double* eigval = new double[n];
        int lwork = n * 3;
        double* work = new double[lwork];
        C_DSYEV('v','u',n,Wp[0],n,eigval,work,lwork);
        delete[] work;

        SharedMatrix Jcopy(new Matrix("Jcopy", n, n));
        double** Jcopyp = Jcopy->pointer();

        C_DCOPY(n*(unsigned long int)n,Wp[0],1,Jcopyp[0],1);

        // Now form Jp^{-1/2} = U(T)*j'^{-1/2}*U,
        // where j'^{-1/2} is the diagonal matrix of the inverse square roots
        // of the eigenvalues, and U is the matrix of eigenvectors of J'
        double max_J = eigval[n-1];

        int nsig = 0;
        for (int ind=0; ind<n; ind++)
        {
            if (eigval[ind] / max_J < tol || eigval[ind] <= 0.0)
                eigval[ind] = 0.0;
            else
            {
                nsig++;
                eigval[ind] = 1.0 / eigval[ind];
            }
            // scale one set of eigenvectors by the diagonal elements j^{-1/2}
            C_DSCAL(n, eigval[ind], Wp[ind], 1);
        }
        delete[] eigval;

        C_DGEMM('T','N',n,n,n,1.0,Jcopyp[0],n,Wp[0],n,0.0,J[0],n);

    }
    metric_->set_name("SO Basis Fitting Inverse (Eig)");
}
void FittingMetric::form_full_inverse()
{
    is_inverted_ = true;
    algorithm_ = "FULL";

    form_fitting_metric();

    pivot();
    for (int h = 0; h < metric_->nirrep(); h++)
    {

        if (metric_->colspi()[h] == 0) continue;

        // Cholesky Decomposition
        double** J = metric_->pointer(h);
        int info = C_DPOTRF('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);
        // Inverse
        info = C_DPOTRI('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);

        for (int A = 0; A < metric_->colspi()[h]; A++)
            for (int B = 0; B < A; B++)
                J[A][B] = J[B][A];
    }
    metric_->set_name("SO Basis Fitting Inverse (Full)");
}
void FittingMetric::form_cholesky_factor()
{
    is_inverted_ = true;
    algorithm_ = "CHOLESKY";

    form_fitting_metric();

    //pivot();
    for (int h = 0; h < metric_->nirrep(); h++)
    {

        if (metric_->colspi()[h] == 0) continue;

        // Cholesky Decomposition
        double** J = metric_->pointer(h);
        C_DPOTRF('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);
    }
    metric_->set_name("SO Basis Cholesky Factor (Full)");
}

void FittingMetric::pivot()
{
    for (int h = 0; h < metric_->nirrep(); h++)
    {

        if (metric_->colspi()[h] == 0) continue;

        double** J = metric_->pointer(h);
        int* P = pivots_->pointer(h);
        int norbs = metric_->colspi()[h];
        double* Temp = new double[norbs];

        // Pivot
        double max;
        int Temp_p;
        int pivot;
        for (int i = 0; i<norbs-1; i++)
        {
            max = 0.0;
            //Where's the pivot diagonal?
            for (int j = i; j<norbs; j++)
                if (max <= fabs(J[j][j]))
                {
                    max = fabs(J[j][j]);
                    pivot = j;
                }

            //Rows
            C_DCOPY(norbs,&J[pivot][0],1,Temp,1);
            C_DCOPY(norbs,&J[i][0],1,&J[pivot][0],1);
            C_DCOPY(norbs,Temp,1,&J[i][0],1);

            //Columns
            C_DCOPY(norbs,&J[0][pivot],norbs,Temp,1);
            C_DCOPY(norbs,&J[0][i],norbs,&J[0][pivot],norbs);
            C_DCOPY(norbs,Temp,1,&J[0][i],norbs);

            Temp_p = P[i];
            P[i] = P[pivot];
            P[pivot] = Temp_p;
        }
        delete[] Temp;

        int* R = rev_pivots_->pointer(h);
        for (int i = 0; i < norbs; i++)
            R[P[i]] = i;
    }
}
*/

} // close namespace panache
