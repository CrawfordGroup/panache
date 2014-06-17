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

#ifndef PANACHE_DFTENSOR_H
#define PANACHE_DFTENSOR_H

#include <fstream>

#include "Matrix.h"
#include "Timing.h"

#include <future>


namespace panache
{

class Molecule;
class BasisSet;

namespace reorder
{
class Orderings;
}

class DFTensor
{

protected:

    /// Debug level
    int debug_;
    /// Print level
    int print_;

    /// Molecule (fo convenience)
    std::shared_ptr<Molecule> molecule_;
    /// Primary basis set
    std::shared_ptr<BasisSet> primary_;
    /// Dealias basis set
    std::shared_ptr<BasisSet> auxiliary_;

    /// Symmetric inverse fitting metric
    SharedMatrix metric_;

    void common_init();
    void build_metric();
    void print_header();

    /// Filename for the matrix on disk
    std::string filename_;

    /// fstream object for the matrix on disk
    std::unique_ptr<std::fstream> matfile_;

    /// C matrix from the caller
    double * Cmo_;

    /// Whether the matrix is the transpose or not (think calling from FORTRAN)
    bool Cmo_trans_;

    /// Number of MO (# of columns of Cmo)
    int nmo_;
    int nmo2_;

    // Holds QC intermediate

    int nso_;
    int naux_;
    int nso2_;
    int nsotri_;

    void OpenFile(void);
    void CloseFile(void);
    void ResetFile(void);

    // buffers and placeholders for transformations
    std::unique_ptr<double[]> qso_; // holds the entire Q matrix
    std::unique_ptr<double[]> qc_;  // holds C(T) Q
    std::unique_ptr<double[]> q_;   // holds a couple Q (packed storage)
    std::unique_ptr<double[]> q_single_;   // holds a single Q (expanded storage)

    double * outbuffer_;
    size_t outbuffersize_;

    bool isinmem_; // Is Q completely in memory?
    int curq_;     // next Q to be (batch) read
    int curbatchq_;     // next Q to be (batch) read
    double * curbatch_; // Where the current batch is

    int GetBatch_Base(int ntoget);


    // Timing stuff
    timing::Timer timer_genqso;  // total time spent in GenQso
    timing::Timer timer_getbatch_qso;  // total time spent in GetBatch_Qso
    timing::Timer timer_getbatch_qmo;  // total time spent in GetBatch_Qmo


    // Threading
    int nthreads_;


    // Prefetching
    std::unique_ptr<double[]> q2_;
    std::future<int> fill_future_;

public:

    DFTensor(std::shared_ptr<BasisSet> primary,
             std::shared_ptr<BasisSet> auxiliary,
             const std::string & filename);

    ~DFTensor();

    int QsoDimensions(int & naux, int & nso2);

    // Calculate the Q matrix
    void GenQso(bool inmem);

    void SetCMatrix(double * cmo, int nmo, bool cmo_is_trans);

    void SetOutputBuffer(double * buf, size_t size);
    int GetBatch_Qso(void);
    int GetBatch_Qmo(void);

    void ResetBatches(void);


/*
    int CalculateERI(double * qso, int qsosize, int shell1, int shell2, int shell3, int shell4, double * outbuffer, int buffersize);

    int CalculateERIMulti(double * qso, int qsosize,
                          int shell1, int nshell1,
                          int shell2, int nshell2,
                          int shell3, int nshell3,
                          int shell4, int nshell4,
                          double * outbuffer, int buffersize);

    void ReorderQ(double * qso, int qsosize, const reorder::Orderings & order);
    void ReorderQ_GAMESS(double * qso, int qsosize);
*/
};

}

#endif //PANACHE_DFTENSOR_H

