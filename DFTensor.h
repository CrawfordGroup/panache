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

    int nso_;
    int naux_;
    int nso2_;

    void OpenFile(void);
    void CloseFile(void);
    void ResetFile(void);

    std::unique_ptr<double[]> qso_;
    bool isinmem_;
    int curq_;

public:

    DFTensor(std::shared_ptr<BasisSet> primary,
             std::shared_ptr<BasisSet> auxiliary);

    ~DFTensor();

    int TensorDimensions(int & d1, int & d2, int & d3);

    // Calculate the Q matrix
    void GenQ(bool inmem);//, double * cmo, int nmo, cmo_is_trans);
    int GetBatch(double * mat, size_t size);


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

