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

#include "Matrix.h"

namespace panache {

class Molecule;
class BasisSet;

class DFTensor {

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

    /// Number of AO primary functions
    int nso_;
    /// Number of MO primary functions
    int nmo_;

    /// Number of grid points
    int naux_;

    void common_init();
    void build_metric();
    void print_header();

public:

    DFTensor(std::shared_ptr<BasisSet> primary,
                   std::shared_ptr<BasisSet> auxiliary,
                   int nso, int nmo); 
    ~DFTensor();

    SharedMatrix Qso();
};

}

#endif //PANACHE_DFTENSOR_H
