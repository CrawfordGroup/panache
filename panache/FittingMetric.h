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

#ifndef PANACHE_FITTINGMETRIC_H
#define PANACHE_FITTINGMETRIC_H

#include <vector>

namespace panache {

class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;
class Molecule;

class FittingMetric {

protected:
    /// Pointer to the auxiliary basis set
    SharedBasisSet aux_;

    /// Number of auxiliary basis functions
    int naux_;

    /// Pointer to the poisson basis set
    SharedBasisSet pois_;
    /// Is the metric poisson?
    bool is_poisson_;

    /// Is the metric inverted or just a J matrix?
    bool is_inverted_;

    /// Range separation omega (0.0 if not used)
    double omega_;

    /// The fitting metric or symmetric inverse
    double * metric_;

    int nthreads_;  //!< Number of threads to use

    /// The indices (per irrep) of pivots
    std::vector<int> pivots_;
    /// The indices (per irrep) of reverse pivots
    std::vector<int> rev_pivots_;

    /// The fitting algorithm selected
    std::string algorithm_;

    /// Fully pivot the fitting metric
    void pivot();

public:
    FittingMetric(const FittingMetric & f) = delete;
    FittingMetric(const FittingMetric && f) = delete;
    FittingMetric & operator=(const FittingMetric && f) = delete;

    /// DF Fitting Metric
    FittingMetric(SharedBasisSet aux, int nthreads);
    // DF Fitting Metric
    FittingMetric(SharedBasisSet aux, double omega, int nthreads);
    /// Poisson Fitting Metric
    //FittingMetric(SharedBasisSet aux, SharedBasisSet pois, bool force_C1 = false);

    /// Destructor
    ~FittingMetric();

    /// What algorithm to use for symmetric inverse?
    std::string get_algorithm() const {return algorithm_; }
    /// Are poisson functions used?
    bool is_poisson() const {return is_poisson_; }
    /// Is the metric inverted?
    bool is_inverted() const {return is_inverted_; }

    /// The fitting metric or symmetric inverse
    double * get_metric() const {return metric_; }

    /// The vector of pivots (for stability) (pivoted->global)
    std::vector<int> get_pivots() const {return pivots_; }

    /// The vector of back pivots (for stability) (global->pivoted)
    std::vector<int> get_reverse_pivots() const {return rev_pivots_; }

    /// The gaussian fitting basis
    SharedBasisSet get_auxiliary_basis() const {return aux_; }
    /// The poisson fitting basis
    SharedBasisSet get_poisson_basis() const {return pois_; }

    /// Build the raw fitting metric (sets up indices to canonical)
    void form_fitting_metric();
    /// Build the eigendecomposed half inverse metric (calls form_fitting_metric)
    void form_eig_inverse(double tol = 1.0E-10);

/*
    /// Build the Cholesky half inverse metric (calls form_fitting_metric)
    void form_cholesky_inverse();
    /// Build the QR half inverse metric (calls form_fitting_metric)
    void form_QR_inverse(double tol = 1.0E-10);
    /// Build the full inverse metric. NOT RECOMMENDED: Numerical stability (calls form_fitting_metric)
    void form_full_inverse();
    /// Build the full inverse metric.
    void form_full_eig_inverse(double tol = 1.0E-10);
    /// Build the full metric's Cholesky factor. RECOMMENDED: Numerical stability
    void form_cholesky_factor();
*/
};

/*!
 * \brief A shared FittingMetric object
 */
typedef std::shared_ptr<FittingMetric> SharedFittingMetric;

}

#endif //PANACHE_FITTINGMETRIC_H
