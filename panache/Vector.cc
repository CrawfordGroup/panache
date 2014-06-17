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

#include <cstring>
#include "Lapack.h"
#include <algorithm>
#include <numeric>
#include "Matrix.h"
#include "Vector.h"
#include "Output.h"

namespace panache {

Vector::Vector()
{
    dimpi_ = NULL;
    nirrep_ = 0;
    name_ = "";
}

Vector::Vector(const Vector& c)
{
    nirrep_ = c.nirrep_;
    dimpi_ = c.dimpi_;
    alloc();
    copy_from(c);
    name_ = c.name_;
}

Vector::Vector(int nirreps, int *dimpi)
    : dimpi_(nirreps)
{
    nirrep_ = nirreps;
    dimpi_ = dimpi;
    alloc();
}

Vector::Vector(int dim)
    : dimpi_(1)
{
    nirrep_ = 1;
    dimpi_[0] = dim;
    alloc();
}

Vector::Vector(const std::string& name, int nirreps, int *dimpi)
    : dimpi_(nirreps)
{
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
    name_ = name;
}

Vector::Vector(const std::string& name, int dim)
    : dimpi_(1)
{
    nirrep_ = 1;
    dimpi_[0] = dim;
    alloc();
    name_ = name;
}

Vector::Vector(const Dimension& v)
{
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
    name_ = v.name();
}

Vector::Vector(const std::string& name, const Dimension& v)
{
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
    name_ = name;
}

Vector::~Vector()
{
    release();
}

SharedVector Vector::create(const std::string &name, const Dimension &dim)
{
    return SharedVector(new Vector(name, dim));
}

void Vector::init(int nirreps, int *dimpi)
{
    dimpi_.init("", nirreps);
    nirrep_ = nirreps;
    dimpi_ = dimpi;
    alloc();
}

void Vector::init(int nirreps, const int *dimpi, const std::string& name)
{
    name_ = name;
    dimpi_.init("", nirreps);
    dimpi_ = dimpi;
    alloc();
}

void Vector::init(const Dimension &v)
{
    name_ = v.name();
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
}

Vector* Vector::clone()
{
    Vector *temp = new Vector(nirrep_, dimpi_);
    temp->copy(this);
    return temp;
}

void Vector::alloc()
{
    if (vector_.size())
        release();

    int total = dimpi_.sum();
    v_.resize(total);

    std::fill(vector_.begin(), vector_.end(), (double*)0);
    std::fill(v_.begin(), v_.end(), 0.0);

    assign_pointer_offsets();
}

void Vector::assign_pointer_offsets()
{
    // Resize just to be sure it's the correct size
    vector_.resize(dimpi_.n(), 0);

    size_t offset = 0;
    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h])
            vector_[h] = &(v_[0]) + offset;
        else
            vector_[h] = NULL;
        offset += dimpi_[h];
    }
}

void Vector::release()
{
    std::fill(vector_.begin(), vector_.end(), (double*)0);
    std::fill(v_.begin(), v_.end(), 0.0);
}

void Vector::copy_from(const Vector& other)
{
    nirrep_ = other.dimpi_.n();
    dimpi_ = other.dimpi_;
    v_ = other.v_;
    assign_pointer_offsets();
}

void Vector::copy(const Vector *rhs)
{
    copy_from(*rhs);
}

void Vector::copy(const Vector &rhs)
{
    copy_from(rhs);
}

void Vector::set(double *vec)
{
    std::copy(vec, vec + dimpi_.sum(), v_.begin());
}

void Vector::zero()
{
    std::fill(v_.begin(), v_.end(), 0.0);
}

double *Vector::to_block_vector() {
    size_t size=0;
    for (int h=0; h<nirrep_; ++h)
        size += dimpi_[h];

    double *temp = new double[size];
    size_t offset = 0;
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<dimpi_[h]; ++i) {
            temp[i+offset] = vector_[h][i];
        }
        offset += dimpi_[h];
    }

    return temp;
}

void Vector::gemv(bool transa, double alpha, Matrix* A, Vector* X, double beta)
{
    char trans = transa ? 't' : 'n';

    for (int h =0; h < nirrep_; ++h) {
        C_DGEMV(trans, A->rowspi_[h], A->colspi_[h], alpha, &(A->matrix_[h][0][0]),
                A->rowspi_[h], &(X->vector_[h][0]), 1, beta, &(vector_[h][0]), 1);
    }
}

double Vector::dot(Vector* X)
{
    double tmp = 0.0;

    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h] != X->dimpi_[h]) {
           output::printf("Vector::dot: Vectors are not of the same size.\n");
            return 0.0;
        }
        for (int i=0; i<dimpi_[h]; ++i) {
            tmp += vector_[h][i] * X->vector_[h][i];
        }
    }

    return tmp;
}

double Vector::norm()
{
    double tmp = 0.0;

    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<dimpi_[h]; ++i) {
            tmp += vector_[h][i] * vector_[h][i];
        }
    }

    return tmp;
}

template<class T>
struct scale_vector {
    typedef T data_type;
    typedef data_type result_type;

    const data_type scalar;
    scale_vector(const data_type& sc) : scalar(sc) {}

    result_type operator()(data_type value) const { value *= scalar; return value; }
};

void Vector::scale(const double& sc)
{
    std::transform(v_.begin(), v_.end(), v_.begin(), scale_vector<double>(sc));
}

void Vector::add(const std::vector<double>& rhs)
{
    size_t min = std::min(rhs.size(), v_.size());
    for (size_t i=0; i<min; ++i)
        v_[i] += rhs[i];
}

IntVector::IntVector() {
    vector_ = NULL;
    dimpi_ = NULL;
    nirrep_ = 0;
    name_ = "";
}

IntVector::IntVector(const IntVector& c) {
    vector_ = NULL;
    nirrep_ = c.nirrep_;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = c.dimpi_[h];
    alloc();
    copy_from(c.vector_);
    name_ = c.name_;
}

IntVector::IntVector(int nirreps, int *dimpi) {
    vector_ = NULL;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
}
IntVector::IntVector(int dim) {
    vector_ = NULL;
    nirrep_ = 1;
    dimpi_ = new int[nirrep_];
    dimpi_[0] = dim;
    alloc();
}
IntVector::IntVector(const std::string& name, int nirreps, int *dimpi) {
    vector_ = NULL;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
    name_ = name;
}
IntVector::IntVector(const std::string& name, int dim) {
    vector_ = NULL;
    nirrep_ = 1;
    dimpi_ = new int[nirrep_];
    dimpi_[0] = dim;
    alloc();
    name_ = name;
}

IntVector::~IntVector() {
    release();
    if (dimpi_)
        delete[] dimpi_;
}

void IntVector::init(int nirreps, int *dimpi)
{
    if (dimpi_) delete[] dimpi_;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
}

void IntVector::alloc() {
    if (vector_)
        release();

    vector_ = (int**)malloc(sizeof(int*) * nirrep_);
    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h]) {
            vector_[h] = new int[dimpi_[h]];
            memset(vector_[h], 0, dimpi_[h] * sizeof(int));
        }
    }
}

void IntVector::release() {
    if (!vector_)
        return;

    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h])
            delete[] (vector_[h]);
    }
    free(vector_);
    vector_ = NULL;
}

void IntVector::copy_from(int **c) {
    size_t size;
    for (int h=0; h<nirrep_; ++h) {
        size = dimpi_[h] * sizeof(int);
        if (size)
            memcpy(&(vector_[h][0]), &(c[h][0]), size);
    }
}

void IntVector::copy(const IntVector *rhs) {
    if (nirrep_ != rhs->nirrep_) {
        release();
        if (dimpi_)
            delete[] dimpi_;
        nirrep_ = rhs->nirrep_;
        dimpi_ = new int[nirrep_];
        for (int h=0; h<nirrep_; ++h)
            dimpi_[h] = rhs->dimpi_[h];
        alloc();
    }
    copy_from(rhs->vector_);
}
void IntVector::copy(const IntVector &rhs) {
    if (nirrep_ != rhs.nirrep_) {
        release();
        if (dimpi_)
            delete[] dimpi_;
        nirrep_ = rhs.nirrep_;
        dimpi_ = new int[nirrep_];
        for (int h=0; h<nirrep_; ++h)
            dimpi_[h] = rhs.dimpi_[h];
        alloc();
    }
    copy_from(rhs.vector_);
}

void IntVector::set(int *vec) {
    int h, i, ij;

    ij = 0;
    for (h=0; h<nirrep_; ++h) {
        for (i=0; i<dimpi_[h]; ++i) {
            vector_[h][i] = vec[ij++];
        }
    }
}

int *IntVector::to_block_vector() {
    size_t size=0;
    for (int h=0; h<nirrep_; ++h)
        size += dimpi_[h];

    int *temp = new int[size];
    size_t offset = 0;
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<dimpi_[h]; ++i) {
            temp[i+offset] = vector_[h][i];
        }
        offset += dimpi_[h];
    }

    return temp;
}

} // end namespace panache
