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

#ifndef PANACHE_SIMPLEMATRIX_H
#define PANACHE_SIMPLEMATRIX_H

#include <memory>

namespace panache {

class SimpleMatrix
{
private:
    double * data_;
    unsigned int nrow_, ncol_;

    void Delete_(void)
    {
        if(data_ != nullptr)
        {
            delete [] data_;
            data_ = nullptr;
        }
    }

public:
    Matrix(const Matrix & rhs) = delete;
    Matrix & operator=(const Matrix & rhs) = delete;

    unsigned int nrow(void) const
    {
        return nrow_;
    }

    unsigned int ncol(void) const
    {
        return ncol_;
    }

    Matrix(unsigned int nrow, unsigned int ncol)
        : nrow_(nrow),ncol_(ncol)
    {
        data_ = new double[nrow*ncol];
    }

    ~Matrix()
    {
        Delete_();
    }

    double & operator() (unsigned int i, unsigned int j)
    {
        return data_[i*ncol_ + j];
    }

    const double & operator() (unsigned int i, unsigned int j) const
    {
        return data_[i*ncol_ + j];
    }

    double * pointer(void)
    {
        return data_;
    }

    double const * pointer(void) const
    {
        return data_;
    }

    void zero(void)
    {
        for(size_t i = 0; i < nrow_*ncol_; i++)
            data_[i] = 0.0;
    }

};

} // close namespace Panache

#endif
