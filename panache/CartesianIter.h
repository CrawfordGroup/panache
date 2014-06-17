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

#ifndef PANACHE_CARTESIANITER_H
#define PANACHE_CARTESIANITER_H

#include <array>

namespace panache {

/** CartesianIter gives the ordering of the Cartesian functions
    that is used in PSI4. */
class CartesianIter
{
private:
    std::array<int, 3> exponents_;
    int l_;
    int bfn_;

public:
    /// Initialize the iterator for the given angular momentum.
    CartesianIter(int l);
    ~CartesianIter();

    /// Start the iteration.
    virtual void start();
    /// Move to the next Cartesian function.
    virtual void next();
    /// Returns nonzero if the iterator currently holds valid data.
    virtual operator int();

    /// Returns the number of Cartesian functions.
    int n() const {
        return ((l_>=0)?((((l_)+2)*((l_)+1))>>1):0);
    }
    /// Returns the x exponent
    int a() const {
        return exponents_[0];
    }
    /// Returns the y exponent
    int b() const {
        return exponents_[1];
    }
    /// Returns the z exponent
    int c() const {
        return exponents_[2];
    }
    /// Return the angular momentum
    int l() const {
        return l_;
    }
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i) const {
        return exponents_[i]; 
    }
    /** Returns the number of the current basis function within the shell.
        This starts at 0 and sequentially increases as next() is called. */
    int bfn() {
        return bfn_;
    }
    std::array<int, 3> exponents() { return exponents_; }
};


/** RedundantCartesianIter objects loop through all possible combinations
    of a given number of axes.  This is used to compute the transformation
    matrices that maps a set of Cartesian functions to another set of
    Cartesian functions in a rotated coordinate system. */
class RedundantCartesianIter {
private:
    int done_;
    int l_;
    int *axis_;

public:
    /// Create a object for the given angular momentum.
    RedundantCartesianIter(int l);
    virtual ~RedundantCartesianIter();

    /// Return the current Cartesian basis function number.
    virtual int bfn();

    /// Initialize the iterator.
    void start();
    /// Move to the next combination of axes.
    void next();
    /// Returns nonzero if the iterator currently hold valid data.
    operator int() {
        return !done_;
    }

    /// The current exponent of x.
    int a() const;
    /// The current exponent of y.
    int b() const;
    /// The current exponent of z.
    int c() const;
    /// The angular momentum.
    int l() const {
        return l_;
    }
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i) const;
    /// Return the i'th axis.
    int axis(int i) const {
        return axis_[i];
    }
};


/** Like RedundantCartesianIter, except a, b, and c are fixed to a given
    value. */
class RedundantCartesianSubIter {
private:
    int done_;
    int l_;
    int e_[3];
    int *axis_;

    // the locations of the z's in the axis array
    int *zloc_;
    // the locations of the y's in the subarray after the z's are removed
    int *yloc_;
    int valid();

    static bool advance(int l, int *loc, int n);

public:
    /// Create a object for the given angular momentum.
    RedundantCartesianSubIter(int l);
    virtual ~RedundantCartesianSubIter();

    /// Return the current Cartesian basis function number.
    virtual int bfn();

    /** Initialize the iterator.  The constraints on a, b, and c are
        given as arguments. */
    void start(int a, int b, int c);
    /// Move to the next combination of axes.
    void next();
    /// Returns nonzero if the iterator currently hold valid data.
    operator int() const {
        return !done_;
    }

    /// The current exponent of x.
    int a() const {
        return e_[0];
    }
    /// The current exponent of y.
    int b() const {
        return e_[1];
    }
    /// The current exponent of z.
    int c() const {
        return e_[2];
    }
    /// The angular momentum.
    int l() const {
        return l_;
    }
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i) {
        return e_[i];
    }
    /// Return the i'th axis.
    int axis(int i) {
        return axis_[i];
    }
};

}

#endif //PANACHE_CARTESIANITER_H
