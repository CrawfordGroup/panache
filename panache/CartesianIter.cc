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

#include "panache/CartesianIter.h"
#include "panache/Exception.h"

namespace panache {

/*****************
 * CartesianIter *
 *****************/

CartesianIter::CartesianIter(int l) :
     exponents_({{0,0,0}}),l_(l), bfn_(0)
{

}


CartesianIter::~CartesianIter()
{

}

void CartesianIter::start()
{
    bfn_ = exponents_[1] = exponents_[2] = 0;
    exponents_[0] = l_;
}


void CartesianIter::next()
{
    if (exponents_[2] < l_ - exponents_[0]) {
        (exponents_[1])--;
        (exponents_[2])++;
    }
    else {
        (exponents_[0])--;
        exponents_[2] = 0;
        exponents_[1] = l_ - exponents_[0];
    }
    bfn_++;
}

CartesianIter::operator int()
{
    return (exponents_[0] >= 0);
}


/**************************
 * RedundantCartesianIter *
 **************************/

RedundantCartesianIter::RedundantCartesianIter(int l) :
    done_(0), l_(l), axis_(0)
{
    l_ = l;
    axis_ = new int[l_];
}


RedundantCartesianIter::~RedundantCartesianIter()
{
    delete[] axis_;
}


int RedundantCartesianIter::bfn()
{
    int i = a();
    int am = l();
    if (am == i)
        return 0;
    else {
        int j = b();
        int c = am - i;
        return ((((c+1)*c)>>1)+c-j);
    }
}


/*****************************
 * RedundantCartesianSubIter *
 *****************************/

RedundantCartesianSubIter::RedundantCartesianSubIter(int l)
{
    l_ = l;
    axis_ = new int[l_];
    zloc_ = new int[l_];
    yloc_ = new int[l_];
}


RedundantCartesianSubIter::~RedundantCartesianSubIter()
{
    delete[] axis_;
    delete[] zloc_;
    delete[] yloc_;
}


void RedundantCartesianSubIter::start(int a, int b, int c)
{
    if (l_ != a + b + c) {
        throw RuntimeError("RedundantCartesianSubIter::start: bad args");
    }

    if (l_==0) {
        done_ = 1;
        return;
    } else {
        done_ = 0;
    }

    e_[0] = a;
    e_[1] = b;
    e_[2] = c;

    int ii=0;
    for (int i=0; i<c; i++,ii++) {
        axis_[ii] = 2;
        zloc_[i] = c-i-1;
    }
    for (int i=0; i<b; i++,ii++) {
        axis_[ii] = 1;
        yloc_[i] = b-i-1;
    }
    for (int i=0; i<a; i++,ii++) axis_[ii] = 0;
}


bool RedundantCartesianSubIter::advance(int l, int *loc, int n)
{
    int maxloc = l-1;
    for (int i=0; i<n; i++) {
        if (loc[i] < maxloc) {
            loc[i]++;
            for (int j=i-1; j>=0; j--) loc[j] = loc[j+1] + 1;
            return true;
        }
        else {
            maxloc = loc[i]-1;
        }
    }
    return false;
}


// This loops through all unique axis vectors that have a
// given total a, b, and c.  It is done by looping through
// all possible positions for z, then y, leaving x to be
// filled in.
void RedundantCartesianSubIter::next()
{
    int currentz = 0;
    int currenty = 0;
    int nz = c();
    int ny = b();

    if (!advance(l(),zloc_,nz)) {
        if (!advance(l()-nz,yloc_,ny)) {
            done_ = 1;
            return;
        }
        else {
            for (int i=0; i<nz; i++) {
                zloc_[i] = nz-i-1;
            }
        }
    }

    int nonz = l()-nz-1;
    for (int i = l()-1; i>=0; i--) {
        if (currentz<nz && zloc_[currentz]==i) {
            axis_[i] = 2;
            currentz++;
        }
        else if (currenty<ny && yloc_[currenty]==nonz) {
            axis_[i] = 1;
            currenty++;
            nonz--;
        }
        else {
            axis_[i] = 0;
            nonz--;
        }
    }
}

int RedundantCartesianSubIter::valid()
{
    int t[3];
    int i;

    for (i=0; i<3; i++)
        t[i] = 0;

    for (i=0; i<l_; i++)
        t[axis_[i]]++;

    return t[0] == e_[0] && t[1] == e_[1] && t[2] == e_[2];
}

int RedundantCartesianSubIter::bfn()
{
    int i = a();
    int am = l();
    if (am == i)
        return 0;
    else {
        int j = b();
        int c = am - i;
        return ((((c+1)*c)>>1)+c-j);
    }
}

void RedundantCartesianIter::start()
{
    if (l_==0)
        done_ = 1;
    else
        done_ = 0;

    for (int i=0; i<l_; i++)
        axis_[i] = 0;
}

void RedundantCartesianIter::next()
{
    for (int i=0; i<l_; i++) {
        if (axis_[i] == 2)
            axis_[i] = 0;
        else {
            axis_[i]++;
            return;
        }
    }
    done_ = 1;
}

int RedundantCartesianIter::l(int axis) const
{
    int i;
    int r = 0;
    for (i=0; i<l_; i++) if (axis_[i]==axis) r++;
    return r;
}

int RedundantCartesianIter::a() const
{
    return l(0);
}

int RedundantCartesianIter::b() const
{
    return l(1);
}

int RedundantCartesianIter::c() const
{
    return l(2);
}

} // close namespace

