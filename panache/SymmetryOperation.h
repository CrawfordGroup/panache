#ifndef PANACHE_SYMMETRYOPERATION_H
#define PANACHE_SYMMETRYOPERATION_H

#include <cstring> // for memset
#include "SymmOps.h"

namespace panache {

// //////////////////////////////////////////////////////////////////

/** The SymmetryOperation class provides a 3 by 3 matrix
    representation of a symmetry operation, such as a rotation or reflection.
*/
class SymmetryOperation {
  private:
    double d[3][3];
    unsigned short bits_;

    void analyze_d();

  public:
    SymmetryOperation();
    SymmetryOperation(const SymmetryOperation &);
    ~SymmetryOperation();

    SymmetryOperation & operator = (SymmetryOperation const & a); // Assignment operator

    /// returns the trace of the transformation matrix
    double trace() const { return d[0][0]+d[1][1]+d[2][2]; }

    /// returns the i'th row of the transformation matrix
    double* operator[](int i) { return d[i]; }

    /// const version of the above
    const double* operator[](int i) const { return d[i]; }

    /** returns a reference to the (i,j)th element of the transformation
        matrix */
    double& operator()(int i, int j) { return d[i][j]; }

    /// const version of the above
    double operator()(int i, int j) const { return d[i][j]; }

    /// zero out the symop
    void zero() { memset(d,0,sizeof(double)*9); }

    /// This operates on this with r (i.e. return r * this).
    SymmetryOperation operate(const SymmetryOperation& r) const;

    /// This performs the transform r * this * r~
    SymmetryOperation transform(const SymmetryOperation& r) const;

    /// Get the bit value.
    unsigned char bit() const { return bits_; }

    /// Set equal to a unit matrix
    void unit() { zero(); d[0][0] = d[1][1] = d[2][2] = 1.0;}

    /// Set equal to E
    void E() { unit(); bits_ = SymmOps::E; }

    /// Set equal to an inversion
    void i() { zero(); d[0][0] = d[1][1] = d[2][2] = -1.0; bits_ = SymmOps::i; }

    /// Set equal to reflection in xy plane
    void sigma_xy() { unit(); d[2][2] = -1.0; bits_ = SymmOps::Sigma_xy; }

    /// Set equal to reflection in xz plane
    void sigma_xz() { unit(); d[1][1] = -1.0; bits_ = SymmOps::Sigma_xz; }

    /// Set equal to reflection in yz plane
    void sigma_yz() { unit(); d[0][0] = -1.0; bits_ = SymmOps::Sigma_yz; }

    /// Set equal to a clockwise rotation by 2pi/n
    void rotation(int n);
    void rotation(double theta);

    /// Set equal to C2 about the x axis
    void c2_x() { i(); d[0][0] = 1.0; bits_ = SymmOps::C2_x; }

    /// Set equal to C2 about the y axis
    void c2_y() { i(); d[1][1] = 1.0; bits_ = SymmOps::C2_y; }

    /// Set equal to C2 about the z axis
    void c2_z() { i(); d[2][2] = 1.0; bits_ = SymmOps::C2_z; }

    void transpose();
};

} // end namespace panache

#endif //PANACHE_SYMMETRYOPERATION_H
