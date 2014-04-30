#ifndef PANACHE_SYMREP_H
#define PANACHE_SYMREP_H

#include <cstring> // for memset

namespace panache {

class SymmetryOperation;

// //////////////////////////////////////////////////////////////////

/** The SymRep class provides an n dimensional matrix representation of a
    symmetry operation, such as a rotation or reflection.  The trace of a
    SymRep can be used as the character for that symmetry operation.  d is
    hardwired to 5x5 since the H irrep in Ih is 5 dimensional.
*/
class SymRep {
  private:
    int n;
    double d[5][5];

  public:
    SymRep(int =0);
    SymRep(const SymmetryOperation&);
    ~SymRep();

    /// Cast to a SymmetryOperation.
    operator SymmetryOperation() const;

    /// returns the trace of the transformation matrix
    double trace() const
    {
      double r=0;
      for (int i=0; i < n; i++)
        r += d[i][i];
      return r;
    }


    /// set the dimension of d
    void set_dim(int i) { n=i; }

    /// returns the i'th row of the transformation matrix
    double* operator[](int i) { return d[i]; }
    /// const version of the above
    const double* operator[](int i) const { return d[i]; }

    /** returns a reference to the (i,j)th element of the transformation
        matrix */
    double& operator()(int i, int j) { return d[i][j]; }
    /// const version of double& operator()(int i, int j)
    double operator()(int i, int j) const { return d[i][j]; }

    /// zero out the symop
    void zero() { memset(d,0,sizeof(double)*25); }

    /// This operates on this with r (i.e. return r * this).
    SymRep operate(const SymRep& r) const;

    /// This performs the transform r * this * r~
    SymRep transform(const SymRep& r) const;

    /// Set equal to a unit matrix
    void unit() {
      zero(); d[0][0] = d[1][1] = d[2][2] = d[3][3] = d[4][4] = 1.0;
    }

    /// Set equal to the identity
    void E() { unit(); }

    /// Set equal to an inversion
    void i() { zero(); d[0][0] = d[1][1] = d[2][2] = d[3][3] = d[4][4] = -1.0;}

    /// Set equal to reflection in xy plane
    void sigma_h();

    /// Set equal to reflection in xz plane
    void sigma_xz();

    /// Set equal to reflection in yz plane
    void sigma_yz();

    /// Set equal to a clockwise rotation by 2pi/n
    void rotation(int n);
    void rotation(double theta);

    /// Set equal to C2 about the x axis
    void c2_x();

    /// Set equal to C2 about the y axis
    void c2_y();

    /// Set equal to C2 about the z axis
    void c2_z();

    /// print the matrix
    // void print(std::ostream& =ExEnv::out0()) const;
};

} // close namespace panache
// //////////////////////////////////////////////////////////////////

#endif //PANACHE_SYMREP_H
