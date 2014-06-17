#ifndef PANACHE_POINTGROUP_H
#define PANACHE_POINTGROUP_H

#include <cstdint>
#include <vector>
#include <map>
#include <memory>

#include "Vector3.h"
#include "CharacterTable.h"

using std::shared_ptr;

namespace panache {

#define NUM_TO_OPERATOR_ID(x) ((x) ? 1<<((x)-1) : 0)
#define SKIP_THIS_OPERATOR(num,bit) ((bit) ? !((1<<((bit)-1)) & (num)) : 0 )


// ///////////////////////////////////////////////////////////
/** The PointGroup class is really a place holder for a CharacterTable.  It
 contains a string representation of the Schoenflies symbol of a point
 group, a frame of reference for the symmetry operation transformation
 matrices, and a point of origin.  The origin is not respected by the
 symmetry operations, so if you want to use a point group with a nonzero
 origin, first translate all your coordinates to the origin and then set
 the origin to zero.  */
class PointGroup {
  private:
    std::string symb;
    Vector3 origin_;
    unsigned char bits_;

  public:
    PointGroup();

    // These 2 constructors do not work right now.
    /** This constructor takes a string containing the Schoenflies symbol
        of the point group as its only argument. */
    PointGroup(const std::string&);
    /** Like the above, but this constructor also takes a point of origin
        as an argument. */
    PointGroup(const std::string&,const Vector3&);

    /** Using the bitwise representation constructor the point group object. */
    PointGroup(unsigned char bits);
    /** Using the bitwise representation constructor the point group object. */
    PointGroup(unsigned char bits, const Vector3&);

    /** The PointGroup KeyVal constructor looks for three keywords:
       symmetry, symmetry_frame, and origin. symmetry is a string
       containing the Schoenflies symbol of the point group.  origin is an
       array of doubles which gives the x, y, and z coordinates of the
       origin of the symmetry frame.  symmetry_frame is a 3 by 3 array of
       arrays of doubles which specify the principal axes for the
       transformation matrices as a unitary rotation.

       For example, a simple input which will use the default origin and
       symmetry_frame ((0,0,0) and the unit matrix, respectively), might
       look like this:

       <pre>
       pointgrp<PointGroup>: (
         symmetry = "c2v"
       )
       </pre>

       By default, the principal rotation axis is taken to be the z axis.
       If you already have a set of coordinates which assume that the
       rotation axis is the x axis, then you'll have to rotate your frame
       of reference with symmetry_frame:

       <pre>
       pointgrp<PointGroup>: (
         symmetry = "c2v"
         symmetry_frame = [
           [ 0 0 1 ]
           [ 0 1 0 ]
           [ 1 0 0 ]
         ]
       )
       </pre>
     */
    // PointGroup(const Ref<KeyVal>&);

    // PointGroup(StateIn&);
    PointGroup(const PointGroup&);
    PointGroup(const shared_ptr<PointGroup>&);
    ~PointGroup();

    PointGroup& operator=(const PointGroup&);

    /// Returns 1 if the point groups are equivalent, 0 otherwise.
    int equiv(const shared_ptr<PointGroup> &, double tol = 1.0e-6) const;

    /// Returns the CharacterTable for this point group.
    CharacterTable char_table() const;
    /// Returns the Schoenflies symbol for this point group.
    std::string symbol() const { return symb; }
    /// Returns the frame of reference for this point group.
    /// Returns the origin of the symmetry frame.
    Vector3& origin() { return origin_; }
    const Vector3& origin() const { return origin_; }

    /// Sets (or resets) the Schoenflies symbol.
    void set_symbol(const std::string&);

    /// Returns the bitwise representation of the point group
    unsigned char bits() const { return bits_; }

    //void print(FILE *out = outfile) const;
};

}

#endif //PANACHE_POINTGROUP_H
