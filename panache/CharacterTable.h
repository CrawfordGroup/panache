#ifndef PANACHE_CHARACTERTABLE_H
#define PANACHE_CHARACTERTABLE_H

#include "IrreducibleRepresentation.h"
#include "SymmetryOperation.h"

using std::string;

namespace panache {


// ///////////////////////////////////////////////////////////
/** The CharacterTable class provides a workable character table
 for all of the non-cubic point groups.  While I have tried to match the
 ordering in Cotton's book, I don't guarantee that it is always followed.
 It shouldn't matter anyway.  Also note that I don't lump symmetry
 operations of the same class together.  For example, in C3v there are two
 distinct C3 rotations and 3 distinct reflections, each with a separate
 character.  Thus symop has 6 elements rather than the 3 you'll find in
 most published character tables. */
class CharacterTable {
    int nt;                              //< order of the princ rot axis
    PointGroups::Groups pg;              //< the class of the point group
    int nirrep_;                         //< the number of irreps in this pg
    IrreducibleRepresentation *gamma_;   //< an array of irreps
    SymmetryOperation *symop;            //< the matrices describing sym ops
    int *_inv;                           //< index of the inverse symop
    string symb;                    //< the Schoenflies symbol for the pg
    unsigned char bits_;                 //< Bitwise representation of the symmetry operations

    /// this fills in the irrep and symop arrays.
    int make_table();
    void common_init();

  public:
    CharacterTable();
    /** This constructor takes the Schoenflies symbol of a point group as
        input. */
    CharacterTable(const string&);
    /** This constructor takes the bitswise representation of a point group as
        input. */
    CharacterTable(unsigned char);

    CharacterTable(const CharacterTable&);
    ~CharacterTable();

    CharacterTable& operator=(const CharacterTable&);

    /// Returns the number of irreps.
    int nirrep() const { return nirrep_; }
    /// Returns the order of the point group
    int order() const { return nirrep_; }
    /// Returns the Schoenflies symbol for the point group
    const string& symbol() const { return symb; }
    /// Returns the i'th irrep.
    IrreducibleRepresentation& gamma(int i) { return gamma_[i]; }
    /// Returns the i'th symmetry operation.
    SymmetryOperation& symm_operation(int i) { return symop[i]; }

    /** Cn, Cnh, Sn, T, and Th point groups have complex representations.
        This function returns 1 if the point group has a complex
        representation, 0 otherwise. */
    int complex() const {
      return 0;
    }

    /// Returns the index of the symop which is the inverse of symop[i].
    int inverse(int i) const { return _inv[i]; }

    int ncomp() const {
      int ret=0;
      for (int i=0; i < nirrep_; i++) {
        int nc = (gamma_[i].complex()) ? 1 : gamma_[i].degen;
        ret += nc;
      }
      return ret;
    }

    /// Returns the irrep component i belongs to.
    int which_irrep(int i) {
      for (int ir=0, cn=0; ir < nirrep_; ir++) {
        int nc = (gamma_[ir].complex()) ? 1 : gamma_[ir].degen;
        for (int c=0; c < nc; c++,cn++)
          if (cn==i)
            return ir;
      }
      return -1;
    }

    /// Returns which component i is.
    int which_comp(int i) {
      for (int ir=0, cn=0; ir < nirrep_; ir++) {
        int nc = (gamma_[ir].complex()) ? 1 : gamma_[ir].degen;
        for (int c=0; c < nc; c++,cn++)
          if (cn==i)
            return c;
      }
      return -1;
    }

    unsigned char bits();

    /// This prints the irrep to the given file, or stdout if none is given.
    // void print(FILE *out=outfile) const;
};

}

#endif //PANACHE_CHARACTERTABLE_H
