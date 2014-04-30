#ifndef PANACHE_IRREDUCIBLEREPRESENTATION_H
#define PANACHE_IRREDUCIBLEREPRESENTATION_H

#include "SymRep.h"

namespace panache {

class CharacterTable;


/** The IrreducibleRepresentation class provides information associated
    with a particular irreducible representation of a point group.  This
    includes the Mulliken symbol for the irrep, the degeneracy of the
    irrep, the characters which represent the irrep, and the number of
    translations and rotations in the irrep.  The order of the point group
    is also provided (this is equal to the number of characters in an
    irrep).  */
class IrreducibleRepresentation {
  friend class CharacterTable;

  private:
    int g;         // the order of the group
    int degen;     // the degeneracy of the irrep
    int nrot_;     // the number of rotations in this irrep
    int ntrans_;   // the number of translations in this irrep
    int complex_;  // true if this irrep has a complex representation
    char *symb;    // mulliken symbol for this irrep
    char *csymb;    // mulliken symbol for this irrep w/o special characters

    SymRep *rep;   // representation matrices for the symops

  public:
    IrreducibleRepresentation();
    IrreducibleRepresentation(const IrreducibleRepresentation&);
    /** This constructor takes as arguments the order of the point group,
     the degeneracy of the irrep, and the Mulliken symbol of the irrep.
     The Mulliken symbol is copied internally. */
    IrreducibleRepresentation(int,int,const char*,const char* =0);

    ~IrreducibleRepresentation();

    IrreducibleRepresentation& operator=(const IrreducibleRepresentation&);

    /// Initialize the order, degeneracy, and Mulliken symbol of the irrep.
    void init(int =0, int =0, const char* =0, const char* =0);

    /// Returns the order of the group.
    int order() const { return g; }

    /// Returns the degeneracy of the irrep.
    int degeneracy() const { return degen; }

    /// Returns the value of complex_.
    int complex() const { return complex_; }

    /// Returns the number of projection operators for the irrep.
    int nproj() const { return degen*degen; }

    /// Returns the number of rotations associated with the irrep.
    int nrot() const { return nrot_; }

    /// Returns the number of translations associated with the irrep.
    int ntrans() const { return ntrans_; }

    /// Returns the Mulliken symbol for the irrep.
    const char * symbol() const { return symb; }

    /** Returns the Mulliken symbol for the irrep without special
        characters.
    */
    const char * symbol_ns() const { return (csymb?csymb:symb); }

    /** Returns the character for the i'th symmetry operation of the point
     group. */
    double character(int i) const {
      return complex_ ? 0.5*rep[i].trace() : rep[i].trace();
    }

    /// Returns the element (x1,x2) of the i'th representation matrix.
    double p(int x1, int x2, int i) const { return rep[i](x1,x2); }

    /** Returns the character for the d'th contribution to the i'th
     representation matrix. */
    double p(int d, int i) const {
      int dc=d/degen; int dr=d%degen;
      return rep[i](dr,dc);
    }

    /** This prints the irrep to the given file, or stdout if none is
     given.  The second argument is an optional string of spaces to offset
     by. */
     //void print(FILE *out=outfile) const;
};

} // close namespace panache

#endif //PANACHE_IRREDUCIBLEREPRESENTATION_H
