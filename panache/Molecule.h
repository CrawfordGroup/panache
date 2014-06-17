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

#ifndef PANACHE_MOLECULE_H
#define PANACHE_MOLECULE_H

#include <vector>
#include <string>
#include <memory>

#include "Vector3.h"

namespace panache {

class Molecule
{
private:
    struct Atom
    {
        std::string symbol;
        double x,y,z;

        Atom(const std::string & symbol, double x, double y, double z)
            : symbol(symbol),x(x),y(y),z(z) 
        {   }

    };
    
    std::vector<Atom> atoms_;

public:
    Molecule();

    // delete other copy constructors, etc
    Molecule(const Molecule& other) = delete;
    Molecule & operator=(const Molecule& other) = delete;
    ~Molecule();

    /**
     * Add an atom to the molecule
     * \param x cartesian coordinate
     * \param y cartesian coordinate
     * \param z cartesian coordinate
     * \param symb atomic symbol to use
     */
    void add_atom(double x, double y, double z, const std::string & symb);
                  

    /// Number of atoms
    unsigned int natom() const { return atoms_.size(); }

    /// x position of atom
    double x(int atom) const;
    /// y position of atom
    double y(int atom) const;
    /// z position of atom
    double z(int atom) const;

    /// Returns a Vector3 with x, y, z position of atom
    Vector3 xyz(int atom) const;

    /// Returns x, y, or z component of 'atom'
    double xyz(int atom, int _xyz) const;

    /// Returns the cleaned up label of the atom (C2 => C, H4 = H)
    std::string symbol(int atom) const;
};

typedef std::shared_ptr<Molecule> SharedMolecule;

}

#endif //PANACHE_MOLECULE_H
