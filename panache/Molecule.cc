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

#include <algorithm>
#include <sstream>

#include "Molecule.h"
#include "Output.h"

using namespace std;

namespace panache {

Molecule::Molecule()
{
}

Molecule::~Molecule()
{
}

void Molecule::add_atom(double x, double y, double z,
                        const std::string & label)
{
    atoms_.push_back(Atom(label, x, y, z));
}

std::string Molecule::symbol(int atom) const
{
    return atoms_[atom].symbol;
}

Vector3 Molecule::xyz(int atom) const
{
    return Vector3({atoms_[atom].x, atoms_[atom].y, atoms_[atom].z});
}

double Molecule::x(int atom) const
{
    return atoms_[atom].x;
}

double Molecule::y(int atom) const
{
    return atoms_[atom].y;
}

double Molecule::z(int atom) const
{
    return atoms_[atom].z;
}

} // end namespace panache


