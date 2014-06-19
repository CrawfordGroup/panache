/*! \file
 *  \brief Simple class to hold molecule information 
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_MOLECULE_H
#define PANACHE_MOLECULE_H

#include <vector>
#include <string>
#include <memory>

#include "Vector3.h"

namespace panache {


/*!
 * \brief Class to hold information about the system
 */
class Molecule
{
private:

    /*!
     * \brief Holds information about an atom or other center
     */ 
    struct Atom
    {
        std::string symbol; //!< Atomic symbol
        Vector3 coord; //!< coordinates

        Atom(const std::string & symbol, double x, double y, double z)
            : symbol(symbol),coord({x,y,z})
        {   }

    };
    
    std::vector<Atom> atoms_;

public:
    // delete other copy constructors, etc
    Molecule(const Molecule& other) = delete;
    Molecule & operator=(const Molecule& other) = delete;

    /*!
     * \brief Default constructor
     */
    Molecule()
    { };


    // no specific destructor needed


    /*!
     * \brief Add an atom to the molecule
     * \param [in] x cartesian coordinate
     * \param [in] y cartesian coordinate
     * \param [in] z cartesian coordinate
     * \param [in] symb atomic symbol to use
     */
    void add_atom(double x, double y, double z, const std::string & symb)
    {
        atoms_.push_back(Atom(symb, x, y, z));
    }
                  

    /*!
     *  \brief Returns the number of atoms
     */
    unsigned int natom() const { return atoms_.size(); }


    /*!
     *  \brief Returns the x-coordinate of an atom
     *  \param [in] atom Number of the atom to get the x-coordinate for
     *  \return The x-coordinate of atom \p atom
     */
    double x(int atom) const
    {
        return atoms_[atom].coord[0];
    }


    /*!
     *  \brief Returns the y-coordinate of an atom
     *  \param [in] atom Number of the atom to get the y-coordinate for
     *  \return The y-coordinate of atom \p atom
     */
    double y(int atom) const
    {
        return atoms_[atom].coord[1];
    }


    /*!
     *  \brief Returns the z-coordinate of an atom
     *  \param [in] atom Number of the atom to get the z-coordinate for
     *  \return The z-coordinate of atom \p atom
     */
    double z(int atom) const
    {
        return atoms_[atom].coord[2];
    }


    /*!
     * \brief Returns a Vector3 with x, y, z position of an atom
     * \param [in] atom Number of the atom to get the coordinates for
     * \return All coordinates of atom \p atom
     */
    Vector3 xyz(int atom) const
    {
        return atoms_[atom].coord;
    }


    /*!
     *  \brief Returns the x, y, or z-coordinate of an atom
     *  \param [in] atom Number of the atom to get the coordinate for
     *  \param [in] _xyz Coordinate number to get (0 = x, 1 = y, 2 = z)
     *  \return The given coordinate of atom \p atom
     */
    double xyz(int atom, int _xyz) const
    {
        return atoms_[atom].coord[_xyz];
    }



    /*!
     * \brief Returns the atomic symbol for an atom
     * \param [in] atom Number of the atom to get the symbol for
     * \return Symbol of atom \p atom
     */
    std::string symbol(int atom) const
    {
        return atoms_[atom].symbol;
    }
};


/*!
 * \brief A shared pointer to a Molecule object
 */
typedef std::shared_ptr<Molecule> SharedMolecule;

}

#endif //PANACHE_MOLECULE_H
