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
#include "Vector.h"
#include "Output.h"
#include "Masses.h"
#include "PhysConst.h"
#include "Element_to_Z.h"

using namespace std;


// the third parameter of from_string() should be
// one of std::hex, std::dec or std::oct
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

// used by 'if_to_invert_axis' and 'inertia_tensor'
#define ZERO 1.0E-14

namespace panache
{

/**
 * Interprets a string as an integer, throwing if it's unsuccesful.
 */
static
int
str_to_int(const std::string& s)
{
    int i;
    std::istringstream iss(s);
    if((iss >> std::dec >> i).fail())
        throw RuntimeError("Unable to convert " + s + " to an integer");
    return i;
}

/**
 * Interprets a string as an double, throwing if it's unsuccesful.
 */
static
double
str_to_double(const std::string& s)
{
    double d;
    std::istringstream iss(s);
    if((iss >> std::dec >> d).fail())
        throw RuntimeError("Unable to convert " + s + " to a double");
    return d;
}

static string string_to_upper_copy(string str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

static string string_to_lower_copy(string str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}


void if_to_invert_axis(const Vector3& v1, int& must_invert, int& should_invert, double& maxproj)
{
    int xyz, nzero;
    double vabs;

    maxproj = 0.0;
    must_invert = 0;
    should_invert = 0;

    nzero = 0;

    for(xyz=0; xyz<3; xyz++)
    {

        vabs = fabs(v1[xyz]);

        if (vabs < ZERO)
            nzero++;

        if (vabs > fabs(maxproj))
        {
            maxproj = v1[xyz];
        }

    }

    if (nzero == 2)
    {
        if (maxproj < 0.0)
            must_invert = 1;
        else
            must_invert = 0;
    }
    else if (nzero < 2)
    {
        if (maxproj < 0.0)
            should_invert = 1;
        else
            should_invert = 0;
    }
}

Molecule::Molecule():
    fix_orientation_(false),
    move_to_com_(true),
    molecular_charge_(0),
    multiplicity_(1),

    units_(Bohr),
    input_units_to_au_(1.0),
//    units_(Angstrom),
//    input_units_to_au_(1.0/pc_bohr2angstroms),

    nunique_(0),
    nequiv_(0),
    equiv_(0),
    multiplicity_specified_(false),
    charge_specified_(false),
    atom_to_unique_(0),
    //old_symmetry_frame_(0)
    reinterpret_coordentries_(true),
    lock_frame_(false),
    full_pg_(PG_C1),
    full_pg_n_(1)
{
}

Molecule::~Molecule()
{
    clear();
    release_symmetry_information();
    //if (old_symmetry_frame_)
    //  delete old_symmetry_frame_;
}

Molecule& Molecule::operator=(const Molecule& other)
{
    // Self assignment is bad
    if (this == &other)
        return *this;

    all_variables_           = other.all_variables_;
    fragments_               = other.fragments_;
    fragment_charges_        = other.fragment_charges_;
    fragment_multiplicities_ = other.fragment_multiplicities_;
    fix_orientation_         = other.fix_orientation_;
    move_to_com_             = other.move_to_com_;
    molecular_charge_        = other.molecular_charge_;
    multiplicity_            = other.multiplicity_;
    units_                   = other.units_;
    input_units_to_au_       = other.input_units_to_au_;
    all_variables_           = other.all_variables_;
    fragment_types_          = other.fragment_types_;
    geometry_variables_      = other.geometry_variables_;
    charge_specified_        = other.charge_specified_;
    multiplicity_specified_  = other.multiplicity_specified_;
    reinterpret_coordentries_= other.reinterpret_coordentries_;
    lock_frame_              = other.lock_frame_;
    zmat_                    = other.zmat_;

    // These are symmetry related variables, and are filled in by the following functions
    pg_             = std::shared_ptr<PointGroup>();
    nunique_        = 0;
    nequiv_         = 0;
    equiv_          = 0;
    atom_to_unique_ = 0;
    symmetry_from_input_ = other.symmetry_from_input_;
    form_symmetry_information();
    full_pg_        = other.full_pg_;
    full_pg_n_      = other.full_pg_n_;

    atoms_.clear();
    // Deep copy the map of variables
    std::vector<std::shared_ptr<CoordEntry> >::const_iterator iter = other.full_atoms_.begin();
    for(; iter != other.full_atoms_.end(); ++iter)
        full_atoms_.push_back((*iter)->clone(full_atoms_, geometry_variables_));

    // This is called here, so that the atoms list is populated
    update_geometry();

    return *this;
}

Molecule::Molecule(const Molecule& other)
{
    *this = other;
}

void Molecule::set_reinterpret_coordentry(bool rc)
{
    reinterpret_coordentries_ = rc;
}

/// Addition
//Molecule Molecule::operator+(const Molecule& other)
//{

//}

///// Subtraction
//Molecule Molecule::operator-(const Molecule& other)
//{

//}

/// Plus equals
void Molecule::operator+=(const Molecule& other)
{
    throw RuntimeError("Empty method?");
}

void Molecule::clear()
{
    lock_frame_ = false;
    atoms_.empty();
    full_atoms_.empty();
}

void Molecule::add_atom(int Z, double x, double y, double z,
                        const char *label, double mass, double charge, int lineno)
{
    lock_frame_ = false;
    Vector3 temp(x, y, z);
    std::string l(label);

    if (atom_at_position2(temp) == -1)
    {
        // Dummies go to full_atoms_, ghosts need to go to both.
        full_atoms_.push_back(std::shared_ptr<CoordEntry>(new CartesianEntry(full_atoms_.size(), Z, charge, mass, l, l,
                              std::shared_ptr<CoordValue>(new NumberValue(x)),
                              std::shared_ptr<CoordValue>(new NumberValue(y)),
                              std::shared_ptr<CoordValue>(new NumberValue(z)))));
        if(strcmp(label, "X") && strcmp(label, "x")) atoms_.push_back(full_atoms_.back());
    }
    else
    {
        throw RuntimeError("Molecule::add_atom: Adding atom on top of an existing atom.");
    }
}

double Molecule::mass(int atom) const
{
    if (atoms_[atom]->mass() != 0.0)
        return atoms_[atom]->mass();

//    if (fabs(atoms_[atom]->Z() - static_cast<int>(atoms_[atom]->Z())) > 0.0)
//        output::printf("WARNING: Obtaining masses from atom with fractional charge...may be incorrect!!!\n");

    return an2masses[static_cast<int>(atoms_[atom]->Z())];
}

std::string Molecule::symbol(int atom) const
{
    return atoms_[atom]->symbol();
}

std::string Molecule::fsymbol(int atom) const
{
    return full_atoms_[atom]->symbol();
}
std::string Molecule::label(int atom) const
{
    return atoms_[atom]->label();
}

int Molecule::atom_at_position1(double *coord, double tol) const
{
    Vector3 b(coord);
    for (int i=0; i < natom(); ++i)
    {
        Vector3 a = xyz(i);
        if (b.distance(a) < tol)
            return i;
    }
    return -1;
}

int Molecule::atom_at_position2(Vector3& b, double tol) const
{
    for (int i=0; i < natom(); ++i)
    {
        Vector3 a = xyz(i);
        if (b.distance(a) < tol)
            return i;
    }
    return -1;
}

Vector3 Molecule::center_of_mass() const
{
    Vector3 ret;
    double total_m;

    ret = 0.0;
    total_m = 0.0;

    for (int i=0; i<natom(); ++i)
    {
        double m = mass(i);
        ret += m * xyz(i);
        total_m += m;
    }

    ret *= 1.0/total_m;

    return ret;
}

Matrix Molecule::distance_matrix() const
{
    Matrix distance("Distances between atoms in Bohr", natom(), natom());

    for (int i=0; i<natom(); ++i)
    {
        for (int j=0; j<=i; ++j)
        {
            distance(i, j) = distance(j, i) = xyz(i).distance(xyz(j));
        }
    }

    return distance;
}

void Molecule::translate(const Vector3& r)
{
    Vector3 temp;
    for (int i=0; i<nallatom(); ++i)
    {
        temp = input_units_to_au_ * full_atoms_[i]->compute();
        temp += r;
        temp = temp/input_units_to_au_;
        full_atoms_[i]->set_coordinates(temp[0], temp[1], temp[2]);
    }
}

void Molecule::move_to_com()
{
    Vector3 com = -center_of_mass();
    translate(com);
}

Matrix Molecule::geometry() const
{
    Matrix geom(natom(), 3);
    for (int i=0; i<natom(); ++i)
    {
        geom(i, 0) = x(i);
        geom(i, 1) = y(i);
        geom(i, 2) = z(i);
    }

    return geom;
}

Matrix Molecule::full_geometry() const
{
    Matrix geom(nallatom(), 3);
    for (int i=0; i<nallatom(); ++i)
    {
        geom(i, 0) = fx(i);
        geom(i, 1) = fy(i);
        geom(i, 2) = fz(i);
    }

    return geom;
}

void Molecule::set_geometry(double** geom)
{
    lock_frame_ = false;
    bool dummy_found = false;
    for (int i=0; i<nallatom(); ++i)
    {
        if(full_atoms_[i]->symbol() == "X")
        {
            dummy_found = true;
            break;
        }
    }
    // We don't track the coordinates of the dummy atoms.  For now, just convert the entries
    // to Cartesians if the entry contains
    if(dummy_found)
    {
        atoms_.clear();
        int count = 0;
        std::vector<int> fragment_changes;
        for(int i = 0; i < fragments_.size(); ++i)
            fragment_changes.push_back(0);
        for (int i=0; i<nallatom(); ++i)
        {
            std::shared_ptr<CoordEntry> at = full_atoms_[i];

            if(at->symbol() == "X")
            {
                // Find out which fragment this atom is removed from, then bail
                bool found = false;
                for(int frag = 0; frag < fragments_.size(); ++frag)
                {
                    if(i >= fragments_[frag].first && i < fragments_[frag].second)
                    {
                        found = true;
                        fragment_changes[frag]++;
                        break;
                    }
                }
                if(!found)
                    throw RuntimeError("Problem converting ZMatrix coordinates to Cartesians."
                                       "Try again, without dummy atoms.");
                continue;
            }

            int entrynum = at->entry_number();
            double zval = at->Z();
            double charge = at->charge();
            double mass = at->mass();
            std::string symbol = at->symbol();
            std::string label = at->label();
            std::shared_ptr<CoordEntry> new_atom(
                new CartesianEntry(entrynum,
                                   zval,
                                   charge,
                                   mass,
                                   symbol,
                                   label,
                                   std::shared_ptr<CoordValue>(new NumberValue(geom[count][0]/input_units_to_au_)),
                                   std::shared_ptr<CoordValue>(new NumberValue(geom[count][1]/input_units_to_au_)),
                                   std::shared_ptr<CoordValue>(new NumberValue(geom[count][2]/input_units_to_au_))
                                  ));
            // Copy over all known basis sets
            const std::map<std::string, std::string>& basissets = at->basissets();
            std::map<std::string, std::string>::const_iterator bs = basissets.begin();
            for(; bs != basissets.end(); ++bs)
                new_atom->set_basisset(bs->second, bs->first);
            atoms_.push_back(new_atom);
            count++;
        }
        full_atoms_.clear();
        for(int i = 0; i < atoms_.size(); ++i)
            full_atoms_.push_back(atoms_[i]);
        // Now change the bounds of each fragment, to reflect the missing dummy atoms
        int cumulative_count = 0;
        for(int frag = 0; frag < fragments_.size(); ++frag)
        {
            fragments_[frag].first -= cumulative_count;
            cumulative_count += fragment_changes[frag];
            fragments_[frag].second -= cumulative_count;
        }
        geometry_variables_.clear();
    }
    else
    {
        for (int i=0; i<natom(); ++i)
        {
            atoms_[i]->set_coordinates(geom[i][0] / input_units_to_au_,
                                       geom[i][1] / input_units_to_au_,
                                       geom[i][2] / input_units_to_au_);
        }
    }
}

void Molecule::set_full_geometry(double** geom)
{
    lock_frame_ = false;
    for (int i=0; i<nallatom(); ++i)
    {
        full_atoms_[i]->set_coordinates(geom[i][0] / input_units_to_au_,
                                        geom[i][1] / input_units_to_au_,
                                        geom[i][2] / input_units_to_au_);
    }
}

void Molecule::set_geometry(const Matrix& geom)
{
    lock_frame_ = false;
    set_geometry(geom.pointer());
}

void Molecule::set_full_geometry(const Matrix& geom)
{
    lock_frame_ = false;
    set_full_geometry(geom.pointer());
}

void Molecule::rotate(const Matrix& R)
{
    Matrix new_geom(natom(), 3);
    Matrix geom = geometry();

    // Multiple the geometry by the rotation matrix.
    new_geom.gemm(false, false, 1.0, geom, R, 0.0);

    set_geometry(new_geom);
}

void Molecule::rotate_full(const Matrix& R)
{
    Matrix new_geom(nallatom(), 3);
    Matrix geom = full_geometry();

    // Multiply the geometry by the rotation matrix.
    new_geom.gemm(false, false, 1.0, geom, R, 0.0);

    set_full_geometry(new_geom);
}

/**
 * Checks whether the user has specified the charge in the options, and returns the appropriate value.
 * @return The charge from the options keywords, if specified.  If not, the value passed to the molecule
 *         specification, which takes the default value provided by liboptions if not specified.
 */
int Molecule::molecular_charge() const
{
    return molecular_charge_;
}

/**
 * Checks whether the user has specified the multiplicity in the options, and returns the appropriate value.
 * @return The multiplicity from the options keywords, if specified.  If not, the value passed to the molecule
 *         specification, which takes the default value provided by liboptions if not specified.
 */
int Molecule::multiplicity() const
{
    return multiplicity_;
}

void Molecule::symmetrize_to_abelian_group(double tol)
{
    // The coordinates are a bit crude, so we symmetrize them
    // First, populate the atom list
    reinterpret_coordentries();
    // Now, redetect the symmetry with a really crude tolerance
    SharedMatrix frame = symmetry_frame(tol);
    // Put it on the center of mass and rotate
    move_to_com();
    rotate_full(*frame.get());
    set_point_group(find_point_group(tol));
    // Clean up the molecule, to make sure it actually has the correct symmetry
    symmetrize();
}

void Molecule::reinterpret_coordentries()
{
    atoms_.clear();
    EntryVectorIter iter;
    for (iter = full_atoms_.begin(); iter != full_atoms_.end(); ++iter)
    {
        (*iter)->invalidate();
    }
    int temp_charge = molecular_charge_;
    int temp_multiplicity = multiplicity_;
    molecular_charge_ = 0;
    multiplicity_    = 1;
    for(int fragment = 0; fragment < fragments_.size(); ++fragment)
    {
        if(fragment_types_[fragment] == Absent)
            continue;
        if(fragment_types_[fragment] == Real)
        {
            molecular_charge_ += fragment_charges_[fragment];
            multiplicity_    += fragment_multiplicities_[fragment] - 1;
        }
        for(int atom = fragments_[fragment].first; atom < fragments_[fragment].second; ++atom)
        {
            full_atoms_[atom]->compute();
            full_atoms_[atom]->set_ghosted(fragment_types_[fragment] == Ghost);
            if(full_atoms_[atom]->symbol() != "X") atoms_.push_back(full_atoms_[atom]);
        }
    }
    // TODO: This is a hack to ensure that set_multiplicity and set_molecular_charge
    // work for single-fragment molecules.
    if (fragments_.size() < 2)
    {
        molecular_charge_ = temp_charge;
        multiplicity_ = temp_multiplicity;
    }

    if(zmat_)
    {
        // Even if the user asked us to lock the frame, we should reorient here for zmatrices
        SharedMatrix frame = symmetry_frame();
        rotate_full(*frame.get());
        move_to_com();
    }

}

void Molecule::update_geometry()
{
    if (fragments_.size() == 0)
        throw RuntimeError("Molecule::update_geometry: There are no fragments in this molecule.");

    // Idempotence condition
    if (lock_frame_)
        return;


    if (reinterpret_coordentries_)
        reinterpret_coordentries();

    if (move_to_com_)
        move_to_com();

    // If the no_reorient command was given, don't reorient
    if (fix_orientation_ == false)
    {
        // Now we need to rotate the geometry to its symmetry frame
        // to align the axes correctly for the point group
        // symmetry_frame looks for the highest point group so that we can align
        // the molecule according to its actual symmetry, rather than the symmetry
        // the the user might have provided.
        SharedMatrix frame = symmetry_frame();
        rotate_full(*frame.get());
    }

    // Recompute point group of the molecule, so the symmetry info is updated to the new frame
    set_point_group(find_point_group());
    set_full_point_group();

    symmetrize(); // Symmetrize the molecule to remove any noise.

    lock_frame_ = true;
}

void Molecule::activate_all_fragments()
{
    lock_frame_ = false;
    for(int i = 0; i < fragment_types_.size(); ++i)
    {
        fragment_types_[i] = Real;
    }
}

int Molecule::nactive_fragments()
{
    int n = 0;
    for(int i = 0; i < fragment_types_.size(); ++i)
    {
        if ( fragment_types_[i] == Real ) n++;
    }
    return n;
}

void Molecule::deactivate_all_fragments()
{
    lock_frame_ = false;
    for(int i = 0; i < fragment_types_.size(); ++i)
    {
        fragment_types_[i] = Absent;
    }
}

void Molecule::set_active_fragment(int fragment)
{
    lock_frame_ = false;
    fragment_types_[fragment - 1] = Real;
}

void Molecule::set_ghost_fragment(int fragment)
{
    lock_frame_ = false;
    fragment_types_[fragment - 1] = Ghost;
}

std::shared_ptr<Molecule> Molecule::extract_subsets(const std::vector<int> &real_list, const std::vector<int> &ghost_list) const
{
    if(ghost_list.size() + real_list.size() > fragments_.size())
        throw RuntimeError("The sum of real- and ghost-atom subsets is greater than the number of subsets");

    std::shared_ptr<Molecule> clone(new Molecule(*this));
    clone->deactivate_all_fragments();
    for(int fragment = 0; fragment < real_list.size(); ++fragment)
    {
        clone->set_active_fragment(real_list[fragment]+1); // The active fragment code subtracts 1
    }
    for(int fragment = 0; fragment < ghost_list.size(); ++fragment)
    {
        clone->set_ghost_fragment(ghost_list[fragment]+1); // The ghost fragment code subtracts 1
    }
    clone->update_geometry();
    return clone;
}

Matrix* Molecule::inertia_tensor() const
{
    int i;
    Matrix* tensor = new Matrix("Inertia Tensor", 3, 3);
    Matrix& temp = *tensor;

    for (i = 0; i < natom(); i++)
    {
        // I(alpha, alpha)
        temp(0, 0) += mass(i) * (y(i) * y(i) + z(i) * z(i));
        temp(1, 1) += mass(i) * (x(i) * x(i) + z(i) * z(i));
        temp(2, 2) += mass(i) * (x(i) * x(i) + y(i) * y(i));

        // I(alpha, beta)
        temp(0, 1) -= mass(i) * x(i) * y(i);
        temp(0, 2) -= mass(i) * x(i) * z(i);
        temp(1, 2) -= mass(i) * y(i) * z(i);
    }

    //    mirror
    temp(1, 0) = temp(0, 1);
    temp(2, 0) = temp(0, 2);
    temp(2, 1) = temp(1, 2);

    // Check the elements for zero and make them a hard zero.
    for (int i=0; i < 3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            if (fabs(tensor->get(i, j)) < ZERO)
                tensor->set(i, j, 0.0);
        }
    }

    return tensor;
}

Vector Molecule::rotational_constants(double zero_tol) const
{

    SharedMatrix pI(inertia_tensor());
    Vector evals(3);
    SharedMatrix eigenvectors(new Matrix(3, 3));
    pI->diagonalize(eigenvectors, evals, ascending);

    // Conversion factor from moments to rotational constants.
    double im2rotconst = pc_h / (8 * pc_pi * pc_pi * pc_c);
    // Add factor to put moments into SI units - give result in wavenumbers.
    im2rotconst /= (pc_bohr2m * pc_bohr2m * pc_amu2kg * 100);

    Vector rot_const(3);
    for (int i=0; i<3; ++i)
    {
        if (evals[i] < zero_tol)
            rot_const[i] = 0.0;
        else
            rot_const[i] = im2rotconst/evals[i];
    }

    /*
      output::printf("\n\tRotational constants (cm^-1) :\n");
      if (rot_const[0] == 0) // linear
        output::printf("\tA = **********  ");
      else               // non-linear
        output::printf("\tA = %10.5lf  ", rot_const[0]);

      if (rot_const[1] == 0) // atom
        output::printf("  B = **********    C = **********  \n");
      else               // molecule
        output::printf("  B = %10.5lf   C = %10.5lf\n", rot_const[1], rot_const[2]);
    */
    return rot_const;
}

RotorType Molecule::rotor_type(double zero_tol) const
{

    Vector rot_const = rotational_constants();

    // Determine degeneracy of rotational constants.
    double tmp, abs, rel;
    int degen = 0;
    for (int i=0; i<2; i++)
    {
        for (int j=i+1; j<3 && degen<2; j++)
        {
            abs = fabs(rot_const[i] - rot_const[j]);
            tmp = (rot_const[i] > rot_const[j]) ? rot_const[i] : rot_const[j];
            if (abs > 1.0E-14)
                rel = abs/tmp;
            else
                rel = 0.0;
            if (rel < zero_tol)
                degen++;
        }
    }
    //output::printf("\tDegeneracy is %d\n", degen);

    // Determine rotor type
    RotorType rotor_type;

    if (natom() == 1)
        rotor_type = RT_ATOM;
    else if (rot_const[0] == 0.0)  // A == 0, B == C
        rotor_type = RT_LINEAR;
    else if (degen == 2)           // A == B == C
        rotor_type = RT_SPHERICAL_TOP;
    else if (degen == 1)           // A  > B == C
        rotor_type = RT_SYMMETRIC_TOP;  // A == B > C
    else
        rotor_type = RT_ASYMMETRIC_TOP; // A != B != C

    return rotor_type;
}

//
// Symmetry
//
bool Molecule::has_inversion(Vector3& origin, double tol) const
{
    for (int i=0; i<natom(); ++i)
    {
        Vector3 inverted = origin-(xyz(i) - origin);
        int atom = atom_at_position2(inverted, tol);
        if (atom < 0 || !atoms_[atom]->is_equivalent_to(atoms_[i]))
        {
            return false;
        }
    }
    return true;
}

bool Molecule::is_plane(Vector3& origin, Vector3& uperp, double tol) const
{
    for (int i=0; i<natom(); ++i)
    {
        Vector3 A = xyz(i)-origin;
        Vector3 Apar = uperp.dot(A)*uperp;
        Vector3 Aperp = A - Apar;
        A = (Aperp- Apar) + origin;
        int atom = atom_at_position2(A, tol);
        if (atom < 0 || !atoms_[atom]->is_equivalent_to(atoms_[i]))
        {
            return false;
        }
    }
    return true;
}

bool Molecule::is_axis(Vector3& origin, Vector3& axis, int order, double tol) const
{
    for (int i=0; i<natom(); ++i)
    {
        Vector3 A = xyz(i) - origin;
        for (int j=1; j<order; ++j)
        {
            Vector3 R = A;
            R.rotate(j*2.0*M_PI/order, axis);
            R += origin;
            int atom = atom_at_position2(R, tol);
            if (atom < 0 || !atoms_[atom]->is_equivalent_to(atoms_[i]))
            {
                return false;
            }
        }
    }
    return true;
}

enum AxisName { XAxis, YAxis, ZAxis };

static AxisName like_world_axis(Vector3& axis, const Vector3& worldxaxis, const Vector3& worldyaxis, const Vector3& worldzaxis)
{
    AxisName like;
    double xlikeness = fabs(axis.dot(worldxaxis));
    double ylikeness = fabs(axis.dot(worldyaxis));
    double zlikeness = fabs(axis.dot(worldzaxis));
    if ((xlikeness - ylikeness) > 1.0e-12 && (xlikeness - zlikeness) > 1.0e-12)
    {
        like = XAxis;
        if (axis.dot(worldxaxis) < 0) axis = - axis;
    }
    else if ((ylikeness - zlikeness) > 1.0e-12)
    {
        like = YAxis;
        if (axis.dot(worldyaxis) < 0) axis = - axis;
    }
    else
    {
        like = ZAxis;
        if (axis.dot(worldzaxis) < 0) axis = - axis;
    }
    return like;
}

void Molecule::is_linear_planar(bool& linear, bool& planar, double tol) const
{
    if (natom() < 3)
    {
        linear = true;
        planar = true;
        return;
    }

    // find three atoms not on the same line
    Vector3 A = xyz(0);
    Vector3 B = xyz(1);
    Vector3 BA = B-A;
    BA.normalize();
    Vector3 CA;

    int i;
    double min_BAdotCA = 1.0;
    for (i=2; i<natom(); ++i)
    {
        Vector3 tmp = xyz(i) - A;
        tmp.normalize();
        if (fabs(BA.dot(tmp)) < min_BAdotCA)
        {
            CA = tmp;
            min_BAdotCA = fabs(BA.dot(tmp));
        }
    }
    if (min_BAdotCA >= 1.0 - tol)
    {
        linear = true;
        planar = true;
        return;
    }

    linear = false;
    if (natom() < 4)
    {
        planar = true;
        return;
    }

    // check for nontrivial planar molecules
    Vector3 BAxCA = BA.cross(CA);
    BAxCA.normalize();
    for (i=2; i<natom(); ++i)
    {
        Vector3 tmp = xyz(i)-A;
        if (fabs(tmp.dot(BAxCA)) > tol)
        {
            planar = false;
            return;
        }
    }
    planar = true;
}

int Molecule::atom_to_unique_offset(int iatom) const
{
    int iuniq = atom_to_unique_[iatom];
    int nequiv = nequiv_[iuniq];
    for (int i=0; i<nequiv; ++i)
    {
        if (equiv_[iuniq][i] == iatom)
            return i;
    }
    throw RuntimeError("Molecule::atom_to_unique_offset: I should've found the atom requested...but didn't.");
    return -1;
}

int Molecule::max_nequivalent() const
{
    int max = 0;
    for (int i=0; i<nunique(); ++i)
        if (max < nequivalent(i))
            max = nequivalent(i);
    return max;
}

std::shared_ptr<Matrix> Molecule::symmetry_frame(double tol)
{
    int i, j;

    Vector3 com = center_of_mass();

    Vector3 worldxaxis(1.0, 0.0, 0.0);
    Vector3 worldyaxis(0.0, 1.0, 0.0);
    Vector3 worldzaxis(0.0, 0.0, 1.0);

    bool linear, planar;
    is_linear_planar(linear, planar, tol);

    bool have_inversion = has_inversion(com, tol);

    // check for C2 axis
    Vector3 c2axis;
    bool have_c2axis = false;
    if (natom() < 2)
    {
        have_c2axis = true;
        c2axis = Vector3(0.0, 0.0, 1.0);
    }
    else if (linear)
    {
        have_c2axis = true;
        c2axis = xyz(1) - xyz(0);
        c2axis.normalize();
    }
    else if (planar && have_inversion)
    {
        // there is a c2 axis that won't be found using the usual
        // algorithm. fine two noncolinear atom-atom vectors (we know
        // that linear == 0)
        Vector3 BA = xyz(1) - xyz(0);
        BA.normalize();
        for (i=2; i<natom(); ++i)
        {
            Vector3 CA = xyz(i) - xyz(0);
            CA.normalize();
            Vector3 BAxCA = BA.cross(CA);
            if (BAxCA.norm() > tol)
            {
                have_c2axis = true;
                BAxCA.normalize();
                c2axis = BAxCA;
                break;
            }
        }
    }
    else
    {
        // loop through pairs of atoms o find c2 axis candidates
        for (i=0; i<natom(); ++i)
        {
            Vector3 A = xyz(i) - com;
            double AdotA = A.dot(A);
            for (j=0; j<=i; ++j)
            {
                // the atoms must be identical
                if (!atoms_[i]->is_equivalent_to(atoms_[j])) continue;
                Vector3 B = xyz(j)-com;
                // the atoms must be the same distance from the com
                if (fabs(AdotA - B.dot(B)) > tol) continue;
                Vector3 axis = A+B;
                // atoms colinear with the com don't work
                if (axis.norm() < tol) continue;
                axis.normalize();
                if (is_axis(com, axis, 2, tol))
                {
                    have_c2axis = true;
                    c2axis = axis;
                    goto symmframe_found_c2axis;
                }
            }
        }
    }
symmframe_found_c2axis:

    AxisName c2like = ZAxis;
    if (have_c2axis)
    {
        // try to make the sign of the axis correspond to one of the
        // world axes
        c2like = like_world_axis(c2axis, worldxaxis, worldyaxis, worldzaxis);
    }

    // check for c2 axis perp to first c2 axis
    Vector3 c2axisperp;
    bool have_c2axisperp = false;
    if (have_c2axis)
    {
        if (natom() < 2)
        {
            have_c2axisperp = true;
            c2axisperp = Vector3(1.0, 0.0, 0.0);
        }
        else if (linear)
        {
            if (have_inversion)
            {
                have_c2axisperp = true;
                c2axisperp = c2axis.perp_unit(Vector3(0.0,0.0,1.0));
            }
        }
        else
        {
            // loop through paris of atoms to find c2 axis candidates
            for (i=0; i<natom(); ++i)
            {
                Vector3 A = xyz(i) - com;
                double AdotA = A.dot(A);
                for (j=0; j<i; ++j)
                {
                    // the atoms must be identical
                    if (!atoms_[i]->is_equivalent_to(atoms_[j])) continue;
                    Vector3 B = xyz(j) - com;
                    // the atoms must be the same distance from the com
                    if (fabs(AdotA - B.dot(B)) > tol) continue;
                    Vector3 axis= A+B;
                    // atoms colinear with the com don't work
                    if (axis.norm() < tol) continue;
                    axis.normalize();
                    // if axis is not perp continue
                    if (fabs(axis.dot(c2axis)) > tol) continue;
                    if (is_axis(com, axis, 2, tol))
                    {
                        have_c2axisperp = true;
                        c2axisperp = axis;
                        goto symmframe_found_c2axisperp;
                    }
                }
            }
        }
    }
symmframe_found_c2axisperp:

    AxisName c2perplike;
    if (have_c2axisperp)
    {
        // try to make the sign of the axis correspond to one of
        // the world axes
        c2perplike = like_world_axis(c2axisperp, worldxaxis, worldyaxis, worldzaxis);

        // try to make c2axis the z axis
        if (c2perplike == ZAxis)
        {
            Vector3 tmpv = c2axisperp;
            tmpv = c2axisperp;
            c2axisperp = c2axis;
            c2axis = tmpv;
            c2perplike = c2like;
            c2like = ZAxis;
        }
        if (c2like != ZAxis)
        {
            if (c2like == XAxis) c2axis = c2axis.cross(c2axisperp);
            else c2axis = c2axisperp.cross(c2axis);
            c2like = like_world_axis(c2axis, worldxaxis, worldyaxis, worldzaxis);
        }
        // try to make c2axisperplike the x axis
        if (c2perplike == YAxis)
        {
            c2axisperp = c2axisperp.cross(c2axis);
            c2perplike = like_world_axis(c2axisperp, worldxaxis, worldyaxis, worldzaxis);
        }
    }

    // Check for vertical plane
    bool have_sigmav = false;
    Vector3 sigmav;
    if (have_c2axis)
    {
        if (natom() < 2)
        {
            have_sigmav = true;
            sigmav = c2axisperp;
        }
        else if (linear)
        {
            have_sigmav = true;
            if (have_c2axisperp)
            {
                sigmav = c2axisperp;
            }
            else
            {
                sigmav = c2axis.perp_unit(Vector3(0.0, 0.0, 1.0));
            }
        }
        else
        {
            // loop through pairs of atoms to find sigma v plane
            // candidates
            for (i=0; i<natom(); ++i)
            {
                Vector3 A = xyz(i) - com;
                double AdotA = A.dot(A);
                // the second atom can equal i because i might be
                // in the plane
                for (j=0; j<=i; ++j)
                {
                    // the atoms must be identical
                    if (!atoms_[i]->is_equivalent_to(atoms_[j])) continue;
                    Vector3 B = xyz(j) - com;
                    // the atoms must be the same distance from the com
                    if (fabs(AdotA - B.dot(B)) > tol) continue;
                    Vector3 inplane = B+A;
                    double norm_inplane = inplane.norm();
                    if (norm_inplane < tol) continue;
                    inplane *= 1.0/norm_inplane;
                    Vector3 perp = c2axis.cross(inplane);
                    double norm_perp = perp.norm();
                    if (norm_perp < tol) continue;
                    perp *= 1.0/norm_perp;
                    if (is_plane(com, perp, tol))
                    {
                        have_sigmav = true;
                        sigmav = perp;
                        goto symmframe_found_sigmav;
                    }
                }
            }
        }
    }

symmframe_found_sigmav:
    if (have_sigmav)
    {
        // try to make the sign of the oop vec correspond to one of
        // the world axes
        int sigmavlike = like_world_axis(sigmav, worldxaxis, worldyaxis, worldzaxis);

        // Choose sigmav to be the world x axis, if possible
        if (c2like == ZAxis && sigmavlike == YAxis)
        {
            sigmav = sigmav.cross(c2axis);
        }
        else if (c2like == YAxis && sigmavlike == ZAxis)
        {
            sigmav = c2axis.cross(sigmav);
        }
    }

    // under certain conditions i need to know if there is any sigma
    // plane
    bool have_sigma = false;
    Vector3 sigma;
    if (!have_inversion && !have_c2axis)
    {
        if (planar)
        {
            // find two noncolinear atom-atom vectors
            // we know that linear==0 since !have_c2axis
            Vector3 BA = xyz(1) - xyz(0);
            BA.normalize();
            for (i=2; i<natom(); ++i)
            {
                Vector3 CA = xyz(i) - xyz(0);
                CA.normalize();
                Vector3 BAxCA = BA.cross(CA);
                if (BAxCA.norm() > tol)
                {
                    have_sigma = true;
                    BAxCA.normalize();
                    sigma = BAxCA;
                    break;
                }
            }
        }
        else
        {
            // loop through pairs of atoms to contruct trial planes
            for (i=0; i<natom(); ++i)
            {
                Vector3 A = xyz(i) - com;
                double AdotA = A.dot(A);
                for (j=0; j<i; ++j)
                {
                    // the atomsmust be identical
                    if (!atoms_[i]->is_equivalent_to(atoms_[j])) continue;
                    Vector3 B = xyz(j)-com;
                    double BdotB = B.dot(B);
                    // the atoms must be the same distance from the com
                    if (fabs(AdotA - BdotB) > tol) continue;
                    Vector3 perp = B-A;
                    double norm_perp = perp.norm();
                    if (norm_perp < tol) continue;
                    perp *= 1.0 / norm_perp;
                    if (is_plane(com, perp, tol))
                    {
                        have_sigma = true;
                        sigma = perp;
                        goto found_sigma;
                    }
                }
            }
        }
    }
found_sigma:

    if (have_sigma)
    {
        // try to make the sign of the oop vec correspond to one of
        // the world axes
        double xlikeness = fabs(sigma.dot(worldxaxis));
        double ylikeness = fabs(sigma.dot(worldyaxis));
        double zlikeness = fabs(sigma.dot(worldzaxis));

        if (xlikeness > ylikeness && xlikeness > zlikeness)
        {
            if (sigma.dot(worldxaxis) < 0) sigma = -sigma;
        }
        else if (ylikeness > zlikeness)
        {
            if (sigma.dot(worldyaxis) < 0) sigma = -sigma;
        }
        else
        {
            if (sigma.dot(worldzaxis) < 0) sigma = -sigma;
        }
    }

    // Find the three axes for the symmetry frame
    Vector3 xaxis = worldxaxis;
    Vector3 yaxis;
    Vector3 zaxis = worldzaxis;
    if (have_c2axis)
    {
        zaxis = c2axis;
        if (have_sigmav)
        {
            xaxis = sigmav;
        }
        else if (have_c2axisperp)
        {
            xaxis = c2axisperp;
        }
        else
        {
            // any axis orthogonal to the zaxis will do
            xaxis = zaxis.perp_unit(zaxis);
        }
    }
    else if (have_sigma)
    {
        zaxis = sigma;
        xaxis = zaxis.perp_unit(zaxis);
    }

#define NOISY_ZERO 1.0e-8
    // Clean up our z axis
    if (fabs(zaxis[0]) < NOISY_ZERO)
        zaxis[0] = 0.0;
    if (fabs(zaxis[1]) < NOISY_ZERO)
        zaxis[1] = 0.0;
    if (fabs(zaxis[2]) < NOISY_ZERO)
        zaxis[2] = 0.0;

    // Clean up our x axis
    if (fabs(xaxis[0]) < NOISY_ZERO)
        xaxis[0] = 0.0;
    if (fabs(xaxis[1]) < NOISY_ZERO)
        xaxis[1] = 0.0;
    if (fabs(xaxis[2]) < NOISY_ZERO)
        xaxis[2] = 0.0;
#undef NOISY_ZERO

    // the y is then -x cross z
    yaxis = -xaxis.cross(zaxis);

//    output::printf("xaxis %20.14lf %20.14lf %20.14lf\n", xaxis[0], xaxis[1], xaxis[2]);
//    output::printf("yaxis %20.14lf %20.14lf %20.14lf\n", yaxis[0], yaxis[1], yaxis[2]);
//    output::printf("zaxis %20.14lf %20.14lf %20.14lf\n", zaxis[0], zaxis[1], zaxis[2]);

    SharedMatrix frame(new Matrix(3, 3));
    for (i=0; i<3; ++i)
    {
        frame->set(0, i,0, xaxis[i]);
        frame->set(0, i,1, yaxis[i]);
        frame->set(0, i,2, zaxis[i]);
    }

    return frame;
}

std::shared_ptr<PointGroup> Molecule::find_highest_point_group(double tol) const
{
    unsigned char pg_bits = 0;

    typedef void (SymmetryOperation::*symm_func)();

    // The order of the next 2 arrays MUST match!
    unsigned char symm_bit[] =
    {
        SymmOps::C2_z,
        SymmOps::C2_y,
        SymmOps::C2_x,
        SymmOps::i,
        SymmOps::Sigma_xy,
        SymmOps::Sigma_xz,
        SymmOps::Sigma_yz
    };

    symm_func ptrs[] =
    {
        &SymmetryOperation::c2_z,
        &SymmetryOperation::c2_y,
        &SymmetryOperation::c2_x,
        &SymmetryOperation::i,
        &SymmetryOperation::sigma_xy,
        &SymmetryOperation::sigma_xz,
        &SymmetryOperation::sigma_yz
    };

    SymmetryOperation symop;

    int matching_atom = -1;
    // Only needs to detect the 8 symmetry operations
    for (int g=0; g<7; ++g)
    {

        symm_func local_ptr = ptrs[g];

        // Call the function pointer
        (symop.*local_ptr)();

        bool found = true;

        for (int i=0; i<natom(); ++i)
        {
            Vector3 op(symop(0,0), symop(1,1), symop(2,2));
            Vector3 pos = xyz(i) * op;

            if ((matching_atom = atom_at_position2(pos, tol)) >= 0)
            {
                if (atoms_[i]->is_equivalent_to(atoms_[matching_atom]) == false)
                {
                    found = false;
                    break;
                }
            }
            else
            {
                found = false;
                break;
            }
        }

        if (found)
        {
            pg_bits |= symm_bit[g];
        }
    }

    std::shared_ptr<PointGroup> pg = std::shared_ptr<PointGroup>(new PointGroup(pg_bits));

    return pg;
}

void Molecule::reset_point_group(const std::string& pgname)
{
    symmetry_from_input_ = string_to_lower_copy(pgname);
    set_point_group(find_point_group());
}


std::shared_ptr<PointGroup> Molecule::find_point_group(double tol) const
{
    std::shared_ptr<PointGroup> pg = find_highest_point_group(tol);
    const std::string user = symmetry_from_input();

    if (!user.empty())
    {
        // Need to handle the cases that the user only provides C2, C2v, C2h, Cs.
        // These point groups need directionality.

        int end = user.length() - 1;

        bool user_specified_direction = false;
        // Did the user provide directionality? If they did, the last letter would be x, y, or z
        if (user[end] == 'X' || user[end] == 'x' || user[end] == 'Y' || user[end] == 'y' || user[end] == 'Z' || user[end] == 'z')
        {
            // Directionality given, assume the user is smart enough to know what they're doing.
            user_specified_direction = true;
        }

        if (symmetry_from_input() != pg->symbol())
        {
            std::shared_ptr<PointGroup> user(new PointGroup((string_to_upper_copy(symmetry_from_input())).c_str()));

            if (user_specified_direction == true)
            {
                // Assume the user knows what they're doing.

                // Make sure user is subgroup of pg
                if ((pg->bits() & user->bits()) != user->bits())
                {
                    std::stringstream err;

                    err << "User specified point group (" << PointGroups::bits_to_full_name(user->bits()) <<
                        ") is not a subgroup of the highest detected point group (" <<
                        PointGroups::bits_to_full_name(pg->bits()) << ")";
                    throw RuntimeError(err.str());
                }
            }
            else
            {
                unsigned char similars[3];
                char count;

                PointGroups::similar(user->bits(), similars, count);

                int type=0;
                bool found = false;
                for (type=0; type < count; ++type)
                {
                    // If what the user specified and the similar type matches the full point group we've got a
                    // match
                    if ((similars[type] & pg->bits()) == similars[type])
                    {
                        found = true;
                        break;
                    }
                }

                if (found)
                {
                    // Construct a point group object using the found similar
                    user = std::shared_ptr<PointGroup>(new PointGroup(similars[type]));
                }
                else
                {
                    std::stringstream err;

                    err << "User specified point group (" << PointGroups::bits_to_full_name(user->bits()) <<
                        ") is not a subgroup of the highest detected point group (" <<
                        PointGroups::bits_to_full_name(pg->bits()) << "). " <<
                        "If this is because the symmetry increased, try to start the calculation " <<
                        "again from the last geometry, after checking any symmetry-dependent input, " <<
                        "such as DOCC.";
                    throw RuntimeError(err.str().c_str());
                }
            }

            // If we make it here, what the user specified is good.
            pg = user;
        }
    }

    return pg;
}

std::shared_ptr<PointGroup> Molecule::point_group() const
{
    if (!pg_)
        throw RuntimeError("Molecule::point_group: Molecular point group has not been set.");
    return pg_;
}

void Molecule::set_point_group(std::shared_ptr<PointGroup> pg)
{
    pg_ = pg;
    // Call this here, the programmer will forget to call it, as I have many times.
    form_symmetry_information();
}

bool Molecule::has_symmetry_element(Vector3& op, double tol) const
{
    for (int i=0; i<natom(); ++i)
    {
        Vector3 result = xyz(i) * op;
        int atom = atom_at_position2(result, tol);

        if (atom != -1)
        {
            if (!atoms_[atom]->is_equivalent_to(atoms_[i]))
                return false;
        }
        else
            return false;
    }

    return true;
}

void Molecule::symmetrize()
{
    Matrix temp(natom(), 3);
    CharacterTable ct = point_group()->char_table();

    // Obtain atom mapping of atom * symm op to atom
    int **atom_map = compute_atom_map(this);

    // Symmetrize the molecule to remove any noise
    for (int atom=0; atom<natom(); ++atom)
    {
        for (int g=0; g<ct.order(); ++g)
        {

            int Gatom = atom_map[atom][g];

            SymmetryOperation so = ct.symm_operation(g);

            // Full so must be used if molecule is not in standard orientation
            temp.add(0, atom, 0, so(0, 0) * x(Gatom) / ct.order());
            temp.add(0, atom, 0, so(0, 1) * y(Gatom) / ct.order());
            temp.add(0, atom, 0, so(0, 2) * z(Gatom) / ct.order());
            temp.add(0, atom, 1, so(1, 0) * x(Gatom) / ct.order());
            temp.add(0, atom, 1, so(1, 1) * y(Gatom) / ct.order());
            temp.add(0, atom, 1, so(1, 2) * z(Gatom) / ct.order());
            temp.add(0, atom, 2, so(2, 0) * x(Gatom) / ct.order());
            temp.add(0, atom, 2, so(2, 1) * y(Gatom) / ct.order());
            temp.add(0, atom, 2, so(2, 2) * z(Gatom) / ct.order());
        }
    }

    // Delete the atom map.
    delete_atom_map(atom_map, this);
    // Set the geometry to ensure z-matrix variables get updated
    set_geometry(temp);

}

void Molecule::release_symmetry_information()
{
    for (int i=0; i<nunique_; ++i)
    {
        delete[] equiv_[i];
    }
    delete[] equiv_;
    delete[] nequiv_;
    delete[] atom_to_unique_;
    nunique_ = 0;
    equiv_   = 0;
    nequiv_  = 0;
    atom_to_unique_ = 0;
}

void Molecule::form_symmetry_information(double tol)
{
    if (equiv_)
        release_symmetry_information();

    if (natom() == 0)
    {
        nunique_ = 0;
        equiv_   = 0;
        nequiv_  = 0;
        atom_to_unique_ = 0;
        //output::printf("No atoms detected, returning\n");fflush(outfile);
        return;
    }

    nequiv_         = new int[natom()];
    atom_to_unique_ = new int[natom()];
    equiv_          = new int*[natom()];

    if (point_group()->symbol() == "c1")
    {
        nunique_ = natom();
        for (int i=0; i<natom(); ++i)
        {
            nequiv_[i] = 1;
            equiv_[i] = new int[1];
            equiv_[i][0] = i;
            atom_to_unique_[i] = i;
        }
        return;
    }

    // The first atom is always unique
    nunique_           = 1;
    nequiv_[0]         = 1;
    equiv_[0]          = new int[1];
    equiv_[0][0]       = 0;
    atom_to_unique_[0] = 0;

    CharacterTable ct  = point_group()->char_table();

    Vector3 ac;
    SymmetryOperation so;
    Vector3 np;

    // Find the equivalent atoms
    int i;
    for (i=1; i<natom(); ++i)
    {
        ac = xyz(i);
        int i_is_unique = 1;
        int i_equiv = 0;

        // Apply all symmetry ops in the group to the atom
        for (int g=0; g<ct.order(); ++g)
        {
            so = ct.symm_operation(g);
            for (int ii=0; ii<3; ++ii)
            {
                np[ii] = 0;
                for (int jj=0; jj<3; ++jj)
                    np[ii] += so(ii, jj) * ac[jj];
            }

            // See if the transformed atom is equivalent to a
            // unique atom
            for (int j=0; j<nunique_; ++j)
            {
                int unique = equiv_[j][0];
                Vector3 aj(xyz(unique));
                if (np.distance(aj) < tol
                        && Z(unique) == Z(i)
                        && fabs(mass(unique)-mass(i)) < tol)
                {
                    i_is_unique = 0;
                    i_equiv = j;
                    break;
                }
            }
        }
        if (i_is_unique)
        {
            nequiv_[nunique_] = 1;
            equiv_[nunique_] = new int[1];
            equiv_[nunique_][0] = i;
            atom_to_unique_[i] = nunique_;
            nunique_++;
        }
        else
        {
            int *tmp = new int[nequiv_[i_equiv]+1];
            memcpy(tmp, equiv_[i_equiv], nequiv_[i_equiv]*sizeof(int));
            delete[] equiv_[i_equiv];
            equiv_[i_equiv] = tmp;
            equiv_[i_equiv][nequiv_[i_equiv]] = i;
            nequiv_[i_equiv]++;
            atom_to_unique_[i] = i_equiv;
        }
    }

    // The first atom in the equiv list is considered the primary
    // unique atom. Just to make things look pretty, make the
    // atom with the most zeros in its x, y, z coordinate the
    // unique atom. Nothing else should rely on this being done.
    double ztol=1.0e-5;
    for (i=0; i<nunique_; ++i)
    {
        int maxzero = 0;
        int jmaxzero = 0;
        for (int j=0; j<nequiv_[i]; ++j)
        {
            int nzero = 0;
            for (int k=0; k<3; ++k)
            {
                double tmp = equiv_[i][j];
                if (fabs(xyz(tmp, k)) < ztol)
                    nzero++;
            }
            if (nzero > maxzero)
            {
                maxzero = nzero;
                jmaxzero = j;
            }
        }
        int tmp = equiv_[i][jmaxzero];
        equiv_[i][jmaxzero] = equiv_[i][0];
        equiv_[i][0] = tmp;
    }
}

std::string Molecule::sym_label()
{
    if (!pg_) set_point_group(find_point_group());
    return pg_->symbol();
}

char** Molecule::irrep_labels()
{
    if (pg_==NULL) set_point_group(find_point_group());
    int nirreps = pg_->char_table().nirrep();
    char **irreplabel = (char **) malloc(sizeof(char *)*nirreps);
    for (int i=0; i<nirreps; i++)
    {
        irreplabel[i] = (char *) malloc(sizeof(char)*5);
        ::memset(irreplabel[i], 0, sizeof(char)*5);
        strcpy(irreplabel[i],pg_->char_table().gamma(i).symbol());
    }
    return irreplabel;
}

Vector3 Molecule::xyz(int atom) const
{
    return input_units_to_au_ * atoms_[atom]->compute();
}

Vector3 Molecule::fxyz(int atom) const
{
    return input_units_to_au_ * full_atoms_[atom]->compute();
}

double Molecule::xyz(int atom, int _xyz) const
{
    return input_units_to_au_ * atoms_[atom]->compute()[_xyz];
}

const double& Molecule::Z(int atom) const
{
    return atoms_[atom]->Z();
}

double Molecule::fZ(int atom) const
{
    return full_atoms_[atom]->Z();
}

double Molecule::x(int atom) const
{
    return input_units_to_au_ * atoms_[atom]->compute()[0];
}

double Molecule::y(int atom) const
{
    return input_units_to_au_ * atoms_[atom]->compute()[1];
}

double Molecule::z(int atom) const
{
    return input_units_to_au_ * atoms_[atom]->compute()[2];
}

double Molecule::fx(int atom) const
{
    return input_units_to_au_ * full_atoms_[atom]->compute()[0];
}

double Molecule::fy(int atom) const
{
    return input_units_to_au_ * full_atoms_[atom]->compute()[1];
}

double Molecule::fz(int atom) const
{
    return input_units_to_au_ * full_atoms_[atom]->compute()[2];
}

double Molecule::charge(int atom) const
{
    return atoms_[atom]->charge();
}

double Molecule::fcharge(int atom) const
{
    return full_atoms_[atom]->charge();
}

int Molecule::true_atomic_number(int atom) const
{
    Element_to_Z Z;
    Z.load_values();
    return (int)Z[atoms_[atom]->symbol()];
}

int Molecule::ftrue_atomic_number(int atom) const
{
    Element_to_Z Z;
    Z.load_values();
    return (int)Z[full_atoms_[atom]->symbol()];
}

const std::shared_ptr<CoordEntry>& Molecule::atom_entry(int atom) const
{
    return atoms_[atom];
}

double Molecule::fmass(int atom) const
{
    return full_atoms_[atom]->mass();
}

std::string Molecule::flabel(int atom) const
{
    return full_atoms_[atom]->label();
}

void Molecule::set_variable(const std::string &str, double val)
{
    lock_frame_ = false;
    geometry_variables_[str] = val;
    try
    {
        update_geometry();
    }
    catch (...)
    {
        // Update geometry might have added some atoms, delete them to be safe.
        atoms_.clear();
    }
}

double Molecule::get_variable(const std::string &str)
{
    if(geometry_variables_.count(str))
    {
        return geometry_variables_[str];
    }
    else
    {
        throw RuntimeError(str + " not known");
    }
}

bool Molecule::is_variable(const std::string &str) const
{
    return find(all_variables_.begin(), all_variables_.end(), str) != all_variables_.end();
}

std::string Molecule::schoenflies_symbol() const
{
    return point_group()->symbol();
}

// RAK, 4-2012, return true if all atoms correctly map onto other atoms
bool Molecule::valid_atom_map(double tol) const
{
    double np[3];
    SymmetryOperation so;
    CharacterTable ct = point_group()->char_table();

    // loop over all centers
    for (int i=0; i < natom(); i++)
    {
        Vector3 ac(xyz(i));

        // For each operation in the pointgroup, transform the coordinates of
        // center "i" and see which atom it maps into
        for (int g=0; g < ct.order(); g++)
        {
            so = ct.symm_operation(g);

            for (int ii=0; ii < 3; ii++)
            {
                np[ii] = 0;
                for (int jj=0; jj < 3; jj++)
                    np[ii] += so(ii,jj) * ac[jj];
            }

            if (atom_at_position1(np, tol) < 0)
                return false;
        }
    }
    return true;
}


// These two declarations are left here as it's not clear that anyone else will use them:

// Function used by set_full_point_group() to find the max. order of a rotational axis.
int matrix_3d_rotation_Cn(Matrix &coord, Vector3 axis, bool reflect, double TOL, int max_Cn_to_check=-1);

// Function used by set_full_point_group() to scan a given geometry and
// determine if an atom is present at a given location.
bool atom_present_in_geom(Matrix & geom, Vector3 & b, double tol);

bool atom_present_in_geom(Matrix & geom, Vector3 & b, double tol)
{
    for (int i=0; i < geom.nrow(); ++i)
    {
        Vector3 a(geom(i,0),geom(i,1),geom(i,2));
        if (b.distance(a) < tol)
            return true;
    }
    return false;
}

// full_pg_n_ is highest order n in Cn.  0 for atoms or infinity.
void Molecule::set_full_point_group(double zero_tol)
{

    // Get cartesian geometry and put COM at origin
    Matrix geom = geometry();
    Vector3 com = center_of_mass();
    for (int i=0; i<natom(); ++i)
    {
        geom.add(i, 0, -com[0]);
        geom.add(i, 1, -com[1]);
        geom.add(i, 2, -com[2]);
    }

    // Get rotor type
    RotorType rotor = rotor_type(zero_tol);
    //output::printf("\t\tRotor type        : %s\n", RotorTypeList[rotor].c_str());

    // Get the D2h point group from Jet and Ed's code: c1 ci c2 cs d2 c2v c2h d2h
    // and ignore the user-specified subgroup in this case.
    std::shared_ptr<PointGroup> pg = find_highest_point_group(zero_tol);
    std::string d2h_subgroup = pg->symbol();
    //std::string d2h_subgroup = point_group()->symbol();
    //output::printf("d2h_subgroup %s \n", d2h_subgroup.c_str());

    // Check inversion
    Vector3 v3_zero(0, 0, 0);
    bool op_i = has_inversion(v3_zero, zero_tol);
    //output::printf("\t\tInversion symmetry: %s\n", (op_i ? "yes" : "no"));

    int i;
    double dot, phi;
    Vector3 x_axis(1,0,0);
    Vector3 y_axis(0,1,0);
    Vector3 z_axis(0,0,1);
    SharedMatrix test_mat;
    Vector3 rot_axis;

    if (rotor == RT_ATOM)   // atoms
    {
        full_pg_ = PG_ATOM;
        full_pg_n_ = 0;
    }
    else if (rotor == RT_LINEAR)   // linear molecules
    {
        if (op_i)
            full_pg_ = PG_Dinfh;
        else
            full_pg_ = PG_Cinfv;
        full_pg_n_ = 0;
    }
    else if (rotor == RT_SPHERICAL_TOP)   // spherical tops
    {
        if (!op_i)   // The only spherical top without inversion is Td.
        {
            full_pg_ = PG_Td;
            full_pg_n_ = 3;
        }
        else   // Oh or Ih ?
        {
            // Oh has a S4 and should be oriented properly already.
            test_mat = geom.matrix_3d_rotation(z_axis, pc_pi/2, true);
            bool op_symm = geom.equal_but_for_row_order(test_mat, zero_tol);
            //output::printf("\t\tS4z : %s\n", (op_symm ? "yes" : "no"));

            if (op_symm)
            {
                full_pg_ = PG_Oh;
                full_pg_n_ = 4;
            }
            else
            {
                full_pg_ = PG_Ih;
                full_pg_n_ = 5;
            }
        }
    }
    else if (rotor == RT_ASYMMETRIC_TOP)   // asymmetric tops cannot exceed D2h, right?
    {

        if (d2h_subgroup == "c1")
        {
            full_pg_ = PG_C1;
            full_pg_n_ = 1;
        }
        else if (d2h_subgroup == "ci")
        {
            full_pg_ = PG_Ci;
            full_pg_n_ = 1;
        }
        else if (d2h_subgroup == "c2")
        {
            full_pg_ = PG_Cn;
            full_pg_n_ = 2;
        }
        else if (d2h_subgroup == "cs")
        {
            full_pg_ = PG_Cs;
            full_pg_n_ = 1;
        }
        else if (d2h_subgroup == "d2")
        {
            full_pg_ = PG_Dn;
            full_pg_n_ = 2;
        }
        else if (d2h_subgroup == "c2v")
        {
            full_pg_ = PG_Cnv;
            full_pg_n_ = 2;
        }
        else if (d2h_subgroup == "c2h")
        {
            full_pg_ = PG_Cnh;
            full_pg_n_ = 2;
        }
        else if (d2h_subgroup == "d2h")
        {
            full_pg_ = PG_Dnh;
            full_pg_n_ = 2;
        }
        else
           output::printf("\t\tWarning: Cannot determine point group.\n");
    }
    else if (rotor == RT_SYMMETRIC_TOP)
    {

        // Find principal axis that is unique and make it z-axis.
        SharedMatrix It(inertia_tensor());
        Vector I_evals(3);
        SharedMatrix I_evects(new Matrix(3, 3));
        It->diagonalize(I_evects, I_evals, ascending);
        // I_evects->print_out();
        // output::printf("I_evals %15.10lf %15.10lf %15.10lf\n", I_evals[0], I_evals[1], I_evals[2]);

        int unique_axis = 1;
        if (fabs(I_evals[0] - I_evals[1]) < zero_tol)
            unique_axis = 2;
        else if (fabs(I_evals[1] - I_evals[2]) < zero_tol)
            unique_axis = 0;

        // Compute angle between unique axis and the z-axis
        // Returned eigenvectors appear to be columns (in Fortan style) ?!
        Vector3 old_axis(I_evects->get(0,unique_axis),
                         I_evects->get(1,unique_axis),
                         I_evects->get(2,unique_axis));

        dot = z_axis.dot(old_axis);
        if (fabs(dot-1) < 1.0e-10)
            phi = 0.0;
        else if (fabs(dot+1) < 1.0e-10)
            phi = pc_pi;
        else
            phi = acos(dot);

        // Rotate geometry to put unique axis on the z-axis, if it isn't already.
        if (fabs(phi) > 1.0e-14)
        {
            rot_axis = z_axis.cross(old_axis);
            test_mat = geom.matrix_3d_rotation(rot_axis, phi, false);
            //output::printf("Rotating by %lf to get principal axis on z-axis.\n", phi);
            geom.copy(test_mat);
        }

        //output::printf("Geometry to analyze - principal axis on z-axis:\n");
        //for (i=0; i<natom(); ++i)
        //output::printf("%20.15lf %20.15lf %20.15lf\n", geom(i,0), geom(i,1), geom(i,2));
        //output::printf("\n");

        // Determine order Cn and Sn of principal axis.
        int Cn_z = matrix_3d_rotation_Cn(geom, z_axis, false, zero_tol);
        //output::printf("\t\tHighest rotation axis (Cn_z) : %d\n", Cn_z);

        int Sn_z = matrix_3d_rotation_Cn(geom, z_axis, true, zero_tol);
        //output::printf("\t\tHighest rotation axis (Sn_z) : %d\n", Sn_z);

        // Check for sigma_h (xy plane).
        bool op_sigma_h = false;
        for (i=0; i<natom(); ++i)
        {
            if (fabs(geom(i,2)) < zero_tol)
                continue; // atom is in xy plane
            else
            {
                Vector3 test_atom(geom(i,0), geom(i,1), -1*geom(i,2));
                if (!atom_present_in_geom(geom, test_atom, zero_tol))
                    break;
            }
        }
        if (i == natom())
            op_sigma_h = true;
        //output::printf("\t\t sigma_h : %s\n", (op_sigma_h ? "yes" : "no"));

        // Rotate one off-axis atom to the yz plane and check for sigma_v's.
        int pivot_atom_i;
        for (i=0; i<natom(); ++i)
        {
            double dist_from_z = sqrt(geom(i,0)*geom(i,0)+geom(i,1)*geom(i,1));
            if (fabs(dist_from_z) > zero_tol)
            {
                pivot_atom_i = i;
                break;
            }
        }
        if (pivot_atom_i == natom())
            throw RuntimeError("Not a linear molecule but could not find off-axis atom.");

        // Rotate around z-axis to put pivot atom in the yz plane
        Vector3 xy_point(geom(pivot_atom_i,0), geom(pivot_atom_i,1), 0);

        xy_point.normalize();
        dot = y_axis.dot(xy_point);
        if (fabs(dot-1) < 1.0e-10)
            phi = 0.0;
        else if (fabs(dot+1) < 1.0e-10)
            phi = pc_pi;
        else
            phi = acos(dot);

        bool is_D = false;
        if (fabs(phi) > 1.0e-14)
        {
            test_mat = geom.matrix_3d_rotation(z_axis, phi, false);
            //output::printf("Rotating by %8.3e to get atom %d in yz-plane.\n", phi, pivot_atom_i+1);
            geom.copy(test_mat);
        }

        // Check for sigma_v (yz plane).
        bool op_sigma_v = false;
        for (i=0; i<natom(); ++i)
        {
            if (fabs(geom(i,0)) < zero_tol)
                continue; // atom is in yz plane
            else
            {
                Vector3 test_atom(-1*geom(i,0), geom(i,1), geom(i,2));
                if (!atom_present_in_geom(geom, test_atom, zero_tol))
                    break;
            }
        }
        if (i == natom())
            op_sigma_v = true;
        //output::printf("\t\tsigma_v : %s\n", (op_sigma_v ? "yes" : "no"));

        //output::printf("geom to analyze - one atom in yz plane\n");
        //for (i=0; i<natom(); ++i)
        //output::printf("%20.15lf %20.15lf %20.15lf\n", geom(i,0), geom(i,1), geom(i,2));
        //output::printf("\n");

        // Check for perpendicular C2's.
        // Loop through pairs of atoms to find c2 axis candidates.
        for (i=0; i<natom(); ++i)
        {
            Vector3 A(geom(i,0), geom(i,1), geom(i,2));
            double AdotA = A.dot(A);
            for (int j=0; j<i; ++j)
            {

                if (Z(i) != Z(j)) continue; // ensure same atomic number

                Vector3 B(geom(j,0), geom(j,1), geom(j,2)); // ensure same distance from com
                if (fabs(AdotA - B.dot(B)) > 1.0e-6) continue; // loose check

                // Use sum of atom vectors as axis if not 0.
                Vector3 axis= A+B;
                if (axis.norm() < 1.0e-12) continue;
                axis.normalize();

                // Check if axis is perpendicular to z-axis.
                if (fabs(axis.dot(z_axis)) > 1.0e-6) continue;

                // Do the thorough check for C2.
                if (matrix_3d_rotation_Cn(geom, axis, false, zero_tol, 2) == 2)
                    is_D = true;
            }
        }
        //output::printf("\t\tperp. C2's :  %s\n", (is_D ? "yes" : "no"));

        // Now assign point groups!  Sn first.
        if (Sn_z == 2 * Cn_z && !is_D)
        {
            full_pg_   = PG_Sn;
            full_pg_n_ = Sn_z;
            return;
        }

        if (is_D)    // has perpendicular C2's
        {
            if (op_sigma_h && op_sigma_v)   // Dnh : Cn, nC2, sigma_h, nSigma_v
            {
                full_pg_   = PG_Dnh;
                full_pg_n_ = Cn_z;
            }
            else if (Sn_z == 2*Cn_z)   // Dnd : Cn, nC2, S2n axis coincident with Cn
            {
                full_pg_   = PG_Dnd;
                full_pg_n_ = Cn_z;
            }
            else                       // Dn : Cn, nC2
            {
                full_pg_   = PG_Dn;
                full_pg_n_ = Cn_z;
            }
        }
        else        // lacks perpendicular C2's
        {
            if (op_sigma_h && Sn_z == Cn_z)  // Cnh : Cn, sigma_h, Sn coincident with Cn
            {
                full_pg_   = PG_Cnh;
                full_pg_n_ = Cn_z;
            }
            else if (op_sigma_v)             // Cnv : Cn, nCv
            {
                full_pg_   = PG_Cnv;
                full_pg_n_ = Cn_z;
            }
            else                             // Cn  : Cn
            {
                full_pg_   = PG_Cn;
                full_pg_n_ = Cn_z;
            }
        }
    } // symmetric top

    return;
}

/*
** @brief Find maximum n in Cn around given axis, i.e., the highest-order rotation axis.
**
** @param coord Matrix    : points to rotate - column dim is 3
** @param axis  Vector3   : axis around which to rotate, does not need to be normalized
** @param bool  reflect   : if true, really look for Sn not Cn
** @returns n
*/
int matrix_3d_rotation_Cn(Matrix &coord, Vector3 axis, bool reflect, double TOL, int max_Cn_to_check)
{
    int max_possible;
    if (max_Cn_to_check == -1)     // default
        max_possible = coord.nrow(); // Check all atoms. In future, make more intelligent.
    else
        max_possible = max_Cn_to_check;

    int Cn = 1; // C1 is there for sure
    SharedMatrix rotated_mat;
    bool present;

    for (int n=2; n<max_possible+1; ++n)
    {
        rotated_mat = coord.matrix_3d_rotation(axis, 2*pc_pi/n, reflect);
        present = coord.equal_but_for_row_order(rotated_mat, TOL);

        if (present)
            Cn = n;
    }
    return Cn;
}


// Return point group name such as D3d or S8 in string form, with the 'n'
// replaced by an integer.
std::string Molecule::full_point_group() const
{

    string pg_with_n = FullPointGroupList[full_pg_];

    // These don't need changes - have no 'n'.
    if (pg_with_n == "D_inf_h" || pg_with_n == "C_inf_v" ||
            pg_with_n == "C1"      || pg_with_n == "Cs"      ||
            pg_with_n == "Ci"      || pg_with_n == "Td"      ||
            pg_with_n == "Oh"      || pg_with_n == "Ih" )
        return pg_with_n;

    stringstream n_integer;
    n_integer << full_pg_n_;

    // Replace 'n'.  It can only appear once.
    size_t start_pos = pg_with_n.find("n");

    pg_with_n.replace(start_pos, n_integer.str().length(), n_integer.str());

    return  pg_with_n;
}


// Free functions
int **compute_atom_map(const Molecule* molecule)
{
    // grab references to the Molecule
    const Molecule& mol = *molecule;

    // create the character table for the point group
    CharacterTable ct = mol.point_group()->char_table();

    int natom = mol.natom();
    int ng = ct.order();
    int **atom_map;
    atom_map = new int*[natom];
    for (int i=0; i < natom; i++)
        atom_map[i] = new int[ng];

    double np[3];
    SymmetryOperation so;

    // loop over all centers
    for (int i=0; i < natom; i++) {
        Vector3 ac(mol.xyz(i));

        // then for each symop in the pointgroup, transform the coordinates of
        // center "i" and see which atom it maps into
        for (int g=0; g < ng; g++) {
            so = ct.symm_operation(g);

            for (int ii=0; ii < 3; ii++) {
                np[ii] = 0;
                for (int jj=0; jj < 3; jj++)
                    np[ii] += so(ii,jj) * ac[jj];
            }

            atom_map[i][g] = mol.atom_at_position1(np, 0.05);
            if (atom_map[i][g] < 0) {
                std::stringstream ss;
                ss << "ERROR: Symmetry operation " << g << " did not map atom " << i+1 << " to another atom\n";
                ss << "  -> attempted to find atom at\n"
                   << np[0] << " " << np[1] << " " <<  np[2] << "\n";
                throw RuntimeError(ss.str());
            }
        }
    }

    return atom_map;
}

int **compute_atom_map(const std::shared_ptr<Molecule> &molecule)
{
    return compute_atom_map(molecule.get());
}

void delete_atom_map(int **atom_map, const Molecule* molecule)
{
    if (atom_map) {
        int natom = molecule->natom();
        for (int i=0; i < natom; i++)
            delete[] atom_map[i];
        delete[] atom_map;
    }
}

void delete_atom_map(int **atom_map, const std::shared_ptr<Molecule> &molecule)
{
    delete_atom_map(atom_map, molecule.get());
}


} // end namespace panache


