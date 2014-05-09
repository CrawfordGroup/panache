#include <cctype>
#include <cmath> // fabs
#include "CharacterTable.h"
#include "Exception.h"
#include "PhysConst.h"
#include "Output.h"

using namespace std;

////////////////////////////////////////////////////////////////////////

namespace panache {


CharacterTable::CharacterTable()
    : nt(0), pg(PointGroups::C1), nirrep_(0), gamma_(0), symop(0), _inv(0), symb(0),
      bits_(0)
{
}

CharacterTable::CharacterTable(const CharacterTable& ct)
    : nt(0), pg(PointGroups::C1), nirrep_(0), gamma_(0), symop(0), _inv(0), symb(0),
      bits_(0)
{
    *this = ct;
}

CharacterTable::~CharacterTable()
{
    if (gamma_) delete[] gamma_; gamma_=0;
    if (symop) delete[] symop; symop=0;
    if (_inv) delete[] _inv; _inv=0;
    nt=nirrep_=0;
}

CharacterTable&
CharacterTable::operator=(const CharacterTable& ct)
{
    nt=ct.nt; pg=ct.pg; nirrep_=ct.nirrep_;

    symb = ct.symb;

    if (gamma_) delete[] gamma_; gamma_=0;
    if (ct.gamma_) {
        gamma_ = new IrreducibleRepresentation[nirrep_];
        for (int i=0; i < nirrep_; i++) {
            gamma_[i].init();
            gamma_[i] = ct.gamma_[i];
        }
    }

    if (symop)
        delete[] symop;
    symop=0;

    if (ct.symop) {
        symop = new SymmetryOperation[nirrep_];
        for (int i=0; i < nirrep_; i++) {
            symop[i] = ct.symop[i];
        }
    }

    if (_inv)
        delete[] _inv;
    _inv=0;

    if (ct._inv) {
        _inv = new int[nirrep_];
        memcpy(_inv,ct._inv,sizeof(int)* nirrep_);
    }

    return *this;
}

/*
void CharacterTable::print(FILE *out) const
{
    if (!nirrep_) return;

    int i;

    output::printf("  point group %s\n\n", symb.c_str());

    for (i=0; i < nirrep_; i++)
        gamma_[i].print(out);

    output::printf("\n  symmetry operation matrices:\n\n");
    for (i=0; i < nirrep_; i++)
        symop[i].print(out);

    output::printf("\n  inverse symmetry operation matrices:\n\n");
    for (i=0; i < nirrep_; i++)
        symop[inverse(i)].print(out);
}
*/

void CharacterTable::common_init()
{
    // first parse the point group symbol, this will give us the order of the
    // point group(g), the type of point group (pg), the order of the principle
    // rotation axis (nt), and the number of irreps (nirrep_)

    if (!symb.length()) {
        throw RuntimeError("CharacterTable::CharacterTable: null point group");
    }

    if (make_table() < 0) {
        throw RuntimeError("CharacterTable::CharacterTable: could not make table");
    }
}

CharacterTable::CharacterTable(const std::string& cpg)
    : nt(0), pg(PointGroups::C1), nirrep_(0), gamma_(0), symop(0), _inv(0), symb(cpg),
      bits_(0)
{
    // Check the symbol coming in
    if (!PointGroups::full_name_to_bits(cpg, bits_)) {
        output::printf("CharacterTable: Invalid point group name: %s\n", cpg.c_str());
        throw RuntimeError("CharacterTable: Invalid point group name provided.");
    }
    common_init();
}

CharacterTable::CharacterTable(unsigned char bits)
    : bits_(bits)
{
    symb = PointGroups::bits_to_basic_name(bits);
    common_init();
}

unsigned char CharacterTable::bits()
{
    return bits_;
}

int CharacterTable::make_table()
{
    int i,j,ei,gi;
    char label[4];

    switch (bits_) {
    case PointGroups::C1:
        nirrep_ = 1;
        nt = 1;
        break;
    case PointGroups::CsX:
    case PointGroups::CsY:
    case PointGroups::CsZ:
    case PointGroups::Ci:
        nirrep_ = 2;
        nt = 1;
        break;
    case PointGroups::C2X:
    case PointGroups::C2Y:
    case PointGroups::C2Z:
        nirrep_ = 2;
        nt = 2;
        break;
    case PointGroups::C2hX:
    case PointGroups::C2hY:
    case PointGroups::C2hZ:
    case PointGroups::C2vX:
    case PointGroups::C2vY:
    case PointGroups::C2vZ:
    case PointGroups::D2:
        nirrep_ = 4;
        nt = 2;
        break;
    case PointGroups::D2h:
        nirrep_ = 8;
        nt = 2;
        break;
    default:
        throw RuntimeError("Should not have receached here!");
    }

    if (!nirrep_) return 0;

    gamma_ = new IrreducibleRepresentation[nirrep_];
    symop = new SymmetryOperation[nirrep_];
    SymmetryOperation so;

    _inv = new int[nirrep_];

    // the angle to rotate about the principal axis
    double theta = (nt) ? 2.0*M_PI/nt : 2.0*M_PI;

    // Handle irreducible representations:
    switch (bits_) {
    case PointGroups::C1:
        // no symmetry case
        gamma_[0].init(1,1,"A");
        gamma_[0].nrot_ = 3;
        gamma_[0].ntrans_ = 3;
        gamma_[0].rep[0][0][0] = 1.0;

        break;

    case PointGroups::CsX: // reflection through the yz plane
    case PointGroups::CsY: // reflection through the xz plane
    case PointGroups::CsZ: // reflection through the xy plane
        gamma_[0].init(2,1,"A'","Ap");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].nrot_=1;
        gamma_[0].ntrans_=2;

        gamma_[1].init(2,1,"A\"","App");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] = -1.0;
        gamma_[1].nrot_=2;
        gamma_[1].ntrans_=1;

        break;

    case PointGroups::Ci:
        // equivalent to S2 about the z axis
        gamma_[0].init(2,1,"Ag");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].nrot_=3;

        gamma_[1].init(2,1,"Au");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] = -1.0;
        gamma_[1].ntrans_=3;

        break;

    case PointGroups::C2X:
    case PointGroups::C2Y:
    case PointGroups::C2Z:
        gamma_[0].init(2,1,"A","A");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].nrot_=1;
        gamma_[0].ntrans_=1;

        gamma_[1].init(2,1,"B","B");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] = -1.0;
        gamma_[1].nrot_=2;
        gamma_[1].ntrans_=2;

        break;

    case PointGroups::C2hX:
    case PointGroups::C2hY:
    case PointGroups::C2hZ:
        gamma_[0].init(4,1,"Ag","Ag");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].rep[2][0][0] = 1.0;
        gamma_[0].rep[3][0][0] = 1.0;
        gamma_[0].nrot_=1;
        gamma_[0].ntrans_=0;

        gamma_[1].init(4,1,"Bg","Bg");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] = -1.0;
        gamma_[1].rep[2][0][0] =  1.0;
        gamma_[1].rep[3][0][0] = -1.0;
        gamma_[1].nrot_=2;
        gamma_[1].ntrans_=0;

        gamma_[2].init(4, 1,"Au","Au");
        gamma_[2].rep[0][0][0] =  1.0;
        gamma_[2].rep[1][0][0] =  1.0;
        gamma_[2].rep[2][0][0] = -1.0;
        gamma_[2].rep[3][0][0] = -1.0;
        gamma_[2].nrot_=0;
        gamma_[2].ntrans_=1;

        gamma_[3].init(4, 1,"Bu","Bu");
        gamma_[3].rep[0][0][0] =  1.0;
        gamma_[3].rep[1][0][0] = -1.0;
        gamma_[3].rep[2][0][0] = -1.0;
        gamma_[3].rep[3][0][0] =  1.0;
        gamma_[3].nrot_=0;
        gamma_[3].ntrans_=2;

        break;

    case PointGroups::C2vX:
    case PointGroups::C2vY:
    case PointGroups::C2vZ:
        gamma_[0].init(4,1,"A1","A1");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].rep[2][0][0] = 1.0;
        gamma_[0].rep[3][0][0] = 1.0;
        gamma_[0].nrot_=0;
        gamma_[0].ntrans_=1;

        gamma_[1].init(4,1,"A2","A2");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] =  1.0;
        gamma_[1].rep[2][0][0] = -1.0;
        gamma_[1].rep[3][0][0] = -1.0;
        gamma_[1].nrot_=1;
        gamma_[1].ntrans_=0;

        gamma_[2].init(4, 1,"B1","B1");
        gamma_[2].rep[0][0][0] =  1.0;
        gamma_[2].rep[1][0][0] = -1.0;
        gamma_[2].rep[2][0][0] =  1.0;
        gamma_[2].rep[3][0][0] = -1.0;
        gamma_[2].nrot_=1;
        gamma_[2].ntrans_=1;

        gamma_[3].init(4, 1,"B2","B2");
        gamma_[3].rep[0][0][0] =  1.0;
        gamma_[3].rep[1][0][0] = -1.0;
        gamma_[3].rep[2][0][0] = -1.0;
        gamma_[3].rep[3][0][0] =  1.0;
        gamma_[3].nrot_=1;
        gamma_[3].ntrans_=1;

        break;

    case PointGroups::D2:
        gamma_[0].init(4,1,"A","A");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].rep[2][0][0] = 1.0;
        gamma_[0].rep[3][0][0] = 1.0;
        gamma_[0].nrot_=0;
        gamma_[0].ntrans_=0;

        gamma_[1].init(4,1,"B1","B1");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] =  1.0;
        gamma_[1].rep[2][0][0] = -1.0;
        gamma_[1].rep[3][0][0] = -1.0;
        gamma_[1].nrot_=1;
        gamma_[1].ntrans_=1;

        gamma_[2].init(4, 1,"B2","B2");
        gamma_[2].rep[0][0][0] =  1.0;
        gamma_[2].rep[1][0][0] = -1.0;
        gamma_[2].rep[2][0][0] =  1.0;
        gamma_[2].rep[3][0][0] = -1.0;
        gamma_[2].nrot_=1;
        gamma_[2].ntrans_=1;

        gamma_[3].init(4, 1,"B3","B3");
        gamma_[3].rep[0][0][0] =  1.0;
        gamma_[3].rep[1][0][0] = -1.0;
        gamma_[3].rep[2][0][0] = -1.0;
        gamma_[3].rep[3][0][0] =  1.0;
        gamma_[3].nrot_=1;
        gamma_[3].ntrans_=1;

        break;

    case PointGroups::D2h:

        gamma_[0].init(8,1,"Ag","Ag");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].rep[2][0][0] = 1.0;
        gamma_[0].rep[3][0][0] = 1.0;
        gamma_[0].rep[4][0][0] = 1.0;
        gamma_[0].rep[5][0][0] = 1.0;
        gamma_[0].rep[6][0][0] = 1.0;
        gamma_[0].rep[7][0][0] = 1.0;
        gamma_[0].nrot_=0;
        gamma_[0].ntrans_=0;

        gamma_[1].init(8,1,"B1g","B1g");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] =  1.0;
        gamma_[1].rep[2][0][0] = -1.0;
        gamma_[1].rep[3][0][0] = -1.0;
        gamma_[1].rep[4][0][0] =  1.0;
        gamma_[1].rep[5][0][0] =  1.0;
        gamma_[1].rep[6][0][0] = -1.0;
        gamma_[1].rep[7][0][0] = -1.0;
        gamma_[1].nrot_=1;
        gamma_[1].ntrans_=0;

        gamma_[2].init(8,1,"B2g","B2g");
        gamma_[2].rep[0][0][0] =  1.0;
        gamma_[2].rep[1][0][0] = -1.0;
        gamma_[2].rep[2][0][0] =  1.0;
        gamma_[2].rep[3][0][0] = -1.0;
        gamma_[2].rep[4][0][0] =  1.0;
        gamma_[2].rep[5][0][0] = -1.0;
        gamma_[2].rep[6][0][0] =  1.0;
        gamma_[2].rep[7][0][0] = -1.0;
        gamma_[2].nrot_=1;
        gamma_[2].ntrans_=0;

        gamma_[3].init(8,1,"B3g","B3g");
        gamma_[3].rep[0][0][0] =  1.0;
        gamma_[3].rep[1][0][0] = -1.0;
        gamma_[3].rep[2][0][0] = -1.0;
        gamma_[3].rep[3][0][0] =  1.0;
        gamma_[3].rep[4][0][0] =  1.0;
        gamma_[3].rep[5][0][0] = -1.0;
        gamma_[3].rep[6][0][0] = -1.0;
        gamma_[3].rep[7][0][0] =  1.0;
        gamma_[3].nrot_=1;
        gamma_[3].ntrans_=0;

        gamma_[4].init(8,1,"Au","Au");
        gamma_[4].rep[0][0][0] =  1.0;
        gamma_[4].rep[1][0][0] =  1.0;
        gamma_[4].rep[2][0][0] =  1.0;
        gamma_[4].rep[3][0][0] =  1.0;
        gamma_[4].rep[4][0][0] = -1.0;
        gamma_[4].rep[5][0][0] = -1.0;
        gamma_[4].rep[6][0][0] = -1.0;
        gamma_[4].rep[7][0][0] = -1.0;
        gamma_[4].nrot_=0;
        gamma_[4].ntrans_=0;

        gamma_[5].init(8,1,"B1u","B1u");
        gamma_[5].rep[0][0][0] =  1.0;
        gamma_[5].rep[1][0][0] =  1.0;
        gamma_[5].rep[2][0][0] = -1.0;
        gamma_[5].rep[3][0][0] = -1.0;
        gamma_[5].rep[4][0][0] = -1.0;
        gamma_[5].rep[5][0][0] = -1.0;
        gamma_[5].rep[6][0][0] =  1.0;
        gamma_[5].rep[7][0][0] =  1.0;
        gamma_[5].nrot_=0;
        gamma_[5].ntrans_=1;

        gamma_[6].init(8,1,"B2u","B2u");
        gamma_[6].rep[0][0][0] =  1.0;
        gamma_[6].rep[1][0][0] = -1.0;
        gamma_[6].rep[2][0][0] =  1.0;
        gamma_[6].rep[3][0][0] = -1.0;
        gamma_[6].rep[4][0][0] = -1.0;
        gamma_[6].rep[5][0][0] =  1.0;
        gamma_[6].rep[6][0][0] = -1.0;
        gamma_[6].rep[7][0][0] =  1.0;
        gamma_[6].nrot_=0;
        gamma_[6].ntrans_=1;

        gamma_[7].init(8,1,"B3u","B3u");
        gamma_[7].rep[0][0][0] =  1.0;
        gamma_[7].rep[1][0][0] = -1.0;
        gamma_[7].rep[2][0][0] = -1.0;
        gamma_[7].rep[3][0][0] =  1.0;
        gamma_[7].rep[4][0][0] = -1.0;
        gamma_[7].rep[5][0][0] =  1.0;
        gamma_[7].rep[6][0][0] =  1.0;
        gamma_[7].rep[7][0][0] = -1.0;
        gamma_[7].nrot_=0;
        gamma_[7].ntrans_=1;

        break;
    }

    // Handle symmetry operations
    symop[0].E();

    switch (bits_) {
    case PointGroups::C1:

        // nothing to do.

        break;

    case PointGroups::Ci:

        symop[1].i();

        break;

    case PointGroups::CsX: // reflection through the yz plane

        symop[1].sigma_yz();

        break;

    case PointGroups::CsY: // reflection through the xz plane

        symop[1].sigma_xz();

        break;

    case PointGroups::CsZ: // reflection through the xy plane

        symop[1].sigma_xy();

        break;

    case PointGroups::C2X:

        symop[1].c2_x();

        break;

    case PointGroups::C2Y:

        symop[1].c2_y();

        break;

    case PointGroups::C2Z:

        symop[1].rotation(2);

        break;

    case PointGroups::C2hX:

        symop[1].c2_x();
        symop[2].i();
        symop[3].sigma_yz();

        break;

    case PointGroups::C2hY:

        symop[1].c2_y();
        symop[2].i();
        symop[3].sigma_xz();

        break;

    case PointGroups::C2hZ:

        symop[1].rotation(2);
        symop[2].i();
        symop[3].sigma_xy();

        break;

    case PointGroups::C2vX:

        symop[1].c2_x();
        symop[2].sigma_xy();
        symop[3].sigma_xz();

        break;

    case PointGroups::C2vY:

        symop[1].c2_y();
        symop[2].sigma_xy();
        symop[3].sigma_yz();

        break;

    case PointGroups::C2vZ:

        symop[1].rotation(2);
        symop[2].sigma_xz();
        symop[3].sigma_yz();

        break;

    case PointGroups::D2:

        symop[1].rotation(2);
        symop[2].c2_y();
        symop[3].c2_x();

        break;

    case PointGroups::D2h:

        symop[1].rotation(2);
        symop[2].c2_y();
        symop[3].c2_x();
        symop[4].i();
        symop[5].sigma_xy();
        symop[6].sigma_xz();
        symop[7].sigma_yz();

        break;

    default:
        return -1;

    }

    // now find the inverse of each symop
    for (gi=0; gi < nirrep_; gi++) {
        int gj;
        for (gj=0; gj < nirrep_; gj++) {
            so = symop[gi].operate(symop[gj]);

            // is so a unit matrix?
            if (fabs(1.0-so[0][0]) < 1.0e-8 &&
                    fabs(1.0-so[1][1]) < 1.0e-8 &&
                    fabs(1.0-so[2][2]) < 1.0e-8) break;
        }

        if (gj==nirrep_) {
            // ExEnv::err0() << indent
            //      << "make_table: uh oh, can't find inverse of " << gi << endl;
            // abort();
            throw RuntimeError("make_table: uh oh, can't find inverse");
        }

        _inv[gi] = gj;
    }

    // Check the bits of the operator make sure they make what
    // we were given.
    unsigned char sym_bits = 0;
    for (int i=0; i<nirrep_; ++i) {
        sym_bits |= symop[i].bit();
    }

    if (sym_bits != bits_)
        throw RuntimeError("make_table: Symmetry operators did not match the point group given.");

    return 0;
}

} // end namespace panache

