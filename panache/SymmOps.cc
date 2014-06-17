#include <cstring> // for memcpy
#include "SymmOps.h"
#include "Exception.h"
#include "Output.h"


namespace panache {
namespace PointGroups
{

void similar(unsigned char bits, unsigned char* sim, char& cnt)
{
    static unsigned char cs[3] = { CsX, CsY, CsZ };
    static unsigned char c2v[3] = { C2vZ, C2vY, C2vX };
    static unsigned char c2h[3] = { C2hZ, C2hY, C2hX };
    static unsigned char c2[3]  = { C2Z, C2Y, C2X };
    static unsigned char d2h[3] = { D2h, 0, 0 };
    static unsigned char d2[3]  = { D2, 0, 0 };
    static unsigned char ci[3]  = { Ci, 0, 0 };
    static unsigned char c1[3]  = { C1, 0, 0 };

    switch (bits) {
    case CsX:
    case CsY:
    case CsZ:
        memcpy(sim, cs, sizeof(char)*3);
        cnt = 3;
        break;
    case C2vZ:
    case C2vY:
    case C2vX:
        memcpy(sim, c2v, sizeof(char)*3);
        cnt = 3;
        break;
    case C2hZ:
    case C2hY:
    case C2hX:
        memcpy(sim, c2h, sizeof(char)*3);
        cnt = 3;
        break;
    case C2Z:
    case C2Y:
    case C2X:
        memcpy(sim, c2, sizeof(char)*3);
        cnt = 3;
        break;
    case D2h:
        memcpy(sim, d2h, sizeof(char)*1);
        cnt = 1;
        break;
    case Ci:
        memcpy(sim, ci, sizeof(char)*1);
        cnt = 1;
        break;
    case C1:
        memcpy(sim, c1, sizeof(char)*1);
        cnt = 1;
        break;
    case D2:
        memcpy(sim, d2, sizeof(char)*1);
        cnt = 1;
        break;
    default:
        throw RuntimeError("Should not have reaced here.");
    }
}

bool full_name_to_bits(const std::string& pg, unsigned char &bits)
{
    bool retvalue = true;

    if (pg == "C1")
        bits = PointGroups::C1;
    else if (pg == "Ci")
        bits = PointGroups::Ci;
    else if (pg == "C2(x)" || pg == "C2x" || pg == "C2_x")
        bits = PointGroups::C2X;
    else if (pg == "C2(y)" || pg == "C2y" || pg == "C2_y")
        bits = PointGroups::C2Y;
    else if (pg == "C2(z)" || pg == "C2z" || pg == "C2_z")
        bits = PointGroups::C2Z;
    else if (pg == "Cs(x)" || pg == "Csx" || pg == "Cs_x")
        bits = PointGroups::CsX;
    else if (pg == "Cs(y)" || pg == "Csy" || pg == "Cs_y")
        bits = PointGroups::CsY;
    else if (pg == "Cs(z)" || pg == "Csz" || pg == "Cs_z")
        bits = PointGroups::CsZ;
    else if (pg == "D2")
        bits = PointGroups::D2;
    else if (pg == "C2v(X)" || pg == "C2vx" || pg == "C2v_x")
        bits = PointGroups::C2vX;
    else if (pg == "C2v(Y)" || pg == "C2vy" || pg == "C2v_y")
        bits = PointGroups::C2vY;
    else if (pg == "C2v(Z)" || pg == "C2vz" || pg == "C2v_z")
        bits = PointGroups::C2vZ;
    else if (pg == "C2h(X)" || pg == "C2hx" || pg == "C2h_x")
        bits = PointGroups::C2hX;
    else if (pg == "C2h(Y)" || pg == "C2hy" || pg == "C2h_y")
        bits = PointGroups::C2hY;
    else if (pg == "C2h(Z)" || pg == "C2hz" || pg == "C2h_z")
        bits = PointGroups::C2hZ;
    else if (pg == "D2h")
        bits = PointGroups::D2h;

    // Ok, the user gave us Cs, C2v, C2h, C2, but no directionality
    else if (pg == "Cs")
        bits = PointGroups::CsX;
    else if (pg == "C2v")
        bits = PointGroups::C2vZ;
    else if (pg == "C2h")
        bits = PointGroups::C2hZ;
    else if (pg == "C2")
        bits = PointGroups::C2Z;

    else
        retvalue = false;

    return retvalue;
}

const char* bits_to_full_name(unsigned char bits)
{
    switch(bits) {
    case PointGroups::C1:
        return "C1";
    case PointGroups::Ci:
        return "Ci";
    case PointGroups::C2X:
        return "C2(x)";
    case PointGroups::C2Y:
        return "C2(y)";
    case PointGroups::C2Z:
        return "C2(z)";
    case PointGroups::CsZ:
        return "Cs(Z)";
    case PointGroups::CsY:
        return "Cs(Y)";
    case PointGroups::CsX:
        return "Cs(X)";
    case PointGroups::D2:
        return "D2";
    case PointGroups::C2vX:
        return "C2v(X)";
    case PointGroups::C2vY:
        return "C2v(Y)";
    case PointGroups::C2vZ:
        return "C2v(Z)";
    case PointGroups::C2hX:
        return "C2h(X)";
    case PointGroups::C2hY:
        return "C2h(Y)";
    case PointGroups::C2hZ:
        return "C2h(Z)";
    case PointGroups::D2h:
        return "D2h";
    default:
        output::printf("Unrecognized point group bits: %d\n", bits);
        throw RuntimeError("Unrecognized point group bits");
    }
}

const char* bits_to_basic_name(unsigned char bits)
{
    switch(bits) {
    case PointGroups::C1:
        return "c1";
    case PointGroups::Ci:
        return "ci";
    case PointGroups::C2X:
    case PointGroups::C2Y:
    case PointGroups::C2Z:
        return "c2";
    case PointGroups::CsZ:
    case PointGroups::CsY:
    case PointGroups::CsX:
        return "cs";
    case PointGroups::D2:
        return "d2";
    case PointGroups::C2vX:
    case PointGroups::C2vY:
    case PointGroups::C2vZ:
        return "c2v";
    case PointGroups::C2hX:
    case PointGroups::C2hY:
    case PointGroups::C2hZ:
        return "c2h";
    case PointGroups::D2h:
        return "d2h";
    default:
        output::printf("Unrecognized point group bits: %d\n", bits);
        throw RuntimeError("Unrecognized point group bits");
    }
}

}  // end namespace PointGroups
} // end namespace panache

