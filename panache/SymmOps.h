#ifndef PANACHE_SYMMOPS_H
#define PANACHE_SYMMOPS_H

#include <string>

namespace panache {

namespace SymmOps {
   enum Operations { E = 0, C2_z = 1, C2_y = 2, C2_x = 4, i = 8, Sigma_xy = 16, Sigma_xz = 32, Sigma_yz = 64, ID = 128 };
}

namespace PointGroups {

enum Groups {
    C1    = SymmOps::E,
    Ci    = SymmOps::E | SymmOps::i,
    C2X   = SymmOps::E | SymmOps::C2_x ,
    C2Y   = SymmOps::E | SymmOps::C2_y ,
    C2Z   = SymmOps::E | SymmOps::C2_z ,
    CsZ   = SymmOps::E | SymmOps::Sigma_xy ,
    CsY   = SymmOps::E | SymmOps::Sigma_xz ,
    CsX   = SymmOps::E | SymmOps::Sigma_yz ,
    D2    = SymmOps::E | SymmOps::C2_x | SymmOps::C2_y | SymmOps::C2_z ,
    C2vX  = SymmOps::E | SymmOps::C2_x | SymmOps::Sigma_xy | SymmOps::Sigma_xz ,
    C2vY  = SymmOps::E | SymmOps::C2_y | SymmOps::Sigma_xy | SymmOps::Sigma_yz ,
    C2vZ  = SymmOps::E | SymmOps::C2_z | SymmOps::Sigma_xz | SymmOps::Sigma_yz ,
    C2hX  = SymmOps::E | SymmOps::C2_x | SymmOps::Sigma_yz | SymmOps::i ,
    C2hY  = SymmOps::E | SymmOps::C2_y | SymmOps::Sigma_xz | SymmOps::i ,
    C2hZ  = SymmOps::E | SymmOps::C2_z | SymmOps::Sigma_xy | SymmOps::i ,
    D2h   = SymmOps::E | SymmOps::C2_x | SymmOps::C2_y | SymmOps::C2_z | SymmOps::i |
            SymmOps::Sigma_xy | SymmOps::Sigma_xz | SymmOps::Sigma_yz
};


void similar(unsigned char bits, unsigned char* sim, char& cnt);
bool full_name_to_bits(const std::string& pg, unsigned char &bits);
const char* bits_to_full_name(unsigned char bits);
const char* bits_to_basic_name(unsigned char bits);


} // end namespace PointGroups
} // end namespace panache

#endif //PANACHE_SYMMOPS_H
