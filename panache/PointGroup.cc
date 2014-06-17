
#include <cstdlib>
#include <cctype>

#include "PointGroup.h"

using std::shared_ptr;

namespace panache {

////////////////////////////////////////////////////////////////////////

PointGroup::PointGroup()
{
    set_symbol("c1");
    origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const std::string& s)
{
    if (PointGroups::full_name_to_bits(s, bits_) == false)
        throw RuntimeError("PointGroup: Unknown point group name provided.");
    set_symbol(PointGroups::bits_to_basic_name(bits_));
    origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const std::string& s,
                       const Vector3& origin)
{
    if (PointGroups::full_name_to_bits(s, bits_) == false)
        throw RuntimeError("PointGroup: Unknown point group name provided.");
    set_symbol(PointGroups::bits_to_basic_name(bits_));
    origin_ = origin;
}

PointGroup::PointGroup(unsigned char bits)
    : bits_(bits)
{
    set_symbol(PointGroups::bits_to_basic_name(bits));
    origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(unsigned char bits, const Vector3 & origin)
    : bits_(bits)
{
    set_symbol(PointGroups::bits_to_basic_name(bits));
    origin_ = origin;
}

PointGroup::PointGroup(const PointGroup& pg)
{
    *this = pg;
}

PointGroup::PointGroup(const shared_ptr<PointGroup>& pg)
{
    *this = *pg.get();
}

PointGroup::~PointGroup()
{
}

PointGroup&
PointGroup::operator=(const PointGroup& pg)
{
    set_symbol(pg.symb);
    origin_ = pg.origin_;
    return *this;
}

void
PointGroup::set_symbol(const std::string& sym)
{
    if (sym.length()) {
        symb = sym;
    } else {
        set_symbol("c1");
    }
}

CharacterTable
PointGroup::char_table() const
{
    CharacterTable ret(bits_);
    return ret;
}

int
PointGroup::equiv(const shared_ptr<PointGroup> &grp, double tol) const
{
    if (symb != grp->symb)
        return 0;

    return 1;
}


//void
//PointGroup::print(FILE *out) const
//{
//    output::printf("PointGroup: %s\n", symb.c_str());
//}

} // end namespace panache
