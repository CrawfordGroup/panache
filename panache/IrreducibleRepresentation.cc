#include <cstdlib>
#include "IrreducibleRepresentation.h"

namespace panache {

/////////////////////////////////////////////////////////////////////////

IrreducibleRepresentation::IrreducibleRepresentation() :
    g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(0), rep(0), csymb(0)
{
}

IrreducibleRepresentation::IrreducibleRepresentation(
    int order, int d, const char *lab, const char *clab) :
    g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(0), rep(0), csymb(0)
{
    init(order,d,lab,clab);
}


IrreducibleRepresentation::IrreducibleRepresentation(
    const IrreducibleRepresentation& ir) :
    g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(0), rep(0), csymb(0)
{
    *this = ir;
}

IrreducibleRepresentation::~IrreducibleRepresentation()
{
    init();
}

IrreducibleRepresentation&
IrreducibleRepresentation::operator=(const IrreducibleRepresentation& ir)
{
    init(ir.g,ir.degen,ir.symb,ir.csymb);

    nrot_ = ir.nrot_;
    ntrans_ = ir.ntrans_;
    complex_ = ir.complex_;

    for (int i=0; i < g; i++)
        rep[i]=ir.rep[i];

    return *this;
}

void
IrreducibleRepresentation::init(int order, int d, const char *lab,
                                const char *clab)
{
    g=order;
    degen=d;
    ntrans_=nrot_=complex_=0;

    free(symb);
    if (lab)
        symb = strdup(lab);
    else
        symb = NULL;

    free(csymb);
    if (clab) csymb = strdup(clab);
    else csymb = 0;

    if (rep) {
        delete[] rep;
        rep=0;
    }

    if (g) {
        rep = new SymRep[g];
        for (int i=0; i < g; i++)
            rep[i].set_dim(d);
    }
}

}

