#include <cmath>
#include "panache/SimpleMatrix.h"
#include "panache/SphericalTransform.h"
#include "panache/SolidHarmonic.h"
#include "panache/BasisFunctionMacros.h" // INT_**

using panache::SimpleMatrix;

namespace panache {

void SphericalTransform::init()
{
//    output::printf("spher\n");
    int cartdim = INT_NCART(l_);
    SimpleMatrix coefmat(cartdim, cartdim);
    coefmat.zero();

    // Compute the solid harmonic matrix elements
    solidharmonic::solidharmonic(l_, coefmat);

//    output::printf("SphericalTransform: l = %d\n", l_);
//    coefmat.print();

    // Go through and grab the values.
    int pureindex = 0;
    int cartindex = 0;

    for (int i=1; i<=(l_-subl_)/2; ++i)
        pureindex += solidharmonic::npure(subl_+2*i);

    for (int p=0; p<solidharmonic::npure(subl_); ++p) {
        cartindex = 0;
//        for (int ii=0; ii<=l_; ++ii) {
//            int a = l_ - ii;
//            for (int jj=0; jj<=ii; ++jj) {
//                int b = ii - jj;
//                int c = jj;
        for (int a=0; a<=l_; ++a) {
            for (int b=0; (a+b)<=l_; ++b) {
                int c = l_ - a -b;

                int cart1 = solidharmonic::icart(a, b, c);
                int cart2 = INT_CARTINDEX(a+b+c, a, b);

//                output::printf("cart1 = %d, p+pureindex=%d\n", cart1, p+pureindex);
                double coef = coefmat(cart1, p+pureindex);

                if (fabs(coef) > 1.0e-16) {
                    components_.push_back(Component_(a, b, c, coef, cart2, p));
                }
                cartindex++;
            }
        }
    }
}

/*
void SphericalTransform::init_inverse()
{
//    output::printf("ispher\n");
    int cartdim = solidharmonic::ncart(l_);
    SimpleMatrix coefmat(cartdim, cartdim);
    coefmat.zero();

    // Compute the solid harmonic matrix elements
    solidharmonic::solidharmonic(l_, coefmat);

    // Invert and transpose the coefficient matrix
    coefmat.invert();
    coefmat.transpose_this();

    // Go through and grab the values.
    int pureindex = 0;
    int cartindex = 0;

    for (int i=1; i<=(l_-subl_)/2; ++i)
        pureindex += solidharmonic::npure(subl_+2*i);

    for (int p=0; p<solidharmonic::npure(subl_); ++p) {
        cartindex = 0;
//        for (int ii=0; ii<=l_; ++ii) {
//            int c = l_ - ii;
//            for (int jj=0; jj<=ii; ++jj) {
//                int b = ii - jj;
//                int a = jj;
        for (int a=0; a<=l_; ++a) {
            for (int b=0; (a+b)<=l_; ++b) {
                int c = l_ - a -b;

                int cart1 = solidharmonic::icart(a, b, c);
                int cart2 = INT_CARTINDEX(a+b+c, a, b);

                double coef = coefmat(cart1, p+pureindex);

                if (fabs(coef) > 1.0e-16) {
                    components_.push_back(Component_(a, b, c, coef, cart2, p));
                }
                cartindex++;
            }
        }
    }
}
*/

} // close namespace
