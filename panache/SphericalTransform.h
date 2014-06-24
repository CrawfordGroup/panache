#ifndef PANACHE_SPHERICALTRANSFORM_H
#define PANACHE_SPHERICALTRANSFORM_H

#include <vector>
#include <array>

namespace panache {


/** This is a base class for a container for a sparse Cartesian to solid
    harmonic basis function transformation. */
class SphericalTransform
{
private:

    struct Component_ {
        int a, b, c;
        double coef;
        int cartindex, pureindex;

        Component_(int a, int b, int c, double coef, int cartindex, int pureindex)
           : a(a),b(b),c(c),coef(coef),cartindex(cartindex),pureindex(pureindex)
        { }
    };


    std::vector<Component_> components_;

    int l_; // The angular momentum this transform is for.
    int subl_;



    SphericalTransform(int l, int subl) : l_(l)
    {
        if (subl == -1)
            subl_ = l;
        else
            subl_ = subl;
    }


    void init();
    //void init_inverse();


public:
    typedef std::vector<Component_>::const_iterator const_component_iterator;

    /// Returns the Cartesian basis function index of component i
    int cartindex(int i) const { return components_[i].cartindex; }
    /// Returns the spherical harmonic basis index of component i
    int pureindex(int i) const { return components_[i].pureindex; }
    /// Returns the transformation coefficient of component i
    double coef(int i) const { return components_[i].coef; }
    /// Returns the Cartesian basis function's x exponent of component i
    int a(int i) const { return components_[i].a; }
    /// Returns the Cartesian basis function's y exponent of component i
    int b(int i) const { return components_[i].b; }
    /// Returns the Cartesian basis function's z exponent of component i
    int c(int i) const { return components_[i].c; }
    /// Returns the number of components in the transformation
    int n() const { return components_.size(); }
    /// Returns the angular momentum
    int l() const { return l_; }

    const_component_iterator cbegin() const { return components_.cbegin(); }
    const_component_iterator cend() const { return components_.cend(); }

    // Inefficiency due to multiple copying shouldn't be a problem with this class
    //   particularly because of return value optimization.
    //   If for whatever reason it is, copy-and-swap idiom can help
    static SphericalTransform Generate(int l, int subl = -1)
    {
        SphericalTransform r(l, subl);
        r.init();
        return r;
    }

    /*
    static SphericalTransform GenerateInverse(int l, int subl = -1)
    {
        SphericalTransform r(l, subl);
        r.init_inverse();
        return r;
    }
    */

};



} // close namespace panache

#endif //PANACHE_SPHERICALTRANSFORM_H
