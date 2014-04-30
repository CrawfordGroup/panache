#ifndef PANACHE_VECTOR3_H
#define PANACHE_VECTOR3_H

#include <cmath>
#include "Exception.h"

namespace panache {

// Forward declaration
class Vector3;
}


panache::Vector3 operator*(double d, const panache::Vector3& x);

namespace panache {

class Vector3
{
private:
    double v_[3];

public:
    Vector3() {
        v_[0] = v_[1] = v_[2] = 0.0;
    }

    Vector3(const double p[3]) {
        v_[0] = p[0];
        v_[1] = p[1];
        v_[2] = p[2];
    }

    Vector3(double d) {
        v_[0] = v_[1] = v_[2] = d;
    }

    Vector3(double x, double y, double z) {
        v_[0] = x;
        v_[1] = y;
        v_[2] = z;
    }

    Vector3(const Vector3& p) {
        v_[0] = p.v_[0];
        v_[1] = p.v_[1];
        v_[2] = p.v_[2];
    }

    Vector3 & operator=(const double *x) {
        v_[0] = x[0];
        v_[1] = x[1];
        v_[2] = x[2];
        return *this;
    }

    Vector3 & operator=(const Vector3& x) {
        v_[0] = x.v_[0];
        v_[1] = x.v_[1];
        v_[2] = x.v_[2];
        return *this;
    }

    Vector3 & operator=(double d) {
        v_[0] = d;
        v_[1] = d;
        v_[2] = d;
        return *this;
    }

    Vector3 & operator+=(const Vector3& x) {
        v_[0] += x.v_[0];
        v_[1] += x.v_[1];
        v_[2] += x.v_[2];
        return *this;
    }

    Vector3 & operator-=(const Vector3& x) {
        v_[0] -= x.v_[0];
        v_[1] -= x.v_[1];
        v_[2] -= x.v_[2];
        return *this;
    }

    Vector3 & operator*=(const Vector3& x) {
        v_[0] *= x.v_[0];
        v_[1] *= x.v_[1];
        v_[2] *= x.v_[2];
        return *this;
    }

    Vector3 & operator*=(double m) {
        v_[0] *= m;
        v_[1] *= m;
        v_[2] *= m;
        return *this;
    }

    Vector3 & operator/=(double m) {
        v_[0] /= m;
        v_[1] /= m;
        v_[2] /= m;
        return *this;
    }


    Vector3 operator+(const Vector3& x) {
        Vector3 result(*this);
        result += x; 
        return result;
    }

    Vector3 operator-(const Vector3& x) const {
        Vector3 result(*this);
        result -= x;
        return result;
    }

    Vector3 operator*(double d) const
    {
        Vector3 result(*this);
        result *= d;
        return result;
    }

    Vector3 operator*(const Vector3& x) const
    {
        Vector3 result(*this);
        result *= x;
        return result;
    }

    Vector3 operator/(double d) const
    {
        Vector3 result(*this);
        result /= d;
        return result;
    }


    Vector3 operator-() const {
        return Vector3(-v_[0], -v_[1], -v_[2]);
    }

    double& operator[](int i) {
        return v_[i];
    }

    const double& operator[](int i) const {
        return v_[i];
    }

    double get(int i) {
        if (i >= 0 && i <= 2)
            return v_[i];
        else
            throw RuntimeError("Bad vector index");
    }

    double dot(const Vector3& x) const {
        return v_[0]*x.v_[0] + v_[1]*x.v_[1] + v_[2]*x.v_[2];
    }

    double distance(const Vector3& x) const
    {
        Vector3 diff = (*this)-x;
        return diff.norm();
    }

    void normalize()
    {
        double n = 1.0/norm();
        (*this) *= n;
    }

    double norm() const {
        return sqrt(this->dot(*this));
    }

    Vector3 cross(const Vector3& x) const
    {
        Vector3 result(v_[1] * x.v_[2] - v_[2] * x.v_[1],
                       v_[2] * x.v_[0] - v_[0] * x.v_[2],
                       v_[0] * x.v_[1] - v_[1] * x.v_[0]);
        return result;
    }


    void rotate(double theta, Vector3& axis)
    {
        Vector3 result;
        Vector3 unitaxis = axis;
        unitaxis.normalize();

        // split into parallel and perpendicular components along axis
        Vector3 parallel = axis * (this->dot(axis) / axis.dot(axis));
        Vector3 perpendicular = (*this) - parallel;

        // form unit vector perpendicular to parallel and perpendicular
        Vector3 third_axis = axis.perp_unit(perpendicular);
        third_axis = third_axis * perpendicular.norm();

        result = parallel + cos(theta) * perpendicular + sin(theta) * third_axis;
        (*this) = result;
    }

    Vector3 perp_unit(const Vector3& v) const
    {
        // try cross product
        Vector3 result = cross(v);
        double resultdotresult = result.dot(result);

        if (resultdotresult < 1.e-16) {
            // cross product is too small to normalize
            // find the largest of this and v
            double dotprodt = this->dot(*this);
            double dotprodv = v.dot(v);
            const Vector3 *d;
            double dotprodd;
            if (dotprodt < dotprodv) {
                d = &v;
                dotprodd = dotprodv;
            }
            else {
                d = this;
                dotprodd = dotprodt;
            }

            // see if d is big enough
            if (dotprodd < 1.e-16) {
                // choose an arbitrary vector, since the biggest vector is small
                result[0] = 1.0;
                result[1] = 0.0;
                result[2] = 0.0;
                return result;
            }
            else {
                // choose a vector prependicular to d
                // choose it in one of the planes xy, xz, yz
                // choose the plane to be that which contains the two largest
                // components of d
                double absd[3];
                absd[0] = fabs(d->v_[0]);
                absd[1] = fabs(d->v_[1]);
                absd[2] = fabs(d->v_[2]);
                int axis0, axis1;
                if ((absd[1] - absd[0]) > 1.0e-12) {
                    axis0 = 1;
                    if ((absd[2] - absd[0]) > 1.0e-12) {
                        axis1 = 2;
                    }
                    else {
                        axis1 = 0;
                    }
                }
                else {
                    axis0 = 0;
                    if ((absd[2] - absd[1]) > 1.0e-12) {
                        axis1 = 2;
                    }
                    else {
                        axis1 = 1;
                    }
                }

                result[0] = 0.0;
                result[1] = 0.0;
                result[2] = 0.0;
                // do the pi/2 rotation in the plane
                result[axis0] = d->v_[axis1];
                result[axis1] = -d->v_[axis0];
            }
            result.normalize();
            return result;
        }
        else {
            // normalize the cross product and return the result
            result *= 1.0/sqrt(resultdotresult);
            return result;
        }
    }
/*
    std::string to_string() const {
        std::stringstream s;
        s << "[ " << v_[0] << ", " << v_[1] << ", " << v_[2] << " ]";
        return s.str();
    }
*/
};

} // end namespace panache

inline panache::Vector3 operator*(double d, const panache::Vector3& x)
{
    return x*d;
}

#endif //PANACHE_VECTOR3_H
