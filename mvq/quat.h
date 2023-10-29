#ifndef QUAT_H
#define QUAT_H

#include "mvq/vec3.h"
#include "mvq/mat4.h"

namespace mvq {

// A quaternion, used to represent rotations around arbitrary angles

class quatf
{
public:

    static quatf unit()
    {
        return quatf(1.0f, 0.0f, 0.0f, 0.0f);
    }

    // Create a quaternion that represents a rotation of
    // angle degrees around the given axis.
    static quatf rotation(vec3 axis, float angle)
    {
        const float l2 = axis.length2();
        if (l2 > 1e-6)
            axis = axis / (l2*l2);

        const float half_angle_rad = deg2rad(angle/2);
        const float sin_a = sin(half_angle_rad);

        return quatf(
            cos(half_angle_rad),
            axis.x() * sin_a,
            axis.y() * sin_a,
            axis.z() * sin_a
        ).normalized();
    }

    // Create a quaternion representation of point p
    // XXX implementation missing?
    static quatf point(vec3 p);

public:

    quatf()                                     { _v[0] = _v[1] = _v[2] = _v[3] = 0.0f; }
    quatf(const quatf& q)                       { _v[0] = q.w(); _v[1] = q.x(); _v[2] = q.y(); _v[3] = q.z(); }
    quatf(vec3f p)                              { _v[0] = 0.0f;  _v[1] = p.x();  _v[2] = p.y();  _v[3] = p.z(); }
    quatf(float ww, float xx, float yy, float zz) { _v[0] = ww;  _v[1] = xx;  _v[2] = yy;  _v[3] = zz; }
    quatf(const float q[4])                     { _v[0] = q[0];  _v[1] = q[1];  _v[2] = q[2];  _v[3] = q[3]; }

    void    set(float ww, float xx, float yy, float zz)     { _v[0] = ww; _v[1] = xx;  _v[2] = yy;  _v[3] = zz; }
    void    get(float& ww, float& xx, float& yy, float& zz) { ww = _v[0]; xx = _v[1];  yy =_v[2];  zz = _v[3]; }

    // Negation
    quatf   operator-() const                   { return quatf(-_v[0], -_v[1], -_v[2], -_v[3]); }
    void    negate()                            { _v[0] = -_v[0]; _v[1] = -_v[1]; _v[2] = -_v[2]; _v[3] = -_v[3]; }

    // Addition, subtraction
    quatf   operator+(const quatf& q) const     { return quatf(_v[0]+q.w(), _v[1]+q.x(), _v[2]+q.y(), _v[3]+q.z()); }
    void    operator+=(const quatf& q)          { _v[0] += q.w(); _v[1] += q.x(); _v[2] += q.y(); _v[3] += q.z(); }

    quatf   operator-(const quatf& q) const     { return quatf(_v[0]-q.w(), _v[1]-q.x(), _v[2]-q.y(), _v[3]-q.z()); }
    void    operator-=(const quatf& q)          { _v[0] -= q.w(); _v[1] -= q.x(); _v[2] -= q.y(); _v[3] = q.z(); }

    quatf   operator*(float f) const            { return quatf(f*_v[0], f*_v[1], f*_v[2], f*_v[3]); }
    void    operator*=(float f)                 { _v[0] *= f;  _v[1] *= f;  _v[2] *= f;  _v[3] *= f; }

    quatf   operator/(float f) const
    {
        float inv_f = 1.0 / f;
        return quatf(_v[0]*inv_f, _v[1]*inv_f, _v[2]*inv_f, _v[3]*inv_f);
    }
    void    operator/=(float f)
    {
        float inv_f = 1.0 / f;
        _v[0] *= inv_f;
        _v[1] *= inv_f;
        _v[2] *= inv_f;
        _v[3] *= inv_f;
    }

    // Length, squared length
    float   length() const                      { return sqrtf(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]); }
    float   length2() const                     { return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]; }

    // Returns length before normalization
    float   normalize()
    {
        float invlen = 1.0f / length();

        _v[0] *= invlen;
        _v[1] *= invlen;
        _v[2] *= invlen;
        _v[3] *= invlen;

        // XXX rounding?
        return (1.0f / invlen);
    }

    quatf   normalized() const
    {
        float invlen = 1.0f / length();
        return quatf(_v[0]*invlen, _v[1]*invlen, _v[2]*invlen, _v[3]*invlen);
    }

    // Calculate the conjugate (i.e. (w, -v) for quaternion (w, v))
    quatf   conjugate() const
    {
        return quatf(_v[0], -_v[1], -_v[2], -_v[3]);
    }

    // Calculate the inverse of a quaterion
    quatf   inverse() const
    {
        return conjugate() / length2();
    }


    // Transform a point (assumes quat represents a rotation, i.e. normalized)
    vec3f   transform(vec3 p) const;

    mat4d   to_rotation_matrix() const 
    {
        quatf qn = normalized();
        float rr  = qn.w();
        float rx = qn.x();
        float ry = qn.y();
        float rz = qn.z();
        return mat4d(
            1 - 2*(ry*ry + rz*rz),   2*(rx*ry - rr*rz),       2*(rx*rz + rr*ry),        0.0,
            2*(rx*ry + rr*rz),       1 - 2*(rx*rx + rz*rz),   2*(ry*rz - rr*rx),        0.0,
            2*(rx*rz - rr*ry),       2*(ry*rz + rr*rx),       1 - 2*(rx*rx + ry*ry),    0.0,
            0.0,                     0.0,                     0.0,                      1.0
        );
    }

public:
    float&  operator[](int i)                   { return _v[i]; }
    float   operator[](int i) const             { return _v[i]; }

    float&  w()                                 { return _v[0]; }
    float&  x()                                 { return _v[1]; }
    float&  y()                                 { return _v[2]; }
    float&  z()                                 { return _v[3]; }

    float   w() const                           { return _v[0]; }
    float   x() const                           { return _v[1]; }
    float   y() const                           { return _v[2]; }
    float   z() const                           { return _v[3]; }

    void    print() const                       { printf("<%.6f; %.6f %.6f %.6f>", _v[0], _v[1], _v[2], _v[3]); }

protected:
    float _v[4];        // w, x, y, z
};

// XXX not checked
inline quatf
operator*(const quatf& q, const quatf& r)
{
    return quatf(
        q.w()*r.w() - q.x()*r.x() - q.y()*r.y() - q.z()*r.z(),
        q.w()*r.x() + q.x()*r.w() + q.y()*r.z() - q.z()*r.y(),
        q.w()*r.y() + q.y()*r.w() + q.z()*r.x() - q.x()*r.z(),
        q.w()*r.z() + q.z()*r.w() + q.x()*r.y() - q.y()*r.x()
    );
}

inline vec3f
quatf::transform(vec3 p) const
{
    quatf r = *this * quatf(p) * conjugate();
    return vec3f(r.x(), r.y(), r.z());
}

//
// Set default quat type
//

typedef quatf   quat; 

} // mvq

#endif
