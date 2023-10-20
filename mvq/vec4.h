#ifndef MVQ_VEC4_H
#define MVQ_VEC4_H

#include <cmath>
#include <cstdio>

#define VEC4ELEMS(v) (v).x(), (v).y(), (v).z(), (v).w()

namespace mvq {

class vec4d;

//
// vec4f
//

class vec4f
{
public:
    vec4f()                                         { _v[0] = _v[1] = _v[2] = _v[3] = 0.0; }
    vec4f(const vec4f& v)                           { _v[0] = v.x(); _v[1] = v.y(); _v[2] = v.z(); _v[3] = v.w(); }
    vec4f(float xx, float yy, float zz, float ww)   { _v[0] = xx;  _v[1] = yy;  _v[2] = zz;  _v[3] = ww; }
    vec4f(const float v[4])                         { _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2]; _v[3] = v[3]; }

    // Implicit cast to vec4d
    inline operator vec4d() const;

    void    set(float xx, float yy, float zz, float ww)       { _v[0] = xx;  _v[1] = yy;  _v[2] = zz;  _v[3] = ww; }
    void    get(float& xx, float& yy, float& zz, float ww)    { xx = _v[0];  yy =_v[1];  zz = _v[2];  ww = _v[3]; }

    // Negation
    vec4f   operator-() const                       { return vec4f(-_v[0], -_v[1], -_v[2], -_v[3]); }
    void    negate()                                { _v[0] = -_v[0]; _v[1] = -_v[1]; _v[2] = -_v[2]; _v[3] = -_v[3]; }

    // Addition, subtraction
    vec4f   operator+(const vec4f& v) const         { return vec4f(_v[0]+v.x(), _v[1]+v.y(), _v[2]+v.z(), _v[3]+v.w()); }
    void    operator+=(const vec4f& v)              { _v[0] += v.x(); _v[1] += v.y(); _v[2] += v.z(); _v[3] += v.w(); }

    vec4f   operator-(const vec4f& v) const         { return vec4f(_v[0]-v.x(), _v[1]-v.y(), _v[2]-v.z(), _v[3]-v.w()); }
    void    operator-=(const vec4f& v)              { _v[0] -= v.x(); _v[1] -= v.y(); _v[2] -= v.z(); _v[3] -= v.w(); }

    vec4f   operator*(float f) const                { return vec4f(f*_v[0], f*_v[1], f*_v[2], f*_v[3]); }
    void    operator*=(float f)                     { _v[0] *= f;  _v[1] *= f;  _v[2] *= f;  _v[3] *= f; }

    vec4f   operator/(float f) const
    {
        float inv_f = 1.0 / f;
        return vec4f(_v[0]*inv_f, _v[1]*inv_f, _v[2]*inv_f, _v[3]*inv_f);
    }
    void    operator/=(float f)
    {
        float inv_f = 1.0 / f;
        _v[0] *= inv_f;
        _v[1] *= inv_f;
        _v[2] *= inv_f;
        _v[3] *= inv_f;
    }

    // Dot product
    float   operator*(const vec4f& v) const     { return _v[0]*v.x() + _v[1]*v.y() + _v[2]*v.z() + _v[3]*v.w(); }

    // Length, squared length
    float   length() const                      { return sqrtf(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]); }
    float   length2() const                     { return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]; }

    // Returns length before normalization
    float   normalize()
    {
        float invlen = 1.0 / length();

        _v[0] *= invlen;
        _v[1] *= invlen;
        _v[2] *= invlen;
        _v[3] *= invlen;

        // XXX rounding?
        return (1.0 / invlen);
    }

    vec4f   normalized() const
    {
        float invlen = 1.0 / length();
        return vec4f(_v[0]*invlen, _v[1]*invlen, _v[2]*invlen, _v[3]*invlen);
    }

    vec4f   reciprocal() const                  { return vec4f(1.0/_v[0], 1.0/_v[1], 1.0/_v[2], 1.0/_v[3]); }

    float&  operator[](int i)                   { return _v[i]; }
    float   operator[](int i) const             { return _v[i]; }

    float&  x()                                 { return _v[0]; }
    float&  y()                                 { return _v[1]; }
    float&  z()                                 { return _v[2]; }
    float&  w()                                 { return _v[3]; }

    float   x() const                           { return _v[0]; }
    float   y() const                           { return _v[1]; }
    float   z() const                           { return _v[2]; }
    float   w() const                           { return _v[3]; }

    void print() const                          { printf("<%.6f %.6f %.6f, %.6f>", _v[0], _v[1], _v[2], _v[3]); }

protected:
    float   _v[4];
};

// Scalar scaling

inline vec4f
operator*(float f, const vec4f& v)
{
    return vec4f(f*v.x(), f*v.y(), f*v.z(), f*v.w());
}

inline vec4f
operator*(double d, const vec4f& v)
{
    return vec4f(d*v.x(), d*v.y(), d*v.z(), d*v.w());
}

//
// vec4d
//

class vec4d
{
public:
    vec4d()                                             { _v[0] = _v[1] = _v[2] = _v[3] = 0.0; }
    vec4d(const vec4d& v)                               { _v[0] = v.x(); _v[1] = v.y(); _v[2] = v.z(); _v[3] = v.w(); }
    vec4d(double xx, double yy, double zz, double ww)   { _v[0] = xx;  _v[1] = yy;  _v[2] = zz;  _v[3] = ww; }
    vec4d(const double v[4])                            { _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2]; _v[3] = v[3]; }

    // Implicit cast to vec4f
    inline operator vec4f() const                       { return vec4f(_v[0], _v[1], _v[2], _v[3]); }

    void    set(double xx, double yy, double zz, float ww)        { _v[0] = xx;  _v[1] = yy;  _v[2] = zz;  _v[3] = ww; }
    void    get(double& xx, double& yy, double& zz, float ww)     { xx = _v[0];  yy =_v[1];  zz = _v[2];  ww = _v[3]; }

    // Negation
    vec4d   operator-() const                           { return vec4d(-_v[0], -_v[1], -_v[2], -_v[3]); }
    void    negate()                                    { _v[0] = -_v[0]; _v[1] = -_v[1]; _v[2] = -_v[2]; _v[3] = -_v[3]; }

    // Addition, subtraction
    vec4d   operator+(const vec4d& v) const             { return vec4d(_v[0]+v.x(), _v[1]+v.y(), _v[2]+v.z(), _v[3]+v.w()); }
    void    operator+=(const vec4d& v)                  { _v[0] += v.x(); _v[1] += v.y(); _v[2] += v.z(); _v[3] += v.w(); }

    vec4d   operator-(const vec4d& v) const             { return vec4d(_v[0]-v.x(), _v[1]-v.y(), _v[2]-v.z(), _v[3]-v.w()); }
    void    operator-=(const vec4d& v)                  { _v[0] -= v.x(); _v[1] -= v.y(); _v[2] -= v.z(); _v[3] -= v.w(); }

    vec4d   operator*(double f) const                   { return vec4d(f*_v[0], f*_v[1], f*_v[2], f*_v[3]); }
    void    operator*=(double f)                        { _v[0] *= f;  _v[1] *= f;  _v[2] *= f;  _v[3] *= f; }

    vec4d   operator/(double f) const
    {
        double inv_f = 1.0 / f;
        return vec4d(_v[0]*inv_f, _v[1]*inv_f, _v[2]*inv_f, _v[3]*inv_f);
    }
    void    operator/=(double f)
    {
        double inv_f = 1.0 / f;
        _v[0] *= inv_f;
        _v[1] *= inv_f;
        _v[2] *= inv_f;
        _v[3] *= inv_f;
    }

    // Dot product
    double   operator*(const vec4d& v) const            { return _v[0]*v.x() + _v[1]*v.y() + _v[2]*v.z() + _v[3]*v.w(); }

    // Length, squared length
    double   length() const                             { return sqrtf(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]); }
    double   length2() const                            { return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]; }

    // Returns length before normalization
    double   normalize()
    {
        double invlen = 1.0 / length();

        _v[0] *= invlen;
        _v[1] *= invlen;
        _v[2] *= invlen;
        _v[3] *= invlen;

        // XXX rounding?
        return (1.0 / invlen);
    }

    vec4d   normalized() const
    {
        double invlen = 1.0 / length();
        return vec4d(_v[0]*invlen, _v[1]*invlen, _v[2]*invlen, _v[3]*invlen);
    }

    vec4d   reciprocal() const                          { return vec4d(1.0/_v[0], 1.0/_v[1], 1.0/_v[2], 1.0/_v[3]); }

    double&  operator[](int i)                          { return _v[i]; }
    double   operator[](int i) const                    { return _v[i]; }

    double&  x()                                        { return _v[0]; }
    double&  y()                                        { return _v[1]; }
    double&  z()                                        { return _v[2]; }
    double&  w()                                        { return _v[3]; }

    double   x() const                                  { return _v[0]; }
    double   y() const                                  { return _v[1]; }
    double   z() const                                  { return _v[2]; }
    double   w() const                                  { return _v[3]; }

    void print() const                                  { printf("<%.12f %.12f %.12f, %.12f>", _v[0], _v[1], _v[2], _v[3]); }

protected:
    double   _v[4];
};

// Scalar scaling

inline vec4d
operator*(double d, const vec4d& v)
{
    return vec4d(d*v.x(), d*v.y(), d*v.z(), d*v.w());
}

inline vec4d
operator*(float f, const vec4d& v)
{
    return vec4d(f*v.x(), f*v.y(), f*v.z(), f*v.w());
}

// Implicit cast

inline
vec4f::operator vec4d() const
{
    return vec4d(_v[0], _v[1], _v[2], _v[3]);
}

//
// Set default vec4 type
//

typedef vec4f   vec4;


} // namespace mvq

#endif
