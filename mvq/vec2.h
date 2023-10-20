#ifndef MVQ_VEC2_H
#define MVQ_VEC2_H

#include <cmath>
#include <cstdio>
#include "mathutil.h"
#include "vec3.h"

#define VEC2ELEMS(v) (v).x(), (v).y()

namespace mvq {

class vec2d;

//
// vec2f
//

class vec2f
{
public:
    vec2f()                                     { _v[0] = 0.0; _v[1] = 0.0; }
    vec2f(const vec2f& v)                       { _v[0] = v.x(); _v[1] = v.y(); }
    vec2f(float xx, float yy)                   { _v[0] = xx; _v[1] = yy; }
    vec2f(const float v[2])                     { _v[0] = v[0]; _v[1] = v[1]; }

    // Implicit cast to vec2d
    inline operator vec2d() const;

    void    set(float xx, float yy)             { _v[0] = xx;  _v[1] = yy; }
    void    get(float& xx, float& yy) const     { xx = _v[0];  yy = _v[1]; }

    // Negation
    vec2f   operator-() const                   { return vec2f(-_v[0], -_v[1]); }
    void    negate()                            { _v[0] = -_v[0]; _v[1] = -_v[1]; }


    // Addition, subtraction
    vec2f   operator+(const vec2f& v) const     { return vec2f(_v[0]+v.x(), _v[1]+v.y()); }
    void    operator+=(const vec2f& v)          { _v[0] += v.x(); _v[1] += v.y(); }

    vec2f   operator-(const vec2f& v) const     { return vec2f(_v[0]-v.x(), _v[1]-v.y()); }
    void    operator-=(const vec2f& v)          { _v[0] -= v.x(); _v[1] -= v.y(); }

    vec2f   operator*(float f) const            { return vec2f(f*_v[0], f*_v[1]); }
    void    operator*=(float f)                 { _v[0] *= f;  _v[1] *= f; }

    vec2f   operator/(float f) const
    {
        float inv_f = 1.0 / f;
        return vec2f(_v[0]*inv_f, _v[1]*inv_f);
    }
    void    operator/=(float f)
    {
        float inv_f = 1.0 / f;
        _v[0] *= inv_f;
        _v[1] *= inv_f;
    }

    // Dot product
    float   operator*(const vec2f& v)           { return _v[0]*v.x() + _v[1]*v.y(); }
    // Cross product: returns 3-vector!
    vec3f   operator^(const vec2f& v)           { return vec3f(0.0, 0.0, _v[0]*v.y() - _v[1]*v.x()); }

    // Length, squared length
    float   length() const                      { return sqrtf(_v[0]*_v[0] + _v[1]*_v[1]); }
    float   length2() const                     { return _v[0]*_v[0] + _v[1]*_v[1]; }

    // Returns length before normalization
    float   normalize()                         { float invl=1.0/length(); _v[0] *= invl; _v[1] *= invl; return 1.0/invl; }
    vec2f   normalized() const                  { float invl=1.0/length(); return vec2f(_v[0]*invl, _v[1]*invl); }

    vec2f   reciprocal() const                  { return vec2f(1.0/_v[0], 1.0/_v[1]); }

    float&  operator[](int i)                   { return _v[i]; }
    float   operator[](int i) const             { return _v[i]; }

    float&  x()                                 { return _v[0]; }
    float&  y()                                 { return _v[1]; }

    float   x() const                           { return _v[0]; }
    float   y() const                           { return _v[1]; }

    void print() const                          { printf("<%.6f %.6f>", _v[0], _v[1]); }

public:
    float _v[2];
};

// Scalar scaling

inline vec2f
operator*(float f, const vec2f& v)
{
    return vec2f(f*v.x(), f*v.y());
}

inline vec2f
operator*(double d, const vec2f& v)
{
    return vec2f(d*v.x(), d*v.y());
}

//
// vec2d
//

class vec2d
{
public:
    vec2d()                                     { _v[0] = 0.0; _v[1] = 0.0; }
    vec2d(const vec2d& v)                       { _v[0] = v.x(); _v[1] = v.y(); }
    vec2d(double xx, double yy)                 { _v[0] = xx; _v[1] = yy; }
    vec2d(const double v[2])                    { _v[0] = v[0]; _v[1] = v[1]; }

    // Implicit cast to vec2f
    inline operator vec2f() const               { return vec2f(_v[0], _v[1]); }

    void    set(double xx, double yy)           { _v[0] = xx;  _v[1] = yy; }
    void    get(double& xx, double& yy) const   { xx = _v[0];  yy = _v[1]; }

    // Negation
    vec2d   operator-() const                   { return vec2d(-_v[0], -_v[1]); }
    void    negate()                            { _v[0] = -_v[0]; _v[1] = -_v[1]; }

    // Addition, subtraction
    vec2d   operator+(const vec2d& v) const     { return vec2d(_v[0]+v.x(), _v[1]+v.y()); }
    void    operator+=(const vec2d& v)          { _v[0] += v.x(); _v[1] += v.y(); }

    vec2d   operator-(const vec2d& v) const     { return vec2d(_v[0]-v.x(), _v[1]-v.y()); }
    void    operator-=(const vec2d& v)          { _v[0] -= v.x(); _v[1] -= v.y(); }

    vec2d   operator*(double f) const           { return vec2d(f*_v[0], f*_v[1]); }
    void    operator*=(double f)                { _v[0] *= f;  _v[1] *= f; }

    vec2d   operator/(double f) const
    {
        double inv_f = 1.0 / f;
        return vec2d(_v[0]*inv_f, _v[1]*inv_f);
    }
    void    operator/=(double f)
    {
        double inv_f = 1.0 / f;
        _v[0] *= inv_f;
        _v[1] *= inv_f;
    }

    // Dot product
    double   operator*(const vec2d& v)          { return _v[0]*v.x() + _v[1]*v.y(); }
    // Cross product: returns 3-vector!
    vec3f   operator^(const vec2d& v)           { return vec3f(0.0, 0.0, _v[0]*v.y() - _v[1]*v.x()); }

    // Length, squared length
    double  length() const                      { return sqrt(_v[0]*_v[0] + _v[1]*_v[1]); }
    double  length2() const                     { return _v[0]*_v[0] + _v[1]*_v[1]; }

    // Returns length before normalization
    double   normalize()                        { double invl=1.0/length(); _v[0] *= invl; _v[1] *= invl; return 1.0/invl; }
    vec2d   normalized() const                  { double invl=1.0/length(); return vec2d(_v[0]*invl, _v[1]*invl); }

    vec2d   reciprocal() const                  { return vec2d(1.0/_v[0], 1.0/_v[1]); }

    double&  operator[](int i)                  { return _v[i]; }
    double   operator[](int i) const            { return _v[i]; }

    double&  x()                                { return _v[0]; }
    double&  y()                                { return _v[1]; }

    double   x() const                          { return _v[0]; }
    double   y() const                          { return _v[1]; }

    void print() const                          { printf("<%.12lf %.12lf>", _v[0], _v[1]); }

public:
    double _v[2];
};

// Scalar scaling

inline vec2d
operator*(double d, const vec2d& v)
{
    return vec2d(d*v.x(), d*v.y());
}

inline vec2d
operator*(float f, const vec2d& v)
{
    return vec2d(f*v.x(), f*v.y());
}

// Implicit cast

inline
vec2f::operator vec2d() const
{
    return vec2d(_v[0], _v[1]);
}

//
// Set default vec2 type
//

typedef vec2f   vec2;


} // namespace mvq

#endif
