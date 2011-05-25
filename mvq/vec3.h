#ifndef MVQ_VEC3_H
#define MVQ_VEC3_H

#include <cmath>
#include <cstdio>

#define VEC3ELEMS(v) (v).x(), (v).y(), (v).z()

namespace mvq {

class vec3d;

//
// vec3f
//

class vec3f
{
public:
    vec3f()                                     { _v[0] = _v[1] = _v[2] = 0.0; }
    vec3f(const vec3f& v)                       { _v[0] = v.x(); _v[1] = v.y(); _v[2] = v.z(); }
    vec3f(float xx, float yy, float zz)         { _v[0] = xx;  _v[1] = yy;  _v[2] = zz; }
    vec3f(const float v[3])                     { _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2]; }

    // Implicit cast to vec3d
    inline operator vec3d() const;

    void    set(float xx, float yy, float zz)   { _v[0] = xx;  _v[1] = yy;  _v[2] = zz; }
    void    get(float& xx, float& yy, float& zz){ xx = _v[0];  yy =_v[1];  zz = _v[2]; }

    // Negation
    vec3f   operator-() const                   { return vec3f(-_v[0], -_v[1], -_v[2]); }
    void    negate()                            { _v[0] = -_v[0]; _v[1] = -_v[1]; _v[2] = -_v[2]; }

    // Addition, subtraction
    vec3f   operator+(const vec3f& v) const     { return vec3f(_v[0]+v.x(), _v[1]+v.y(), _v[2]+v.z()); }
    void    operator+=(const vec3f& v)          { _v[0] += v.x(); _v[1] += v.y(); _v[2] += v.z(); }

    vec3f   operator-(const vec3f& v) const     { return vec3f(_v[0]-v.x(), _v[1]-v.y(), _v[2]-v.z()); }
    void    operator-=(const vec3f& v)          { _v[0] -= v.x(); _v[1] -= v.y(); _v[2] -= v.z(); }

    vec3f   operator*(float f) const            { return vec3f(f*_v[0], f*_v[1], f*_v[2]); }
    void    operator*=(float f)                 { _v[0] *= f;  _v[1] *= f;  _v[2] *= f; }

    vec3f   operator/(float f) const
    {
        register float inv_f = 1.0 / f;
        return vec3f(_v[0]*inv_f, _v[1]*inv_f, _v[2]*inv_f);
    }
    void    operator/=(float f)
    {
        register float inv_f = 1.0 / f;
        _v[0] *= inv_f;
        _v[1] *= inv_f;
        _v[2] *= inv_f;
    }

    // Dot product
    float   operator*(const vec3f& v) const     { return _v[0]*v.x() + _v[1]*v.y() + _v[2]*v.z(); }
    // Cross product
    vec3f   operator^(const vec3f& w) const
    {
        return vec3f(_v[1]*w.z() - _v[2]*w.y(), _v[2]*w.x() - _v[0]*w.z(), _v[0]*w.y()-_v[1]*w.x());
    }

    // Length, squared length
    float   length() const                      { return sqrtf(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]); }
    float   length2() const                     { return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]; }

    // Returns length before normalization
    float   normalize()
    {
        float invlen = 1.0 / length();

        _v[0] *= invlen;
        _v[1] *= invlen;
        _v[2] *= invlen;

        // XXX rounding?
        return (1.0 / invlen);
    }

    vec3f   normalized() const
    {
        float invlen = 1.0 / length();
        return vec3f(_v[0]*invlen, _v[1]*invlen, _v[2]*invlen);
    }

    vec3f   reciprocal() const                  { return vec3f(1.0/_v[0], 1.0/_v[1], 1.0/_v[2]); }

    float&  operator[](int i)                   { return _v[i]; }
    float   operator[](int i) const             { return _v[i]; }

    float&  x()                                 { return _v[0]; }
    float&  y()                                 { return _v[1]; }
    float&  z()                                 { return _v[2]; }

    float   x() const                           { return _v[0]; }
    float   y() const                           { return _v[1]; }
    float   z() const                           { return _v[2]; }

    void print() const                          { printf("<%.6f %.6f %.6f>", _v[0], _v[1], _v[2]); }

protected:
    float   _v[3];
};

// Scalar scaling

inline vec3f
operator*(float f, const vec3f& v)
{
    return vec3f(f*v.x(), f*v.y(), f*v.z());
}

inline vec3f
operator*(double d, const vec3f& v)
{
    return vec3f(d*v.x(), d*v.y(), d*v.z());
}

//
// vec3d
//

class vec3d
{
public:
    vec3d()                                         { _v[0] = _v[1] = _v[2] = 0.0; }
    vec3d(const vec3d& v)                           { _v[0] = v.x(); _v[1] = v.y(); _v[2] = v.z(); }
    vec3d(double xx, double yy, double zz)          { _v[0] = xx;  _v[1] = yy;  _v[2] = zz; }
    vec3d(const double v[3])                        { _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2]; }

    // Implicit cast to vec3f
    inline operator vec3f() const                   { return vec3f(_v[0], _v[1], _v[2]); }

    void    set(double xx, double yy, double zz)    { _v[0] = xx;  _v[1] = yy;  _v[2] = zz; }
    void    get(double& xx, double& yy, double& zz) { xx = _v[0];  yy =_v[1];  zz = _v[2]; }

    // Negation
    vec3d   operator-() const                       { return vec3d(-_v[0], -_v[1], -_v[2]); }
    void    negate()                                { _v[0] = -_v[0]; _v[1] = -_v[1]; _v[2] = -_v[2]; }

    // Addition, subtraction
    vec3d   operator+(const vec3d& v) const         { return vec3d(_v[0]+v.x(), _v[1]+v.y(), _v[2]+v.z()); }
    void    operator+=(const vec3d& v)              { _v[0] += v.x(); _v[1] += v.y(); _v[2] += v.z(); }

    vec3d   operator-(const vec3d& v) const         { return vec3d(_v[0]-v.x(), _v[1]-v.y(), _v[2]-v.z()); }
    void    operator-=(const vec3d& v)              { _v[0] -= v.x(); _v[1] -= v.y(); _v[2] -= v.z(); }

    vec3d   operator*(double f) const               { return vec3d(f*_v[0], f*_v[1], f*_v[2]); }
    void    operator*=(double f)                    { _v[0] *= f;  _v[1] *= f;  _v[2] *= f; }

    vec3d   operator/(double f) const
    {
        register double inv_f = 1.0 / f;
        return vec3d(_v[0]*inv_f, _v[1]*inv_f, _v[2]*inv_f);
    }
    void    operator/=(double f)
    {
        register double inv_f = 1.0 / f;
        _v[0] *= inv_f;
        _v[1] *= inv_f;
        _v[2] *= inv_f;
    }

    // Dot product
    double  operator*(const vec3d& v) const         { return _v[0]*v.x() + _v[1]*v.y() + _v[2]*v.z(); }
    // Cross product
    vec3d   operator^(const vec3d& w) const
    {
        return vec3d(_v[1]*w.z() - _v[2]*w.y(), _v[2]*w.x() - _v[0]*w.z(), _v[0]*w.y()-_v[1]*w.x());
    }

    // Length, squared length
    double  length() const                          { return sqrtf(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]); }
    double  length2() const                         { return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]; }

    // Returns length before normalization
    double  normalize()
    {
        double invlen = 1.0 / length();

        _v[0] *= invlen;
        _v[1] *= invlen;
        _v[2] *= invlen;

        // XXX rounding?
        return (1.0 / invlen);
    }

    vec3d   normalized() const
    {
        double invlen = 1.0 / length();
        return vec3d(_v[0]*invlen, _v[1]*invlen, _v[2]*invlen);
    }

    vec3d   reciprocal() const                      { return vec3d(1.0/_v[0], 1.0/_v[1], 1.0/_v[2]); }

    double& operator[](int i)                       { return _v[i]; }
    double  operator[](int i) const                 { return _v[i]; }

    double& x()                                     { return _v[0]; }
    double& y()                                     { return _v[1]; }
    double& z()                                     { return _v[2]; }

    double  x() const                               { return _v[0]; }
    double  y() const                               { return _v[1]; }
    double  z() const                               { return _v[2]; }

    void print() const                              { printf("<%.12f %.12f %.12f>", _v[0], _v[1], _v[2]); }

protected:
    double   _v[3];
};

// Scalar scaling

inline vec3d
operator*(double d, const vec3d& v)
{
    return vec3d(d*v.x(), d*v.y(), d*v.z());
}

inline vec3d
operator*(float f, const vec3d& v)
{
    return vec3d(f*v.x(), f*v.y(), f*v.z());
}

// Implicit cast

inline
vec3f::operator vec3d() const
{
    return vec3d(_v[0], _v[1], _v[2]);
}

//
// Set default vec3 type
//

typedef vec3f   vec3;

} // namespace mvq

#endif
