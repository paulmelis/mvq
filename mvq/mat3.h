#ifndef MVQ_MAT3F_H
#define MVQ_MAT3F_H

#include <cmath>
#include <cstdio>
#include <cstring>
#include "vec3.h"

namespace mvq {

//
// mat3f
//

class mat3f
{
public:
    static mat3f identity()
    {
        return mat3f(1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0);
    }

public:
    mat3f()
    {
        d[0][0] = d[1][0] = d[2][0] = 0.0;
        d[0][1] = d[1][1] = d[2][1] = 0.0;
        d[0][2] = d[1][2] = d[2][2] = 0.0;
    }

    mat3f(float m00, float m01, float m02,
        float m10, float m11, float m12,
        float m20, float m21, float m22)
    {
        d[0][0] = m00;
        d[0][1] = m01;
        d[0][2] = m02;
        d[1][0] = m10;
        d[1][1] = m11;
        d[1][2] = m12;
        d[2][0] = m20;
        d[2][1] = m21;
        d[2][2] = m22;
    }

    const char *__str__() const
    {
        static char    t[1024], s[128];
        static char     *p;
        t[0] = '\0';
        p = &t[0];
        for (int r = 0; r < 3; r++)
        {
            sprintf(s, "|%.6f, %.6f, %.6f|\n", d[r][0], d[r][1], d[r][2]);
            strcat(p, s);
            p += strlen(s);
        }
        return t;
    }

    const char *__repr__() const
    {
        return __str__();
    }

public:

    float operator()(int r, int c) const
    {
        return d[r][c];
    }

    float &operator()(int r, int c)
    {
        return d[r][c];
    }

    mat3f operator+(const mat3f& q) const
    {
        mat3f res;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                res.d[r][c] = d[r][c] + q.d[r][c];
        return res;
    }

    mat3f operator-(const mat3f& q) const
    {
        mat3f res;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                res.d[r][c] = d[r][c] - q.d[r][c];
        return res;
    }

    void operator+=(const mat3f& q)
    {
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                d[r][c] += q.d[r][c];
    }

    void operator-=(const mat3f& q)
    {
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                d[r][c] -= q.d[r][c];
    }

public:
    vec3f operator*(const vec3f& v) const
    {
        return vec3f(
            d[0][0]*v.x() + d[0][1]*v.y() + d[0][2]*v.z(),
            d[1][0]*v.x() + d[1][1]*v.y() + d[1][2]*v.z(),
            d[2][0]*v.x() + d[2][1]*v.y() + d[2][2]*v.z()
        );
    }

protected:
    float   d[3][3];        // row, column
};

/*
    column-vector times row-vector -> mat3f

        w.x  w.y  w.z
   v.x   *    .    .
   v.y   .    .
   v.z   .         .

*/

inline mat3f
matrix_product(const vec3f& v, const vec3f& w)
{
    return mat3f(
        v.x()*w.x(), v.x()*w.y(), v.x()*w.z(),
        v.y()*w.x(), v.y()*w.y(), v.y()*w.z(),
        v.z()*w.x(), v.z()*w.y(), v.z()*w.z()
    );
}

//
// mat3d
//

class mat3d
{
public:
    static mat3d identity()
    {
        return mat3d(1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0);
    }

public:
    mat3d()
    {
        d[0][0] = d[1][0] = d[2][0] = 0.0;
        d[0][1] = d[1][1] = d[2][1] = 0.0;
        d[0][2] = d[1][2] = d[2][2] = 0.0;
    }

    mat3d(double m00, double m01, double m02,
        double m10, double m11, double m12,
        double m20, double m21, double m22)
    {
        d[0][0] = m00;
        d[0][1] = m01;
        d[0][2] = m02;
        d[1][0] = m10;
        d[1][1] = m11;
        d[1][2] = m12;
        d[2][0] = m20;
        d[2][1] = m21;
        d[2][2] = m22;
    }

    const char *__str__() const
    {
        static char    t[1024], s[128];
        static char     *p;
        t[0] = '\0';
        p = &t[0];
        for (int r = 0; r < 3; r++)
        {
            sprintf(s, "|%.6f, %.6f, %.6f|\n", d[r][0], d[r][1], d[r][2]);
            strcat(p, s);
            p += strlen(s);
        }
        return t;
    }

    const char *__repr__() const
    {
        return __str__();
    }

public:

    double operator()(int r, int c) const
    {
        return d[r][c];
    }

    double &operator()(int r, int c)
    {
        return d[r][c];
    }

    mat3d operator+(const mat3d& q) const
    {
        mat3d res;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                res.d[r][c] = d[r][c] + q.d[r][c];
        return res;
    }

    mat3d operator-(const mat3d& q) const
    {
        mat3d res;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                res.d[r][c] = d[r][c] - q.d[r][c];
        return res;
    }

    void operator+=(const mat3d& q)
    {
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                d[r][c] += q.d[r][c];
    }

    void operator-=(const mat3d& q)
    {
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                d[r][c] -= q.d[r][c];
    }

public:
    vec3f operator*(const vec3f& v) const
    {
        return vec3f(
            d[0][0]*v.x() + d[0][1]*v.y() + d[0][2]*v.z(),
            d[1][0]*v.x() + d[1][1]*v.y() + d[1][2]*v.z(),
            d[2][0]*v.x() + d[2][1]*v.y() + d[2][2]*v.z()
        );
    }

protected:
    double   d[3][3];        // row, column
};

/*
    column-vector times row-vector -> mat3d

        w.x()  w.y()  w.z()
   v.x()   *    .    .
   v.y()   .    .
   v.z()   .         .

*/

inline mat3d
matrix_product(const vec3d& v, const vec3d& w)
{
    return mat3d(
        v.x()*w.x(), v.x()*w.y(), v.x()*w.z(),
        v.y()*w.x(), v.y()*w.y(), v.y()*w.z(),
        v.z()*w.x(), v.z()*w.y(), v.z()*w.z()
    );
}

//
// Set default mat3 type
//

typedef mat3d   mat3;

} // namespace mvq

#endif
