#ifndef MVQ_MAT4_H
#define MVQ_MAT4_H

#include "vec3.h"
#include "vec4.h"
#include "mathutil.h"

#include <cstdio>

namespace mvq {

//
// Note: assumes pre-multiplication order, i.e. v*M!
//

class mat4d
{
public:

    //
    // static members
    //

    inline static mat4d identity();
    inline static mat4d diagonal(double e00, double e11, double e22, double e33);

    inline static mat4d translate(double x, double y, double z);
    inline static mat4d translate(vec3 t);
    inline static mat4d scale(double x, double y, double z);
    inline static mat4d scale(vec3 s);

    inline static mat4d rotate_x(double angle);
    inline static mat4d	rotate_y(double angle);
    inline static mat4d	rotate_z(double angle);

    inline static mat4d perspective(double fovy, double aspect, double z_near, double z_far);
    inline static mat4d look_at(vec3d eye, vec3d center, vec3d up);
// XXX untested as of yet!
    inline static mat4d frustum(double left, double right, double bottom, double top, double near, double far);

public:

    //
    // constructors
    //

    // zero
    mat4d();

    // homogenous
    mat4d(double h);

    // anything
    mat4d(double e00, double e01, double e02, double e03,
          double e10, double e11, double e12, double e13,
          double e20, double e21, double e22, double e23,
          double e30, double e31, double e32, double e33);

    mat4d(float *m);
    mat4d(double *m);

    // copy
    mat4d(const mat4d& m);

    //
    // instance members
    //

    // initialization
    inline void     set(float m[]);
    inline void     set(double m[]);

    // setting basic values

    inline void     set_zero();
    inline void     set_identity();
    inline void     set_diagonal(double e00, double e11, double e22, double e33);

    // element access

    double&         operator()(int row, int col)        { return d[row][col]; }
    double          operator()(int row, int col) const  { return d[row][col]; }

    void            set(int row, int col, double value) { d[row][col] = value; }
    double          get(int row, int col) const         { return d[row][col]; }

    const double    *ptr() const                        { return (double*)d; }

    // operations

    inline mat4d    operator-() const;
    inline void     negate();

    friend mat4d    operator*(const mat4d& m, double d);
    friend mat4d    operator*(double d, const mat4d& m);
    friend mat4d    operator/(const mat4d& m, double d);

    inline void     operator*=(double s);
    inline void     operator/=(double s);

    friend mat4d    operator+(const mat4d& m1, const mat4d& m2);
    friend mat4d    operator-(const mat4d& m1, const mat4d& m2);

    friend vec4d    operator*(const vec4d v, const mat4d& m);

    friend mat4d    operator*(const mat4d& m1, const mat4d& m2);

    inline mat4d    transpose() const;
    inline mat4d    inverse() const;

    // vec3 transforms

    // assumes p is a point (v.w == 1)
    inline vec3f    transform(float x, float y, float z) const;
    inline vec3d    transform(double x, double y, double z) const;
    inline vec3f    transform(const vec3f p) const;
    inline vec3d    transform(const vec3d p) const;
    // assumes v is a vector (v.w == 0)
    inline vec3f    vtransform(float x, float y, float z) const;
    inline vec3d    vtransform(double x, double y, double z) const;
    inline vec3f    vtransform(const vec3f v) const;
    inline vec3d    vtransform(const vec3d v) const;

    // For transforming normal vectors, i.e. the transform (M^-1)^T
    // (followed by using vtransform())
    inline  mat4d   normal_transform() const;

    // tests

    inline bool     is_identity() const;
    inline bool     is_zero() const;

    void dump() const
    {
        printf("\n");
        for (int r = 0; r < 4; r++)
            printf("|  %.3f %.3f %.3f %.3f |\n", d[r][0], d[r][1], d[r][2], d[r][3]);
        printf("\n");
    }

protected:
    double  d[4][4];
};

//
// static methods
//

inline mat4d
mat4d::identity()
{
    return mat4d(
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    );
}

inline mat4d
mat4d::diagonal(double e00, double e11, double e22, double e33)
{
    return mat4d(
        e00, 0.0, 0.0, 0.0,
        0.0, e11, 0.0, 0.0,
        0.0, 0.0, e22, 0.0,
        0.0, 0.0, 0.0, e33);
}

inline mat4d
mat4d::translate(double x, double y, double z)
{
    return mat4d(
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        x,   y,   z,   1.0);
}

inline mat4d
mat4d::translate(vec3 t)
{
    return mat4d(
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        t.x(), t.y(), t.z(), 1.0);
}

/*
inline mat4d
mat4d::translation(vec3d v)
{
    return mat4d(1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, v.x,v.y,v.z,1.0);
}
*/

inline mat4d
mat4d::scale(double x, double y, double z)
{
    return mat4d(
        x,   0.0, 0.0, 0.0,
        0.0, y,   0.0, 0.0,
        0.0, 0.0, z,   0.0,
        0.0, 0.0, 0.0, 1.0);
}

inline mat4d
mat4d::scale(vec3 s)
{
    return mat4d(
        s.x(),  0.0,    0.0,    0.0,
        0.0,    s.y(),  0.0,    0.0,
        0.0,    0.0,    s.z(),  0.0,
        0.0,    0.0,    0.0,    1.0); 
}

/*
inline mat4d
mat4d::scaling(vec3d v)
{
    return mat4d(v.x,0.0,0.0,0.0, 0.0,v.y,0.0,0.0, 0.0,0.0,v.z,0.0, 0.0,0.0,0.0,1.0);
}
*/

inline mat4d
mat4d::rotate_x(double angle)
{
    return mat4d(
        1.0, 	0.0,		  0.0,		 	0.0,
        0.0, 	cos(angle),   sin(angle),	0.0,
        0.0, 	-sin(angle),  cos(angle),	0.0,
        0.0, 	0.0,		  0.0,		 	1.0
    );
}

inline mat4d
mat4d::rotate_y(double angle)
{
    return mat4d(
        cos(angle),		0.0,	-sin(angle),	0.0,
        0.0,			1.0,	0.0,			0.0,
        sin(angle),		0.0,	cos(angle),		0.0,
        0.0,			0.0,	0.0,			1.0
    );
}

inline mat4d
mat4d::rotate_z(double angle)
{
    return mat4d(
        cos(angle),		sin(angle),		0.0,	0.0,
        -sin(angle),	cos(angle),		0.0,	0.0,
        0.0,			0.0,			1.0,	0.0,
        0.0,			0.0,			0.0,	1.0
    );
}

inline mat4d
mat4d::perspective(double fovy, double aspect, double z_near, double z_far)
{
    double f = atan(0.5*fovy);
    return mat4d(
        f/aspect,   0.0,    0.0,    0.0,
        0.0,        f,      0.0,    0.0,
        0.0,        0.0,    (z_far+z_near)/(z_near-z_far), -1.0,
        0.0,        0.0,    2.0*z_far*z_near/(z_near-z_far), 0.0
    );
}

inline mat4d
mat4d::look_at(vec3d eye, vec3d center, vec3d up)
{
    vec3d f(center-eye);
    f.normalize();
    vec3d s(f^up);
    s.normalize();
    vec3d u(s^f);
    u.normalize();

    return mat4d::translate(-eye.x(), -eye.y(), -eye.z())
        *
        mat4d(
            s[0],   u[0],   -f[0],  0.0,
            s[1],   u[1],   -f[1],  0.0,
            s[2],   u[2],   -f[2],  0.0,
            0.0,    0.0,    0.0,    1.0
        );
}

// XXX untested as of yet!
inline mat4d
mat4d::frustum(double left, double right, double bottom, double top, double near, double far)
{
    const double A = (right+left) / (right-left);
    const double B = (top+bottom) / (top-bottom);
    const double C = -(far+near) / (far-near);
    const double D = -(2.0*far*near) / (far-near);

    return mat4d(
        2.0*near*far/(right-left),  0.0,                        0.0,    0.0,
        0.0,                        2*near*far*(top-bottom),    0.0,    0.0,
        A,                          B,                          C,      -1.0,
        0.0,                        0.0,                        D,      0.0
    );
}

//
// constructors
//

inline
mat4d::mat4d()
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            d[r][c] = 0.0;
}

inline
mat4d::mat4d(double h)
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            d[r][c] = h;
}

inline
mat4d::mat4d(
        double e00, double e01, double e02, double e03,
        double e10, double e11, double e12, double e13,
        double e20, double e21, double e22, double e23,
        double e30, double e31, double e32, double e33)
{
    d[0][0] = e00;
    d[0][1] = e01;
    d[0][2] = e02;
    d[0][3] = e03;

    d[1][0] = e10;
    d[1][1] = e11;
    d[1][2] = e12;
    d[1][3] = e13;

    d[2][0] = e20;
    d[2][1] = e21;
    d[2][2] = e22;
    d[2][3] = e23;

    d[3][0] = e30;
    d[3][1] = e31;
    d[3][2] = e32;
    d[3][3] = e33;
}

inline
mat4d::mat4d(float *m)
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            d[r][c] = m[4*r+c];
}

inline
mat4d::mat4d(double *m)
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            d[r][c] = m[4*r+c];
}

inline
mat4d::mat4d(const mat4d& m)
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            d[r][c] = m.d[r][c];
}

//
// instance methods
//

inline void
mat4d::set(float m[])
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            d[r][c] = m[r*4+c];
}

inline void
mat4d::set(double m[])
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            d[r][c] = m[r*4+c];
}


inline void
mat4d::set_zero()
{
    int r, c;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            d[r][c] = 0.0;
}

inline void
mat4d::set_identity()
{
    set_zero();
    d[0][0] = d[1][1] = d[2][2] = d[3][3] = 1.0;
}

inline void
mat4d::set_diagonal(double e00, double e11, double e22, double e33)
{
    set_zero();
    d[0][0] = e00;
    d[1][1] = e11;
    d[2][2] = e22;
    d[3][3] = e33;
}

//
// operations
//

inline mat4d
mat4d::operator-() const
{
    return mat4d(
        -d[0][0], -d[0][1], -d[0][2], -d[0][3],
        -d[1][0], -d[1][1], -d[1][2], -d[1][3],
        -d[2][0], -d[2][1], -d[2][2], -d[2][3],
        -d[3][0], -d[3][1], -d[3][2], -d[3][3]
    );
}

inline void
mat4d::negate()
{
    int r, c;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            d[r][c] = -d[r][c];
}

// m * d -> m
inline mat4d
operator*(const mat4d& m, double d)
{
    return mat4d(
        m.d[0][0]*d, m.d[0][1]*d, m.d[0][2]*d, m.d[0][3]*d,
        m.d[1][0]*d, m.d[1][1]*d, m.d[1][2]*d, m.d[1][3]*d,
        m.d[2][0]*d, m.d[2][1]*d, m.d[2][2]*d, m.d[2][3]*d,
        m.d[3][0]*d, m.d[3][1]*d, m.d[3][2]*d, m.d[3][3]*d
    );
}

// d * m -> m
inline mat4d
operator*(double d, const mat4d& m)
{
    return mat4d(
        m.d[0][0]*d, m.d[0][1]*d, m.d[0][2]*d, m.d[0][3]*d,
        m.d[1][0]*d, m.d[1][1]*d, m.d[1][2]*d, m.d[1][3]*d,
        m.d[2][0]*d, m.d[2][1]*d, m.d[2][2]*d, m.d[2][3]*d,
        m.d[3][0]*d, m.d[3][1]*d, m.d[3][2]*d, m.d[3][3]*d
    );
}

// m / d -> m
inline mat4d
operator/(const mat4d& m, double d)
{
    return mat4d(
        m.d[0][0]/d, m.d[0][1]/d, m.d[0][2]/d, m.d[0][3]/d,
        m.d[1][0]/d, m.d[1][1]/d, m.d[1][2]/d, m.d[1][3]/d,
        m.d[2][0]/d, m.d[2][1]/d, m.d[2][2]/d, m.d[2][3]/d,
        m.d[3][0]/d, m.d[3][1]/d, m.d[3][2]/d, m.d[3][3]/d
    );
}


inline void
mat4d::operator*=(double s)
{
    int r, c;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            d[r][c] *= s;
}

inline void
mat4d::operator/=(double s)
{
    int r, c;
    double invs = 1.0 / s;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            d[r][c] *= invs;
}

// m1 + m2 -> m3
inline mat4d
operator+(const mat4d& m1, const mat4d& m2)
{
    mat4d    res;
    int     r, c;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            res(r, c) = m1(r, c) + m2(r, c);

    return res;
}

// m1 - m2 -> m3
inline mat4d
operator-(const mat4d& m1, const mat4d& m2)
{
    mat4d    res;
    int     r, c;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            res(r, c) = m1(r, c) - m2(r, c);

    return res;
}

// v * m -> v
inline vec4d
operator*(const vec4d v, const mat4d& m)
{
    vec4d res;

    res.x() = v.x()*m(0,0) + v.y()*m(1,0) + v.z()*m(2,0) + v.w()*m(3,0);
    res.y() = v.x()*m(0,1) + v.y()*m(1,1) + v.z()*m(2,1) + v.w()*m(3,1);
    res.z() = v.x()*m(0,2) + v.y()*m(1,2) + v.z()*m(2,2) + v.w()*m(3,2);
    res.w() = v.x()*m(0,3) + v.y()*m(1,3) + v.z()*m(2,3) + v.w()*m(3,3);

    return res;
}

// m1 * m2 -> m
inline mat4d
operator*(const mat4d& m1, const mat4d& m2)
{
    mat4d    res;
    double  sum;
    int     r, c, i;

    for (r = 0; r < 4; r++)
    {
        for (c = 0; c < 4; c++)
        {
            sum = 0.0;
            for (i = 0; i < 4; i++)
                sum += m1.d[r][i] * m2.d[i][c];
            res.d[r][c] = sum;
        }
    }

    return res;
}

inline mat4d
mat4d::transpose() const
{
    return mat4d(
        d[0][0], d[1][0], d[2][0], d[3][0],
        d[0][1], d[1][1], d[2][1], d[3][1],
        d[0][2], d[1][2], d[2][2], d[3][2],
        d[0][3], d[1][3], d[2][3], d[3][3]
    );
}

mat4d
mat4d::inverse() const
{
    mat4d m(*this);					// a copy of this matrix to work on
    mat4d mi = mat4d::identity();

    int d_;      // shadows mat4d::d otherwise
    int r, c;

    double max;
    int maxrow;
    double pivot, temp;
    double  f;

    for (d_ = 0; d_ < 4; d_++)
    {
        //printf("=============================================\n");
        //printf ("d_ = %d\n", d_);

        // find absolute largest element (the pivot) in column d, row d and below only

        max = fabs(m.d[d_][d_]);
        maxrow = d_;

        for (r = d_+1; r < 4; r++)
        {
            if (fabs(m.d[r][d_]) > max)
            {
                max = fabs(m.d[r][d_]);
                maxrow = r;
            }
        }

        // check for zero pivot

        if (mvq::zero(max))
        {
            // Singular matrix, return zero matrix
            return mat4d();
        }

        // swap row containing max value with row d, if necessary

        if (maxrow != d_)
        {
            //printf("(d_=%d) swapping row %d and %d\n", d_, d_, maxrow);

            // elements of m in columns left of column d (in row d and lower) are zero,
            // no need to swap.
            // this is not the case for mi, of course

            // swap rows in m
            for (c = d_; c < 4; c++)
            {
                temp = m.d[maxrow][c];
                m.d[maxrow][c] = m.d[d_][c];
                m.d[d_][c] = temp;
            }

            // swap rows in mi
            for (c = 0; c < 4; c++)
            {
                temp = mi.d[maxrow][c];
                mi.d[maxrow][c] = mi.d[d_][c];
                mi.d[d_][c] = temp;
            }

            //printf("m now:\n");
            //m.dump();
            //printf("mi now:\n");
            //mi.dump();
        }

        pivot = m.d[d_][d_];

        // make all elements of m in column d (except row d) zero, by
        // subtracting row d scaled by the appropiate amount. similarly alter rows of mi.
        // when the element under inspection is already zero don't touch that row.

        for (r = 0; r < 4; r++)
        {
            // if value is already zero, do nothing
            if (r != d_ && !mvq::zero(m.d[r][d_]))
            {
                // calculate scale factor
                f = m.d[r][d_] / pivot;

                // subtract row
                m.d[r][d_] = 0.0;
                for (c = d_+1; c < 4; c++)
                    m.d[r][c] -= f*m.d[d_][c];
                for (c = 0; c < 4; c++)
                    mi.d[r][c] -= f*mi.d[d_][c];
            }

        }

        //printf("m after making elements in col d zero:\n");
        //m.dump();
        //printf("mi:\n");
        //mi.dump();

        // divide row d by pivot

        for (c = d_+1; c < 4; c++)
            m.d[d_][c] /= pivot;
        for (c = 0; c < 4; c++)
            mi.d[d_][c] /= pivot;

        m.d[d_][d_] = 1.0;

        //printf("m after dividing by pivot:\n");
        //m.dump();
        //printf("mi after dividing by pivot:\n");
        //mi.dump();
    }

    return mi;
}

inline bool
mat4d::is_identity() const
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
        {
            if (r != c && !mvq::zero(d[r][c]))
                return false;
            else if (r == c && !mvq::equal(d[r][c], 1.0))
                return false;
        }

    return true;
}

inline bool
mat4d::is_zero() const
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            if (!mvq::zero(d[r][c]))
                return false;

    return true;
}

// transform homogenous point (i.e. assumes p.w == 1) AND normalizes!
// PRE-MULTIPLIES VECTOR WITH MATRIX, i.e. v*M!
inline vec3f
mat4d::transform(float x, float y, float z) const
{
    float inv_d = 1.0 / (x*d[0][3] + y*d[1][3] + z*d[2][3] + d[3][3]);

    return vec3f((x*d[0][0] + y*d[1][0] + z*d[2][0] + d[3][0]) * inv_d,
                 (x*d[0][1] + y*d[1][1] + z*d[2][1] + d[3][1]) * inv_d,
                 (x*d[0][2] + y*d[1][2] + z*d[2][2] + d[3][2]) * inv_d);
}

inline vec3d
mat4d::transform(double x, double y, double z) const
{
    double inv_d = 1.0 / (x*d[0][3] + y*d[1][3] + z*d[2][3] + d[3][3]);

    return vec3d((x*d[0][0] + y*d[1][0] + z*d[2][0] + d[3][0]) * inv_d,
                 (x*d[0][1] + y*d[1][1] + z*d[2][1] + d[3][1]) * inv_d,
                 (x*d[0][2] + y*d[1][2] + z*d[2][2] + d[3][2]) * inv_d);
}


inline vec3f
mat4d::transform(vec3f p) const
{
    float inv_d = 1.0 / (p.x()*d[0][3] + p.y()*d[1][3] + p.z()*d[2][3] + d[3][3]);

    return vec3f((p.x()*d[0][0] + p.y()*d[1][0] + p.z()*d[2][0] + d[3][0]) * inv_d,
                 (p.x()*d[0][1] + p.y()*d[1][1] + p.z()*d[2][1] + d[3][1]) * inv_d,
                 (p.x()*d[0][2] + p.y()*d[1][2] + p.z()*d[2][2] + d[3][2]) * inv_d);
}


inline vec3d
mat4d::transform(vec3d p) const
{
    double inv_d = 1.0 / (p.x()*d[0][3] + p.y()*d[1][3] + p.z()*d[2][3] + d[3][3]);

    return vec3d((p.x()*d[0][0] + p.y()*d[1][0] + p.z()*d[2][0] + d[3][0]) * inv_d,
                 (p.x()*d[0][1] + p.y()*d[1][1] + p.z()*d[2][1] + d[3][1]) * inv_d,
                 (p.x()*d[0][2] + p.y()*d[1][2] + p.z()*d[2][2] + d[3][2]) * inv_d);
}

// transform vector, i.e. disregard translation (assumes v.w == 0)

inline vec3f
mat4d::vtransform(float x, float y, float z) const
{
    return vec3f(
        x*d[0][0] + y*d[1][0] + z*d[2][0],
        x*d[0][1] + y*d[1][1] + z*d[2][1],
        x*d[0][2] + y*d[1][2] + z*d[2][2]
    );
}

inline vec3d
mat4d::vtransform(double x, double y, double z) const
{
    return vec3d(
        x*d[0][0] + y*d[1][0] + z*d[2][0],
        x*d[0][1] + y*d[1][1] + z*d[2][1],
        x*d[0][2] + y*d[1][2] + z*d[2][2]
    );
}

inline vec3f
mat4d::vtransform(const vec3f v) const
{
    return vec3f(
        v.x()*d[0][0] + v.y()*d[1][0] + v.z()*d[2][0],
        v.x()*d[0][1] + v.y()*d[1][1] + v.z()*d[2][1],
        v.x()*d[0][2] + v.y()*d[1][2] + v.z()*d[2][2]
    );
}

inline vec3d
mat4d::vtransform(const vec3d v) const
{
    return vec3d(
        v.x()*d[0][0] + v.y()*d[1][0] + v.z()*d[2][0],
        v.x()*d[0][1] + v.y()*d[1][1] + v.z()*d[2][1],
        v.x()*d[0][2] + v.y()*d[1][2] + v.z()*d[2][2]
    );
}


inline  mat4d
mat4d::normal_transform() const
{
    return inverse().transpose();
}

//
// Set default mat4 type
//

typedef mat4d   mat4;

} // namespace mvq

#endif
