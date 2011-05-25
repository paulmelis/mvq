#ifndef MVQ_BOUNDS_H
#define MVQ_BOUNDS_H

#include "vec3.h"
#include "mathutil.h"

using namespace mvq;

namespace mvq {

/*
 * Bounding rectangle (a 2D bounding "box")
 */

class brect
{
public:

    brect()
    {
        m_min[0] = m_min[1] = mvq::infinity;
        m_max[0] = m_max[1] = -mvq::infinity;

        m_empty = true;
    }

    brect(double min_x, double min_y, double max_x, double max_y)
    {
        m_min[0] = min_x;
        m_min[1] = min_y;

        m_max[0] = max_x;
        m_max[1] = max_y;

        m_empty = false;
    }

    void set(double min_x, double min_y, double max_x, double max_y)
    {
        m_min[0] = min_x;
        m_min[1] = min_y;

        m_max[0] = max_x;
        m_max[1] = max_y;

        m_empty = false;
    }

    void update(double x, double y)
    {
        if (x < m_min[0]) m_min[0] = x;
        if (y < m_min[1]) m_min[1] = y;
        if (x > m_max[0]) m_max[0] = x;
        if (y > m_max[1]) m_max[1] = y;

        m_empty = false;
    }

    void update(const double p[2])
    {
        update(p[0], p[1]);
    }

    void update(const brect& b)
    {
        update(b.min_x(), b.min_y());
        update(b.min_x(), b.max_y());
        update(b.max_x(), b.min_y());
        update(b.max_x(), b.max_y());
    }

    bool includes(double x, double y) const
    {
        return (x >= m_min[0] && x <= m_max[0] && y >= m_min[1] && y <= m_max[1]);

    }

    bool includes(const double p[2]) const
    {
        return (p[0] >= m_min[0] && p[0] <= m_max[0] && p[1] >= m_min[1] && p[1] <= m_max[1]);

    }

    void enlarge(double e)
    {
        m_min[0] -= e;
        m_min[1] -= e;
        m_max[0] += e;
        m_max[1] += e;
    }

    bool    empty() const   { return m_empty; }

    double  min_x() const    { return m_min[0]; }
    double  min_y() const    { return m_min[1]; }
    double  max_x() const    { return m_max[0]; }
    double  max_y() const    { return m_max[1]; }

    double  center_x() const { return 0.5*(m_max[0]+m_min[0]); }
    double  center_y() const { return 0.5*(m_max[1]+m_min[1]); }

    double  size_x() const   { return m_max[0]-m_min[0]; }
    double  size_y() const   { return m_max[1]-m_min[1]; }

protected:

    bool    m_empty;
    double  m_min[2];
    double  m_max[2];
};

/*
 * forward class declaration
 */

class bsphere;

/*
 * bbox
 */

class bbox
{
public:

    static bbox combine(const bbox& b1, const bbox& b2);

public:

    bbox()
    {
        bounds[0][0] = bounds[0][1] = bounds[0][2] = mvq::infinity;
        bounds[1][0] = bounds[1][1] = bounds[1][2] = -mvq::infinity;

        m_empty = true;
    }

    bbox(const bbox& other)
    {
        bounds[0] = other.min();
        bounds[1] = other.max();
        m_empty = other.is_empty();
    }

    bbox& operator=(const bbox& other)
    {
        bounds[0] = other.min();
        bounds[1] = other.max();
        m_empty = other.is_empty();
        return *this;
    }

    bbox(mvq::vec3d min, mvq::vec3d max)
    {
        bounds[0] = min;
        bounds[1] = max;

        m_empty = false;
    }

    bbox(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z)
    {
        bounds[0][0] = min_x;
        bounds[0][1] = min_y;
        bounds[0][2] = min_z;

        bounds[1][0] = max_x;
        bounds[1][1] = max_y;
        bounds[1][2] = max_z;

        m_empty = false;
    }

    // bbox centered at (x,y,z) with given size in all three dimensions
    bbox(double x, double y, double z, double size)
    {
        double	hsize = size / 2.0;

        bounds[0][0] = x - hsize;
        bounds[0][1] = y - hsize;
        bounds[0][2] = z - hsize;

        bounds[1][0] = x + hsize;
        bounds[1][1] = y + hsize;
        bounds[1][2] = z + hsize;

        m_empty = false;
    }

public:

    void set(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z)
    {
        bounds[0][0] = min_x;
        bounds[0][1] = min_y;
        bounds[0][2] = min_z;

        bounds[1][0] = max_x;
        bounds[1][1] = max_y;
        bounds[1][2] = max_z;

        m_empty = false;
    }

    inline void     update(double x, double y, double z);
    inline void     update(const double p[3]);
    inline void     update(const vec3f p);
    inline void     update(const vec3d p);
    inline void     update(const bbox& b);
    inline void     update(const bsphere& b);

    inline bool     includes(double x, double y, double z) const;
    inline bool     includes(const double p[3]) const;
    inline bool     includes(const vec3d p) const;
    inline bool     includes(const bbox& b) const;
    inline bool     includes(const bsphere& b) const;

    inline bool     intersects(const bbox& b) const;
    inline bool     intersects(const bsphere& b) const;

    inline void     clear();
    inline void     clear_to_single_point(const vec3d p);

    inline void     enlarge(double e);

    inline bool     is_empty() const { return m_empty; }

public:

    inline bsphere 	as_bsphere() const;

public:

    bool    empty()     const   { return m_empty; }

    double  min_x()     const   { return bounds[0][0]; }
    double  min_y()     const   { return bounds[0][1]; }
    double  min_z()     const   { return bounds[0][2]; }
    double  max_x()     const   { return bounds[1][0]; }
    double  max_y()     const   { return bounds[1][1]; }
    double  max_z()     const   { return bounds[1][2]; }

    const vec3d min()   const   { return bounds[0]; }
    const vec3d max()   const   { return bounds[1]; }

    vec3d&  min()               { return bounds[0]; }
    vec3d&  max()               { return bounds[1]; }

    double  center_x()  const   { return 0.5*(bounds[1][0]+bounds[0][0]); }
    double  center_y()  const   { return 0.5*(bounds[1][1]+bounds[0][1]); }
    double  center_z()  const   { return 0.5*(bounds[1][2]+bounds[0][2]); }

    vec3d   center()    const   { return vec3d(0.5*(bounds[1][0]+bounds[0][0]), 0.5*(bounds[1][1]+bounds[0][1]), 0.5*(bounds[1][2]+bounds[0][2])); }

    double  size_x()    const   { return bounds[1][0]-bounds[0][0]; }
    double  size_y()    const   { return bounds[1][1]-bounds[0][1]; }
    double  size_z()    const   { return bounds[1][2]-bounds[0][2]; }

    vec3d   size()      const   { return bounds[1] - bounds[0]; }

    double	max_size()  const   { if (size_x() < size_y()) return (size_y() > size_z() ? size_y() : size_z());
                                    else return (size_x() > size_z() ? size_x() : size_y()); }

    vec3d       &operator[](short idx)      { return bounds[idx]; }
    const vec3d operator[](short idx) const { return bounds[idx]; }

protected:

    vec3d   bounds[2];
    bool    m_empty;
};

/*
 * bsphere
 */

class bsphere
{
public:

    bsphere()
    {
        m_center = vec3d(0.0, 0.0, 0.0);
        m_radius = 0.0;

        m_empty = true;
    }

    bsphere(vec3d center, double radius)
    {
        m_center = center;
        m_radius = radius;

        m_empty = false;
    }

    bsphere(double cx, double cy, double cz, double r)
    {
        m_center = vec3d(cx, cy, cz);
        m_radius = r;

        m_empty = false;
    }

public:

    inline void    	update(double x, double y, double z);
    inline void    	update(const double p[3]);
    inline void    	update(const vec3f p);
    inline void    	update(const vec3d p);
    inline void    	update(const bsphere& b);
    inline void		update(const bbox& b);

    inline bool		includes(double x, double y, double z) const;
    inline bool		includes(const double p[3]) const;
    inline bool		includes(const vec3d p) const;
    inline bool		includes(const bsphere& b) const;

    inline bool		intersects(const bsphere& b) const;
    inline bool		intersects(const bbox& b) const;

public:

    inline bbox     as_bbox() const;

public:

    bool    empty() 	const   { return m_empty; }

    double 	min_x() 	const	{ return m_center.x() - m_radius; }
    double 	min_y() 	const	{ return m_center.y() - m_radius; }
    double 	min_z() 	const	{ return m_center.z() - m_radius; }
    double 	max_x() 	const	{ return m_center.x() + m_radius; }
    double 	max_y() 	const	{ return m_center.y() + m_radius; }
    double 	max_z() 	const	{ return m_center.z() + m_radius; }

    //vec3d    min() const     { return vec3(m_min[0], m_min[1], m_min[2]); }
    //vec3d    max() const     { return vec3(m_max[0], m_max[1], m_max[2]); }

    double  center_x() 	const 	{ return m_center.x(); }
    double  center_y() 	const 	{ return m_center.y(); }
    double  center_z() 	const 	{ return m_center.z(); }

    vec3d   center() 	const  	{ return m_center; }

    double	radius() 	const	{ return m_radius; }

    double  size_x() 	const   { return m_radius*2.0; }
    double  size_y() 	const   { return m_radius*2.0; }
    double  size_z() 	const   { return m_radius*2.0; }

    vec3d   size() 		const	{ return vec3d(m_radius*2.0, m_radius*2.0, m_radius*2.0); }

    void    enlarge(double e)   { m_radius += e; }

protected:

    bool	m_empty;
    vec3d	m_center;
    double	m_radius;
};


/*
 * bbox
 *
 * inline methods
 */

inline bbox
bbox::combine(const bbox& b1, const bbox& b2)
{
    double  v1, v2;
    bbox    res;

    for (short i = 0; i < 3; i++)
    {
        v1 = b1.min()[i];
        v2 = b2.min()[i];
        if (v1 < v2) { res.min()[i] = v1; } else { res.min()[i] = v2; }

        v1 = b1.max()[i];
        v2 = b2.max()[i];
        if (v1 > v2) { res.max()[i] = v1; } else { res.max()[i] = v2; }
    }

    return res;
}

// update

inline void
bbox::update(double x, double y, double z)
{
    if (x < bounds[0][0]) bounds[0][0] = x;
    if (y < bounds[0][1]) bounds[0][1] = y;
    if (z < bounds[0][2]) bounds[0][2] = z;

    if (x > bounds[1][0]) bounds[1][0] = x;
    if (y > bounds[1][1]) bounds[1][1] = y;
    if (z > bounds[1][2]) bounds[1][2] = z;

    m_empty = false;
}

inline void
bbox::update(const double p[3])
{
    update(p[0], p[1], p[2]);
}

inline void
bbox::update(const vec3f p)
{
    update(p.x(), p.y(), p.z());
}

inline void
bbox::update(const vec3d p)
{
    update(p.x(), p.y(), p.z());
}

inline void
bbox::update(const bbox& b)
{
    double x[2] = { b.min_x(), b.max_x() };
    double y[2] = { b.min_y(), b.max_y() };
    double z[2] = { b.min_z(), b.max_z() };

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                update(x[i], y[j], z[k]);
}

inline void
bbox::update(const bsphere& b)
{
    double x[2] = { b.min_x(), b.max_x() };
    double y[2] = { b.min_y(), b.max_y() };
    double z[2] = { b.min_z(), b.max_z() };

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                update(x[i], y[j], z[k]);
}

// includes

inline bool
bbox::includes(double x, double y, double z) const
{
    if (m_empty)
        return false;

    return (x >= bounds[0][0] && x <= bounds[1][0]
            &&
            y >= bounds[0][1] && y <= bounds[1][1]
            &&
            z >= bounds[0][2] && z <= bounds[1][2]);

}

inline bool
bbox::includes(const double p[3]) const
{
    if (m_empty)
        return false;

    return (p[0] >= bounds[0][0] && p[0] <= bounds[1][0]
            &&
            p[1] >= bounds[0][1] && p[1] <= bounds[1][1]
            &&
            p[2] >= bounds[0][2] && p[2] <= bounds[1][2]);

}

inline bool
bbox::includes(const vec3d p) const
{
    if (m_empty)
        return false;

    return (p.x() >= bounds[0][0] && p.x() <= bounds[1][0]
            &&
            p.y() >= bounds[0][1] && p.y() <= bounds[1][1]
            &&
            p.z() >= bounds[0][2] && p.z() <= bounds[1][2]);

}

inline bool
bbox::includes(const bbox& b) const
{
    if (m_empty || b.empty())
        return false;

    return (
        b.min_x() >= min_x() && b.max_x() <= max_x()
        &&
        b.min_y() >= min_y() && b.max_y() <= max_y()
        &&
        b.min_z() >= min_z() && b.max_z() <= max_z()
    );
}

inline bool
bbox::includes(const bsphere& b) const
{
    if (m_empty || b.empty())
        return false;

    return (includes(b.center()+b.radius()*vec3d(-1,-1,-1))
            &&
            includes(b.center()+b.radius()*vec3d(1,1,1)));
}

// intersects

inline bool
bbox::intersects(const bbox& b) const
{
    if (m_empty || b.empty())
        return false;

    // XXX check this one!
    return (
        ((b.min_x() >= min_x() && b.min_x() <= max_x())
        ||
        (b.max_x() >= min_x() && b.max_x() <= max_x()))
        &&
        ((b.min_y() >= min_y() && b.min_y() <= max_y())
        ||
        (b.max_y() >= min_y() && b.max_y() <= max_y()))
        &&
        ((b.min_z() >= min_z() && b.min_z() <= max_z())
        ||
        (b.max_z() >= min_z() && b.max_z() <= max_z()))
    );
}

inline bool
bbox::intersects(const bsphere& b) const
{
    if (m_empty || b.empty())
        return false;

    if (b.center().x() < min_x())
    {
        if (min_x() - b.center().x() >= b.radius())
            return false;
    }
    else if (b.center().x() > max_x())
    {
        if (b.center().x() - max_x() >= b.radius())
            return false;
    }
    // else: min_x() <= b.center().x <= max_x()

    if (b.center().y() < min_y())
    {
        if (min_y() - b.center().y() >= b.radius())
            return false;
    }
    else if (b.center().y() > max_y())
    {
        if (b.center().y() - max_y() >= b.radius())
            return false;
    }
    // else: min_y() <= b.center().yx <= max_y()

    if (b.center().z() < min_z())
    {
        if (min_z() - b.center().z() >= b.radius())
            return false;
    }
    else if (b.center().z() > max_z())
    {
        if (b.center().z() - max_z() >= b.radius())
            return false;
    }
    // else: min_z() <= b.center().z <= max_z()

    return true;
}

inline void
bbox::clear()
{
    bounds[0][0] = bounds[0][1] = bounds[0][2] = mvq::infinity;
    bounds[1][0] = bounds[1][1] = bounds[1][2] = -mvq::infinity;

    m_empty = true;
}

inline void
bbox::clear_to_single_point(const vec3d p)
{
    bounds[0][0] = bounds[1][0] = p.x();
    bounds[0][1] = bounds[1][1] = p.y();
    bounds[0][2] = bounds[1][2] = p.z();
    m_empty = false;
}

inline void
bbox::enlarge(double e)
{
    bounds[0][0] -= e;
    bounds[0][1] -= e;
    bounds[0][2] -= e;
    bounds[1][0] += e;
    bounds[1][1] += e;
    bounds[1][2] += e;
}


// as_bsphere

inline bsphere
bbox::as_bsphere() const
{
    if (m_empty)
        return bsphere();

    double max_size = size_x();
    if (size_y() > max_size) max_size = size_y();
    if (size_z() > max_size) max_size = size_z();

    return bsphere(center(), 0.5*sqrt(3.0)*max_size);
}

/*
 * bsphere
 *
 * inline methods
 */

// update

inline void
bsphere::update(double x, double y, double z)
{
    // XXX we can probably generate a tighter bound
    vec3d 	p = vec3d(x, y, z);

    if (m_center.x()*m_center.x() + m_center.y()*m_center.y() + m_center.z()*m_center.z() > m_radius*m_radius)
    {
        vec3d	cp = m_center - vec3d(x, y, z);
        vec3d 	q = m_center + m_radius*cp.normalized();

        // new center, halfway between p and q
        m_center = 0.5*(p + q);

        // new radius
        m_radius = (p - m_center).length();
    }

    m_empty = false;
}

inline void
bsphere::update(const double p[3])
{
    update(p[0], p[1], p[2]);
}

inline void
bsphere::update(const vec3f p)
{
    update(p.x(), p.y(), p.z());
}

inline void
bsphere::update(const vec3d p)
{
    update(p.x(), p.y(), p.z());
}

inline void
bsphere::update(const bsphere& b)
{
    // XXX can probably generate a tighter bound

    vec3d	cp = (b.center() - m_center).normalized();
    vec3d	q = b.center() + cp * b.radius();

    if (!includes(q))
    {
        vec3d 	r = m_center - cp * m_radius;

        m_center = 0.5*(q + r);
        m_radius = (r - m_center).length();
    }

    m_empty = false;
}

inline void
bsphere::update(const bbox& b)
{
    update(b.as_bsphere());
}

// includes

inline bool
bsphere::includes(double x, double y, double z) const
{
    if (m_empty)
        return false;

    vec3d	d = m_center - vec3d(x, y, z);

    return (d.x()*d.x() + d.y()*d.y() + d.z()*d.z() <= m_radius*m_radius);
}

inline bool
bsphere::includes(const double p[3]) const
{
    if (m_empty)
        return false;

    vec3d	d = m_center - vec3d(p);

    return (d.x()*d.x() + d.y()*d.y() + d.z()*d.z() <= m_radius*m_radius);
}

inline bool
bsphere::includes(const vec3d p) const
{
    if (m_empty)
        return false;

    vec3d	d = m_center - p;

    return (d.x()*d.x() + d.y()*d.y() + d.z()*d.z() <= m_radius*m_radius);
}

inline bool
bsphere::includes(const bsphere& b) const
{
    if (m_empty || b.empty())
        return false;

    vec3d	d = b.center() - m_center;

    if (d.length2() < mvq::epsilon)
        return b.radius() <= m_radius;

    vec3d	p = b.center() + d.normalized() * b.radius();

    return includes(p);
}

// intersects

inline bool
bsphere::intersects(const bsphere& b) const
{
    if (m_empty || b.empty())
        return false;

    double d = (m_center - b.center()).length();

    return (d < m_radius + b.radius());
}

inline bool
bsphere::intersects(const bbox& b) const
{
    if (m_empty || b.empty())
        return false;

    return b.intersects(*this);
}

// as_bbox

inline bbox
bsphere::as_bbox() const
{
    if (m_empty)
        return bbox();

    return bbox(m_center.x()-m_radius, m_center.y()-m_radius, m_center.z()-m_radius,
        m_center.x()+m_radius, m_center.y()+m_radius, m_center.z()+m_radius);
}

}	// namespace mvq

#endif
