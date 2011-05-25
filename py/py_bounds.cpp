#include <sstream>
#include <ostream>
#include <boost/python.hpp>
#include "mvq.h"
#include "bounds.h"

using namespace boost::python;
using namespace mvq;

static std::string
strwrapper_bbox(const bbox b)
{
    std::ostringstream ss;
    ss.precision(6);
    if (b.empty())
        ss << "(bbox: empty)";
    else
        ss 	<< "(bbox: "
            << b.min_x() << "," << b.min_y() << "," << b.min_z() << "; "
            << b.max_x() << "," << b.max_y() << "," << b.max_z() << "; "
            << "center=" << b.center_x() << "," << b.center_y() << "," << b.center_z() << "; "
            << "size=" << b.size_x() << "," << b.size_y() << "," << b.size_z()
            << ")";
    return ss.str();
}

static std::string
strwrapper_bsphere(const bsphere b)
{
    std::ostringstream ss;
    ss.precision(6);
    if (b.empty())
        ss << "(bsphere: empty)";
    else
        ss 	<< "(bsphere: "
            << "center=" << b.center_x() << "," << b.center_y() << "," << b.center_z() << "; "
            << "radius=" << b.radius()
            << ")";
    return ss.str();
}

void
init_bounds()
{
    class_<bbox>("bbox")

        .def(init<double, double, double, double, double, double>())
        .def(init<double, double, double, double>())

        .def("set", &bbox::set)

        .def("update", (void (bbox::*)(double,double,double)) &bbox::update)
        //.def("update", (void (bbox::*)(const double[3])) &bbox::update)
        .def("update", (void (bbox::*)(const vec3d)) &bbox::update)
        .def("update", (void (bbox::*)(const bbox&)) &bbox::update)
        .def("update", (void (bbox::*)(const bsphere&)) &bbox::update)

        .def("includes", (bool (bbox::*)(double,double,double) const) &bbox::includes)
        //.def("includes", (bool (bbox::*)(const double[3]) const) &bbox::includes)
        .def("includes", (bool (bbox::*)(const vec3d) const) &bbox::includes)
        .def("includes", (bool (bbox::*)(const bbox&) const) &bbox::includes)
        .def("includes", (bool (bbox::*)(const bsphere&) const) &bbox::includes)

        .def("intersects", (bool (bbox::*)(const bbox&) const) &bbox::intersects)
        .def("intersects", (bool (bbox::*)(const bsphere&) const) &bbox::intersects)


        .def("as_bsphere", &bbox::as_bsphere)

        .add_property("empty", &bbox::empty)

        .add_property("min_x", &bbox::min_x)
        .add_property("min_y", &bbox::min_y)
        .add_property("min_z", &bbox::min_z)
        .add_property("max_x", &bbox::max_x)
        .add_property("max_y", &bbox::max_y)
        .add_property("max_z", &bbox::max_z)

        .add_property("min", &bbox::min)
        .add_property("max", &bbox::max)

        .add_property("center_x", &bbox::center_x)
        .add_property("center_y", &bbox::center_y)
        .add_property("center_z", &bbox::center_z)

        .add_property("center", &bbox::center)

        .add_property("size_x", &bbox::size_x)
        .add_property("size_y", &bbox::size_y)
        .add_property("size_z", &bbox::size_z)

        .add_property("size", &bbox::size)
        .add_property("max_size", &bbox::max_size)

        //.def(self_ns::str(self))
        .def("__str__", strwrapper_bbox)
        .def("__repr__", strwrapper_bbox)
    ;

    class_<bsphere>("bsphere")

        .def(init<double, double, double, double>())
        .def(init<vec3d, double>())


        .def("update", (void (bsphere::*)(double,double,double)) &bsphere::update)
        //.def("update", (void (bsphere::*)(const double[3])) &bsphere::update)
        .def("update", (void (bsphere::*)(const vec3d)) &bsphere::update)
        .def("update", (void (bsphere::*)(const bsphere&)) &bsphere::update)
        .def("update", (void (bsphere::*)(const bbox&)) &bsphere::update)

        .def("includes", (bool (bsphere::*)(double,double,double) const) &bsphere::includes)
        //.def("includes", (bool (bsphere::*)(const double[3]) const) &bsphere::includes)
        .def("includes", (bool (bsphere::*)(const vec3d) const) &bsphere::includes)
        .def("includes", (bool (bsphere::*)(const bsphere&) const) &bsphere::includes)
        //.def("includes", (bool (bsphere::*)(const bbox&) const) &bsphere::includes)

        .def("intersects", (bool (bsphere::*)(const bsphere&) const) &bsphere::intersects)
        .def("intersects", (bool (bsphere::*)(const bbox&) const) &bsphere::intersects)


        .def("as_bbox", &bsphere::as_bbox)


        .add_property("empty", &bsphere::empty)

        .add_property("min_x", &bsphere::min_x)
        .add_property("min_y", &bsphere::min_y)
        .add_property("min_z", &bsphere::min_z)
        .add_property("max_x", &bsphere::max_x)
        .add_property("max_y", &bsphere::max_y)
        .add_property("max_z", &bsphere::max_z)

        //.add_property("min", &bsphere::min)
        //.add_property("max", &bsphere::max)

        .add_property("center_x", &bsphere::center_x)
        .add_property("center_y", &bsphere::center_y)
        .add_property("center_z", &bsphere::center_z)

        .add_property("center", &bsphere::center)
        .add_property("radius", &bsphere::radius)

        .add_property("size_x", &bsphere::size_x)
        .add_property("size_y", &bsphere::size_y)
        .add_property("size_z", &bsphere::size_z)

        .add_property("size", &bsphere::size)

        //.def(self_ns::str(self))
        .def("__str__", strwrapper_bsphere)
        .def("__repr__", strwrapper_bsphere)

    ;
}
