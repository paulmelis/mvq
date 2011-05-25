#include <sstream>
#include <ostream>
#include <boost/python.hpp>
#include "mvq.h"

using namespace boost::python;
using namespace mvq;

static std::string
strwrapper(const mat4d m)
{
    std::ostringstream ss;
    ss.precision(6);
    for (int r = 0; r < 4; r++)
    {
        ss << "| ";
        ss << m(r,0) << " " << m(r,1) << " " << m(r,2) << " " << m(r,3);
        ss << "|";
        if (r <= 2)
            ss << std::endl;
    }
    return ss.str();
}

void
init_mat4d()
{
    class_<mat4d>("mat4d")

        .def("identity", &mat4d::identity)
        .staticmethod("identity")

        .def("diagonal", &mat4d::diagonal)
        .staticmethod("diagonal")

        .def("translate", &mat4d::translate)
        .staticmethod("translate")

        .def("scale", &mat4d::scale)
        .staticmethod("scale")

        .def("rotate_x", &mat4d::rotate_x)
        .staticmethod("rotate_x")
        .def("rotate_y", &mat4d::rotate_y)
        .staticmethod("rotate_y")
        .def("rotate_z", &mat4d::rotate_z)
        .staticmethod("rotate_z")

        .def(init<double>())
        /*
        XXX why isn't this one enabled?
        .def(init<double,double,double,double,
                double,double,double,double,
                double,double,double,double,
                double,double,double,double>())
        */

        .def("set", &mat4d::set)
        .def("get", &mat4d::get)

        .def("set_zero", &mat4d::set_zero)
        .def("set_identity", &mat4d::set_identity)
        .def("set_diagonal", &mat4d::set_diagonal)

        .def("negate", &mat4d::negate)
        .def(-self)

        .def("transpose", &mat4d::transpose)
        .def("inverse", &mat4d::inverse)

        .def("is_zero", &mat4d::is_zero)
        .def("is_identity", &mat4d::is_identity)

        .def(self + self)
        .def(self - self)
        .def(self * self)
        //.def(self ^ self)

        .def(self * double())
        .def(double() * self)
        .def(self / double())

        .def(self *= double())
        .def(self /= double())

        .def(other<vec4d>() * self)

        .def("transform", (vec3f (mat4d::*)(float, float, float) const) &mat4d::transform)
        .def("transform", (vec3d (mat4d::*)(double, double, double) const) &mat4d::transform)
        .def("transform", (vec3f (mat4d::*)(const vec3f) const) &mat4d::transform)
        .def("transform", (vec3d (mat4d::*)(const vec3d) const) &mat4d::transform)
        .def("vtransform", (vec3f (mat4d::*)(float, float, float) const) &mat4d::vtransform)
        .def("vtransform", (vec3d (mat4d::*)(double, double, double) const) &mat4d::vtransform)
        .def("vtransform", (vec3f (mat4d::*)(const vec3f) const) &mat4d::vtransform)
        .def("vtransform", (vec3d (mat4d::*)(const vec3d) const) &mat4d::vtransform)
        .def("normal_transform", &mat4d::normal_transform)

        //.def(self_ns::str(self))
        .def("__str__", strwrapper)
        .def("__repr__", strwrapper)
    ;
}
