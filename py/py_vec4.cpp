#include <sstream>
#include <ostream>
#include <boost/python.hpp>
#include "vec4.h"
#include "mat4.h"

using namespace boost::python;
using namespace mvq;

static std::string
strwrapper(const vec4f v)
{
    std::ostringstream ss;
    ss.precision(6);
    ss << "(" << v.x << "," << v.y << "," << v.z << "," << v.w << ")";
    return ss.str();
}

static std::string
strwrapper(const vec4d v)
{
    std::ostringstream ss;
    ss.precision(12);
    ss << "(" << v.x << "," << v.y << "," << v.z << "," << v.w << ")";
    return ss.str();
}

void
init_vec4()
{
    class_<vec4f>("vec4f")

        .def(init<double, double, double, double>())

        .def("set", &vec4f::set)

        .def_readwrite("x", &vec4f::x)
        .def_readwrite("y", &vec4f::y)
        .def_readwrite("z", &vec4f::z)
        .def_readwrite("w", &vec4f::w)

        .def("length", &vec4f::length)
        .def("abs", &vec4f::abs)

        .def("homogenize", &vec4f::homogenize)
        .def("homogenized", &vec4f::homogenized)

        .def("negate", &vec4f::negate)
        .def(-self)

        .def("zeroLength", &vec4f::zeroLength)

        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self / double())

        .def(self += self)
        .def(self -= self)
        .def(self *= double())
        .def(self /= double())

        .def(self * double())
        .def(double() * self)

        //.def(self_ns::str(self))
        .def("__str__", (std::string (*) (const vec4f)) &strwrapper)
        .def("__repr__", (std::string (*) (const vec4f)) &strwrapper)
    ;

    class_<vec4d>("vec4d")

        .def(init<double, double, double, double>())

        .def("set", &vec4d::set)

        .def_readwrite("x", &vec4d::x)
        .def_readwrite("y", &vec4d::y)
        .def_readwrite("z", &vec4d::z)
        .def_readwrite("w", &vec4d::w)

        .def("length", &vec4d::length)
        .def("abs", &vec4d::abs)

        .def("homogenize", &vec4d::homogenize)
        .def("homogenized", &vec4d::homogenized)

        .def("negate", &vec4d::negate)
        .def(-self)

        .def("zeroLength", &vec4d::zeroLength)


        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self / float())

        .def(self += self)
        .def(self -= self)
        .def(self *= float())
        .def(self /= float())

        .def(self * float())
        .def(float() * self)

        //.def(self_ns::str(self))
        .def("__str__", (std::string (*) (const vec4d)) &strwrapper)
        .def("__repr__", (std::string (*) (const vec4d)) &strwrapper)
    ;
}
