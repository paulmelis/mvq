#include <sstream>
#include <ostream>
#include <boost/python.hpp>
#include "mvq.h"
#include "vec3.h"

using namespace boost::python;
using namespace mvq;

static std::string
strwrapper(const vec3f v)
{
    std::ostringstream ss;
    ss.precision(6);
    ss << "(" << v.x << "," << v.y << "," << v.z << ")";
    return ss.str();
}

static std::string
strwrapper(const vec3d v)
{
    std::ostringstream ss;
    ss.precision(12);
    ss << "(" << v.x << "," << v.y << "," << v.z << ")";
    return ss.str();
}

void
init_vec3()
{
    class_<vec3f>("vec3f")

        .def(init<float, float, float>())

        .def("set", &vec3f::set)

        .def_readwrite("x", &vec3f::x)
        .def_readwrite("y", &vec3f::y)
        .def_readwrite("z", &vec3f::z)

        .def("length", &vec3f::length)
        .def("abs", &vec3f::abs)

        .def("normalize", &vec3f::normalize)
        .def("normalized", &vec3f::normalized)

        .def("negate", &vec3f::negate)
        .def(-self)

        .def("zeroLength", &vec3f::zeroLength)

        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self / float())
        .def(self ^ self)

        .def(self += self)
        .def(self -= self)
        .def(self *= float())
        .def(self /= float())

        .def(self * float())
        .def(float() * self)

        //.def(self_ns::str(self))
        .def("__str__", (std::string (*) (const vec3f)) &strwrapper)
        .def("__repr__", (std::string (*) (const vec3f)) &strwrapper)
    ;

    class_<vec3d>("vec3d")

        .def(init<double, double, double>())

        .def("set", &vec3d::set)

        .def_readwrite("x", &vec3d::x)
        .def_readwrite("y", &vec3d::y)
        .def_readwrite("z", &vec3d::z)

        .def("length", &vec3d::length)
        .def("abs", &vec3d::abs)

        .def("normalize", &vec3d::normalize)
        .def("normalized", &vec3d::normalized)

        .def("negate", &vec3d::negate)
        .def(-self)

        .def("zeroLength", &vec3d::zeroLength)

        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self / double())
        .def(self ^ self)

        .def(self += self)
        .def(self -= self)
        .def(self *= double())
        .def(self /= double())

        .def(self * double())
        .def(double() * self)

        //.def(self_ns::str(self))
        .def("__str__", (std::string (*) (const vec3d)) &strwrapper)
        .def("__repr__", (std::string (*) (const vec3d)) &strwrapper)
    ;
}
