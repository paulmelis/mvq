#include <sstream>
#include <ostream>
#include <boost/python.hpp>
#include "mvq.h"
#include "vec2.h"

using namespace boost::python;
using namespace mvq;

static std::string
strwrapper(const vec2f v)
{
    std::ostringstream ss;
    ss.precision(6);
    ss << "(" << v.x << "," << v.y << ")";
    return ss.str();
}

static std::string
strwrapper(const vec2d v)
{
    std::ostringstream ss;
    ss.precision(12);
    ss << "(" << v.x << "," << v.y << ")";
    return ss.str();
}

void
init_vec2()
{
    class_<vec2f>("vec2f")

        .def(init<float, float>())

        .def("set", &vec2f::set)

        .def_readwrite("x", &vec2f::x)
        .def_readwrite("y", &vec2f::y)

        .def("length", &vec2f::length)
        .def("length2", &vec2f::length2)

        .def("normalize", &vec2f::normalize)
        .def("normalized", &vec2f::normalized)

        .def("negate", &vec2f::negate)
        .def(-self)

        .def("zeroLength", &vec2f::zeroLength)

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
        .def("__str__", (std::string (*) (const vec2f)) &strwrapper)
        .def("__repr__", (std::string (*) (const vec2f)) &strwrapper)
    ;

    class_<vec2d>("vec2d")

        .def(init<double, double>())

        .def("set", &vec2d::set)

        .def_readwrite("x", &vec2d::x)
        .def_readwrite("y", &vec2d::y)

        .def("length", &vec2d::length)
        .def("abs", &vec2d::abs)

        .def("normalize", &vec2d::normalize)
        .def("normalized", &vec2d::normalized)

        .def("negate", &vec2d::negate)
        .def(-self)

        .def("zeroLength", &vec2d::zeroLength)

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
        .def("__str__", (std::string (*) (const vec2d)) &strwrapper)
        .def("__repr__", (std::string (*) (const vec2d)) &strwrapper)
    ;
}
