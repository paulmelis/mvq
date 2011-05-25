#include <boost/python.hpp>
#include "mathutil.h"

using namespace boost::python;

void init_vec2();
void init_vec3();
void init_vec4();

void init_mat4d();

void init_bounds();

BOOST_PYTHON_MODULE(mvq)
{
    def("deg2rad", &mvq::deg2rad);
    def("rad2deg", &mvq::rad2deg);

    def("equal", &mvq::equal);
    def("zero", &mvq::zero);

    init_vec2();
    init_vec3();
    init_vec4();

    init_mat4d();

    init_bounds();
}
