#include "mvq/quat.h"

int main()
{
    mvq::quatf q(1, 2, 3, 4);

    q.print();              printf("\n");
    q.conjugate().print();  printf("\n");
    q.inverse().print();    printf("\n");

    // Cross-check against Quaternions.jl
    mvq::quatf::rotation(mvq::vec3(0,1,0), 90.0f).print();  printf("\n");
    printf("^ 0.7071067811865476, 0.0, 0.7071067811865475, 0.0\n");

    mvq::quatf::rotation(mvq::vec3(1,1,1), 120.0f).print();  printf("\n");
    printf("^ 0.5000000000000001, 0.5, 0.5, 0.5\n");

    mvq::quatf r = mvq::quatf::rotation(mvq::vec3(0,0,1), 90.0f);
    r.print();              printf("\n");
    mvq::vec3 p(1, 0, 0);
    r.transform(p).print(); printf("\n");

    r = mvq::quatf::rotation(mvq::vec3(1,1,1), 180.0f);
    r.print();      printf("\n");
    printf("^ 6.123233995736766e-17, 0.5773502691896258, 0.5773502691896258, 0.5773502691896258\n");
    p.set(1, 1, 0);
    r.transform(p).print(); printf("\n");

    mvq::quatf q1 = mvq::quatf::rotation(mvq::vec3(0,1,0), 90.0f);
    mvq::mat4 R = q1.to_rotation_matrix();
    R.dump();

    return 0;
}