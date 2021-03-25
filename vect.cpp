#include <math.h>

template<typename A, typename B>
struct vect {
    A x;
    B y;
};

template<typename A, typename B, typename C>
struct vect3D {
    A a0;
    B a1;
    C b1;
};

double modulVect(vect<double, double> r) {
    return (sqrt(r.x * r.x + r.y * r.y));
}

double modulVect3D(vect3D<double,double,double> r) {
    return (sqrt(r.a0 * r.a0 + r.a1 * r.a1 + r.b1*r.b1));
}

double sgn(double num) {
    return (-2 * (double)signbit(num) + 1);
}