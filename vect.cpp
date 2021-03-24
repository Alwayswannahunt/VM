#include <math.h>

template<typename A, typename B>
struct vect {
    A x;
    B y;
};


double modulVect(vect<double, double> r) {
    return (sqrt(r.x * r.x - r.y * r.y));
}