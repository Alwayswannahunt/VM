
// 2D vectors

template<typename A, typename B>
struct vect {
    A x;
    B y;
};

template<typename A, typename B>
vect<A, B> operator * (double a, vect<A, B> r) {
    return { a * r.x, a * r.y };
}

template<typename A, typename B>
vect<A, B> operator * (vect<A, B> r, double a) {
    return{ a * r.x, a * r.y };
}

template<typename A, typename B>
vect<A, B> operator + (vect<A, B> r_1, vect<A, B> r_2) {
    return{ r_1.x + r_2.x, r_1.y + r_2.y };
}

template<typename A, typename B>
vect<A, B> operator - (vect<A, B> r_1, vect<A, B> r_2) {
    return{ r_1.x - r_2.x, r_1.y - r_2.y };
}

double modulVect(vect<double, double>);

// 3D vectors for approximation

template<typename A, typename B, typename C>
struct vect3D {
    A a0;
    B a1;
    C b1;
};

template<typename A, typename B, typename C>
vect3D<A, B, C> operator + (vect3D<A, B, C> r_1, vect3D<A, B, C> r_2) {
    return{ r_1.a0 + r_2.a0 , r_1.a1 + r_2.a1 , r_1.b1 + r_1.b1 };
};

template<typename A, typename B, typename C>
vect3D<A, B, C> operator - (vect3D<A, B, C> r_1, vect3D<A, B, C> r_2) {
    return{ r_1.a0 - r_2.a0 , r_1.a1 - r_2.a1 , r_1.b1 - r_1.b1 };
};

template<typename A, typename B, typename C>
vect3D<A, B, C> operator * (vect3D<A, B, C> r, double a) {
    return { r.a0 * a,r.a1 * a,r.b1 * a };
};

template<typename A, typename B, typename C>
vect3D<A, B, C> operator * (double a, vect3D<A, B, C> r) {
    return { r.a0 * a,r.a1 * a,r.b1 * a };
};