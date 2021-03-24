template<typename A, typename B> struct vect {
    A x;
    B y;
};

template<typename A, typename B> vect<A, B> operator * (double a, vect<A, B> r) {
    return { a * r.x, a * r.y };
}

template<typename A, typename B> vect<A, B> operator * (vect<A, B> r, double a) {
    return{ a * r.x, a * r.y };
}

template<typename A, typename B> vect<A, B> operator + (vect<A, B> r_1, vect<A, B> r_2) {
    return{ r_1.x + r_2.x, r_1.y + r_2.y };
}

template<typename A, typename B> vect<A, B> operator - (vect<A, B> r_1, vect<A, B> r_2) {
    return{ r_1.x - r_2.x, r_1.y - r_2.y };
}

double modulVect(vect<double, double> r) {
    return (sqrt(r.x * r.x - r.y * r.y));
}