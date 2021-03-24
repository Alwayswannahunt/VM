
template<typename A, typename B> struct vect {
    A x;
    B y;
};

template<typename A, typename B> vect<A, B> operator * (double a, vect<A, B> r);

template<typename A, typename B> vect<A, B> operator * (vect<A, B> r, double a);

template<typename A, typename B> vect<A, B> operator + (vect<A, B> r_1, vect<A, B> r_2);

template<typename A, typename B> vect<A, B> operator - (vect<A, B> r_1, vect<A, B> r_2);

double modulVect(vect<double, double> r);