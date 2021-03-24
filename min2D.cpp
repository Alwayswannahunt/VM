
#include <iostream>
#include <math.h>
#include "vect.h"

double sgn(double lol) {
    return(2 * (double)(lol > 0) - 1); // опять та самая апроксимация (ну а что! По двум точкам прямую)   Недавно нашёл функцию signbit, надо было через неё.
}

double f(vect<double, double> r) {   // функция минимум которой ищем
    double x = r.x;
    double y = r.y;
    return ((x * y - (cos(x)/ sin(x))) * (x * y - 1.0 / tan(x)) + sin(x * y) * sin(x * y));
}

double f1D(vect<double, double> r0, double lambda, vect<double, double> grad) {   // функция на линии градиента
    double x = r0.x + lambda * grad.x;
    double y = r0.y + lambda * grad.y;
    return ((x * y - 1.0 / tan(x)) * (x * y - 1.0 / tan(x)) + sin(x * y) * sin(x * y));
}


double derivative(   // поиск первой производной по лямбда
        double (*func)(vect<double, double>, double, vect<double, double>),
        vect<double, double> r,
        double lambda, vect<double,
        double> grad)
{
    double const delta = 0.0002;
    return ((func( r,lambda + delta, grad) - func(r, lambda - delta, grad)) / (2 * delta));
}

double secondDerivative(    // вторая производная по лямбда
    double (*func)(vect<double, double>, double, vect<double, double>),
    vect<double, double> r,
    double lambda,
    vect<double, double> grad)
{
    double const delta = 0.0002;
    return ( ( func(r, lambda + delta, grad) + func(r, lambda - delta, grad) - 2.0*func(r, lambda, grad)) / (delta * delta) );
}

vect<double, double> gradFunc(double (*func)(vect<double, double>), vect<double, double> r) { // поиск градиента функции в точке
    double const delta = 0.0002;
    double x = (f({ r.x + delta,r.y }) - f({ r.x - delta,r.y })) / (2 * delta);
    double y = (f({ r.x ,r.y + delta }) - f({ r.x ,r.y - delta })) / (2 * delta);
    return { x, y };
}


double parabolaMin(double (*func)(vect<double, double>, double, vect<double, double>), vect<double, double> r, double prev_lambda, vect<double, double> grad) {
    double scnd_derv = secondDerivative(func, r, prev_lambda, grad);
    return (prev_lambda - derivative(func,r,prev_lambda,grad)/(scnd_derv*sgn(scnd_derv)));
}



vect<double, double> min2D(vect<double, double> r0) {

    double const h = 0.14;  // а есть логика выбора h?
    vect<double, double> prev_r = r0; // методом пристального вглядывания надо выбирать по разные стороны от 0, например (1,1); (-1,-1)
    vect<double, double> r_1 = { 0,0 };
    vect<double, double> r_2 = { 0,0 };

        double prev_lambda = 0;
        double lambda = 0;
        vect<double, double> grad = gradFunc(f,prev_r);
        do
        {
            prev_lambda = lambda;
            lambda = parabolaMin(f1D, prev_r, prev_lambda, grad);
        } while (f1D(prev_r,prev_lambda, grad) > f1D(prev_r, lambda, grad));
        r_2 = prev_r + lambda * grad;
        // обожаю алгоритмы такие, кстати, интересно сколько я буду это фиксить и искать баги?))) (багов не было, если вам интересно)
        prev_r = 4*h * prev_r;

        do
    {
        r_1 = r_2;
        prev_lambda = 0;
        lambda = 0;
        grad = gradFunc(f, prev_r);
        do
        {
            prev_lambda = lambda;
            lambda = parabolaMin(f1D, prev_r, prev_lambda, grad);
        } while (f1D(prev_r, prev_lambda, grad) > f1D(prev_r, lambda, grad));
        r_2 = prev_r + lambda * grad;
        prev_r = r_2 + h * (r_2 - r_1) * sgn(f(r_1) - f(r_2));
    } while (f(r_1) > f(r_2)  || (modulVect(r_2-r_1) > 0.0001));

    return(r_1);
}