#include <iostream>
#include <math.h>

double f(double x, double a) {   // функция корни которой ищем
    return (a/x - exp(x));
}


double derivative(double(* func) (double, double), double x, double a) { // поиск производной от функции в точке
    double const delta = 0.00002;
    return ((func(x+delta, a) - func(x-delta, a))/(2*delta));
}

double secondDerivative(double(*func) (double, double), double x, double a) {
    double const delta = 0.00002;
    return ((func(x + delta, a) + func(x - delta, a) - 2*func(x, a)) / (delta * delta));
}

double taylorRoot(double(*func) (double, double), double prev_x,double a) { // метод Ньютона
    return (prev_x - func(prev_x, a)/derivative(func, prev_x , a) );
}

double parabolaRoot(double(*func) (double, double), double prev_x, double a) {
    double derv_f = derivative(func, prev_x, a);
    double secnd_derv_f = secondDerivative(func, prev_x, a);
    double f = func(prev_x, a);
    double sgn_f = 2*(double)(derv_f > 0) - 1; // функция сигн через приведение сравнения с нулём.
    return (prev_x - (sgn_f * 2 * f / (sgn_f * derv_f + sqrt(derv_f * derv_f - 2 * f * secnd_derv_f))));
}




double root()     // root of a/x = exp(x)
{
    double prev_x; // точка для поиска корня
    double a; // параметр нашего уравнения
    double (*p_f)(double, double) = f; // указатель на функцию
    std::cout.precision(6); // вывод 6-ти знаков после запятой
    std::cin >> a;
    prev_x = a / 2;   // попадает ровно в корень пачему-тааа

    double x = prev_x; // новая точка поиска
    if ((derivative(f, prev_x, a) * derivative(f, prev_x, a)) > 0.0000000001) // проверка вдруг угадили в корень (везёт же!)
    {
        do {
            prev_x = x;
            x = taylorRoot(p_f, prev_x, a);
        } while (((derivative(f, prev_x, a) * derivative(f, prev_x, a)) > 0.00000001) && ((prev_x - x) * (prev_x - x) > 0.00000001)); // 8 порядков
        do {
            prev_x = x;
            x = parabolaRoot(p_f, prev_x, a);
        } while ((prev_x - x) * (prev_x - x) > 0.000000000001); // 12 порядков даёт точность в 6 знаке
    }
    return x;
}

