#include <math.h>
#include "vect.h"

double sgn(double num) {
	return (-2 * signbit(num) + 1);
}


double differenceOfAproximation(	// сама наша функция
	vect<double,double>* real_point,  // точки по которым аппроксимируем (первая координата - параметр, вторая - значения корня при этом параметре)
	int length,	// длина массива real_point
	vect3D<double,double,double> parametrs,
	vect3D<double, double, double> prev_parametrs)
{
	double result = 0;
	for (int i = 0; i < length; i++)
	{
		double Q = (1 + parametrs.b1 * real_point[i].x);
		double P = parametrs.a0 + parametrs.a1 * real_point[i].x + parametrs.b1 * real_point[i].x * real_point[i].x / 2;
		double prev_Q = (1 + prev_parametrs.b1 * real_point[i].x);
		result = (real_point[i].y * Q - P) / prev_Q;
	}
	return result;
}


double differenceOfAproximation1D( // одномерное выражение для нашей функции
	vect<double, double> * real_point,  // точки по которым аппроксимируем (первая координата - параметр, вторая - значения корня при этом параметре)
	int length,	// длина массива real_point
	vect3D<double, double, double> parametrs,
	vect3D<double, double, double> prev_parametrs,
	vect3D<double, double, double> grad,
	double lambda) // точка
{
	return (differenceOfAproximation(real_point, length, parametrs + lambda * grad, prev_parametrs));
}

vect3D<double, double, double> gradFunc(
	double (*p_func)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>),
	vect3D<double, double, double> grad_point,
	vect<double, double>* real_point,  // точки по которым аппроксимируем (первая координата - параметр, вторая - значения корня при этом параметре)
	int length,	// длина массива real_point
	vect3D<double, double, double> parametrs,
	vect3D<double, double, double> prev_parametrs
)
{
	const double delta = 0.00002;
	double x = (p_func(real_point, length, parametrs + delta * vect3D<double, double, double>{1, 0, 0}, prev_parametrs) - p_func(real_point, length, parametrs - delta * vect3D<double, double, double>{1, 0, 0}, prev_parametrs)) / (2 * delta);
	double y = (p_func(real_point, length, parametrs + delta * vect3D<double, double, double>{0, 1, 0}, prev_parametrs) - p_func(real_point, length, parametrs - delta * vect3D<double, double, double>{0, 1, 0}, prev_parametrs)) / (2 * delta);
	double z = (p_func(real_point, length, parametrs + delta * vect3D<double, double, double>{0, 0, 1}, prev_parametrs) - p_func(real_point, length, parametrs - delta * vect3D<double, double, double>{0, 0, 1}, prev_parametrs)) / (2 * delta);
	return{x,y,z};
}

double derivative1D(
	double (*p_func)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double),
	vect<double, double>* real_point,  // точки по которым аппроксимируем (первая координата - параметр, вторая - значения корня при этом параметре)
	int length,	// длина массива real_point
	vect3D<double, double, double> parametrs,
	vect3D<double, double, double> prev_parametrs,
	vect3D<double, double, double> grad,
	double lambda) // точка в которой ищем производную
{
	const double delta = 0.00002;
	return ((p_func(real_point, length, parametrs, prev_parametrs, grad, lambda + delta) - p_func(real_point, length, parametrs, prev_parametrs, grad, lambda - delta)) / (2 * delta));
}

double secondDerivative1D(
	double (*p_func)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double),
	vect<double, double>* real_point,  // точки по которым аппроксимируем (первая координата - параметр, вторая - значения корня при этом параметре)
	int length,	// длина массива real_point
	vect3D<double, double, double> parametrs,
	vect3D<double, double, double> prev_parametrs,
	vect3D<double, double, double> grad,
	double lambda) // точка в которой ищем производную)
{
	const double delta = 0.00002;
	return ((p_func(real_point, length, parametrs, prev_parametrs, grad, lambda + delta) + p_func(real_point, length, parametrs, prev_parametrs, grad, lambda - delta) - 2 * p_func(real_point, length, parametrs, prev_parametrs, grad, lambda)) / (delta * delta));
}

double parabolaMin(
	double (*p_func_frst_derv)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double),
	double (*p_func_secnd_derv)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double),
	vect<double, double>* real_point,  // точки по которым аппроксимируем (первая координата - параметр, вторая - значения корня при этом параметре)
	int length,	// длина массива real_point
	vect3D<double, double, double> parametrs,
	vect3D<double, double, double> prev_parametrs,
	vect3D<double, double, double> grad,
	double lambda,
	double prev_lambda)
{
	double scnd_derv  =  p_func_secnd_derv(real_point, length, parametrs, prev_parametrs, grad, lambda);

	return(prev_lambda - p_func_frst_derv(real_point,length,parametrs,prev_parametrs,grad,lambda) / (scnd_derv*signbit(scnd_derv)) );
}

vect3D<double, double, double> approximationOfRoot(vect<double, double>* real_point, int length) {
	return{ 1,1,1 };
}