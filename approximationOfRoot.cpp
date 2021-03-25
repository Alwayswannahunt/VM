#include <math.h>
#include "vect.h"
#include <iostream>



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
		double P = (parametrs.a0 + parametrs.a1 * real_point[i].x);
		double prev_Q = (1 + prev_parametrs.b1 * real_point[i].x);
		result += ((real_point[i].y * Q - P) / prev_Q) * ((real_point[i].y * Q - P) / prev_Q);
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
	vect<double, double>* real_point,  // точки по которым аппроксимируем (первая координата - параметр, вторая - значения корня при этом параметре)
	int length,	// длина массива real_point
	vect3D<double, double, double> parametrs, // точка в которой ищется градиент
	vect3D<double, double, double> prev_parametrs
)
{
	const double delta = 0.0002;
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
	const double delta = 0.0002;
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
	const double delta = 0.0002;
	return ((p_func(real_point, length, parametrs, prev_parametrs, grad, lambda + delta) + p_func(real_point, length, parametrs, prev_parametrs, grad, lambda - delta) - 2 * p_func(real_point, length, parametrs, prev_parametrs, grad, lambda)) / (delta * delta));
}

double parabolaMin(
	double (*p_func_frst_derv)(double (*p_func)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double), vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double),
	double (*p_func_secnd_derv)(double (*p_func)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double), vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double),
	double (*p_func)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double),
	vect<double, double>* real_point,  // точки по которым аппроксимируем (первая координата - параметр, вторая - значения корня при этом параметре)
	int length,	// длина массива real_point
	vect3D<double, double, double> parametrs,
	vect3D<double, double, double> prev_parametrs,
	vect3D<double, double, double> grad,
	double lambda,
	double prev_lambda)
{
	double scnd_derv  =  p_func_secnd_derv(p_func,real_point, length, parametrs, prev_parametrs, grad, lambda);

	return(prev_lambda - p_func_frst_derv(p_func ,real_point,length,parametrs,prev_parametrs,grad,lambda) / (scnd_derv*sgn(scnd_derv)) );
}

vect3D<double, double, double> approximationOfRoot(vect<double, double>* real_point, int length) {
	double (*p_diff_aprox)(vect<double, double>*, int, vect3D<double,double,double>, vect3D<double, double, double>) = differenceOfAproximation; // указатель на функцию разности апроксимации и реальных точек
	double (*p_diff_aprox_1D)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>,double) = differenceOfAproximation1D; // указатель на функции разности апроксимации и реальных точек через градиент и лямбда
	double (*p_frst_derv_1D)(double (*p_func)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double),vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double) = derivative1D;
	double (*p_scnd_derv_1D)(double (*p_func)(vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double),vect<double, double>*, int, vect3D<double, double, double>, vect3D<double, double, double>, vect3D<double, double, double>, double) = secondDerivative1D;
	double const h =0.0001;  // а есть логика выбора h?
	vect3D<double, double, double> prev_r = {0.2,0.2,0.2};
	vect3D<double, double, double> prev_r_1 = { 0,0,0 };
	vect3D<double, double, double> prev_r_2 = { 0,0,0 };
	vect3D<double, double,double> r_1 = { 0,0,0 };
	vect3D<double, double,double> r_2 = { 0,0,0 };
	// код прямо копирка из второго задания с учётом многомерности, логику не менял!
	double prev_lambda = 0;
	double lambda = 0;
	vect3D<double, double,double> grad = gradFunc(p_diff_aprox,real_point,length,prev_r,prev_r)/modulVect3D(gradFunc(p_diff_aprox, real_point, length, prev_r, prev_r));
	do
	{
		prev_lambda = lambda;
		lambda = parabolaMin(p_frst_derv_1D,p_scnd_derv_1D, p_diff_aprox_1D,real_point,length,r_1, prev_r,grad, lambda, prev_lambda);
	} while (p_diff_aprox_1D(real_point,length,prev_r, prev_r, grad, prev_lambda) > p_diff_aprox_1D(real_point, length, prev_r,prev_r, grad, lambda));
	prev_r_2 = r_2;
	r_2 = prev_r + prev_lambda * grad;
	prev_r = 3 * h * prev_r;

	do
	{
		prev_r_1 = prev_r_2;
		r_1 = r_2;
		prev_lambda = 0;
		lambda = 0;
		grad = gradFunc(p_diff_aprox, real_point, length, prev_r, prev_r) / modulVect3D(gradFunc(p_diff_aprox, real_point, length, prev_r, prev_r));
		do
		{
			prev_lambda = lambda;
			lambda = parabolaMin(p_frst_derv_1D, p_scnd_derv_1D, p_diff_aprox_1D, real_point, length, r_1, prev_r, grad, lambda, prev_lambda);
		} while (p_diff_aprox_1D(real_point, length, prev_r, prev_r, grad, prev_lambda) > p_diff_aprox_1D(real_point, length, prev_r, prev_r, grad, lambda));
		prev_r_2 = r_2;
		r_2 = prev_r + prev_lambda * grad;
		prev_r = r_2 + h * (r_2 - r_1) * sgn(p_diff_aprox(real_point,length,r_1,prev_r_1) - p_diff_aprox(real_point, length, r_2, prev_r_2));
	} while ((p_diff_aprox(real_point, length, r_1, prev_r_1) > p_diff_aprox(real_point, length, r_2, prev_r_2)));
	// кстати ctrl+c & ctrl+v из второй задачки, я тут подумал что модуль второй производной для минимума пораболы используется же только в одномере, так что нормальный вариант!
	std::cout << p_diff_aprox(real_point, length, r_1, prev_r_1)<< '\n';
	return{ r_1};
}