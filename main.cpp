#include "vect.h"
#include "root.h"
#include "min2D.h"
#include "approximationOfRoot.h"
#include <iostream>
#include <fstream>

int main() {
	std::cout << "Hello User! I search root for a/x - exp(x).\n" << "Type a: ";
	std::cout << "I found "<< root(1) << '\n';
	std::cout.precision(6);
	std::cout << "Hello User! I search minimum of (xy - ctg(x))^2 + sin^2(xy).\n";
	vect<double, double> r1 = min2D({ 1,1 });
	std::cout << "I found min of (xy - ctg(x))^2 + sin^2(xy) in " << r1.x << " " << r1.y << '\n';
	vect<double, double> r2 = min2D({ -1,-1 });
	std::cout << "I found min of (xy - ctg(x))^2 + sin^2(xy) in " << r2.x << " " << r2.y << '\n';
	std::cout << "I decided that ";
	if (f(r1) < f(r2)) std::cout << "min of (xy - ctg(x))^2 + sin^2(xy) in " << r1.x << " " << r1.y << "\n" << " Value of func " << f(r1) << '\n';
	else std::cout << "min of (xy - ctg(x))^2 + sin^2(xy) in " << r2.x << " " << r2.y << "\n" << " Value of func " << f(r2) << '\n';
	// логика третьего задания:
	//
	std::ofstream file;
	file.open("output.txt");
	vect<double, double>* point;
	int length = 20;
	point = new vect<double, double>[length];
	for (int i = 0; i < length; i++)
	{
		point[i] = { 1 + 5*(double)i , root(1 + 5 * (double)i) };
		file << point[i].x << " " << point[i].y << '\n';
	}
	vect3D<double, double, double> aprox_param = approximationOfRoot(point, length);
	std::cout << "coeff of aproximation: a0 - " << aprox_param.a0 << ", a1 - " << aprox_param.a1 << ", b1 - " << aprox_param.b1 << ". That`s all!";
	file.close();
	return 1;
}