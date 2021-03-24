#include "min2D.h"
#include "root.h"
#include <iostream>

int main() {
	std::cout << "Hello User! I search root for a/x - exp(x).\n" << "Type a: ";
	std::cout << "I found "<< root() << '\n';
	std::cout.precision(6);
	std::cout << "Hello User! I search minimum of (xy - ctg(x))^2 + sin^2(xy).\n";
	vect<double, double> r1 = min2D({ 1,1 });
	std::cout << "I found min of (xy - ctg(x))^2 + sin^2(xy) in " << r1.x << " " << r1.y << '\n';
	vect<double, double> r2 = min2D({ -1,-1 });
	std::cout << "I found min of (xy - ctg(x))^2 + sin^2(xy) in " << r2.x << " " << r2.y << '\n';
	std::cout << "I decided that ";
	if (f(r1) < f(r2)) std::cout << "min of (xy - ctg(x))^2 + sin^2(xy) in " << r1.x << " " << r1.y << "\n" << " Value of func " << f(r1);
	else std::cout << "min of (xy - ctg(x))^2 + sin^2(xy) in " << r2.x << " " << r2.y << "\n" << " Value of func " << f(r2);

	return 1;
}