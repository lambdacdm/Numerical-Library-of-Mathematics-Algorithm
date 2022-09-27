#include <iostream>
#include "../lib/mathalgorithm.h"
using namespace mal;

/* minimal inclusion:
#include "../lib/ode.h"
using namespace ode; */

int main()
{
    function<double(double, double)> f = [](double x, double y) { return y; };
    std::cout << DSolve(f, {0., 1.}, 1.) << std::endl;
    system("pause");
    return 0;
}