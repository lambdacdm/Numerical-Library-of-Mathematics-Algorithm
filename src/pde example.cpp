#include <iostream>
#include "../lib/mathalgorithm.h"
using namespace mal;

/* minimal inclusion:
#include "../lib/pde.h"
using namespace pde; */

int main()
{
    function<double(double)> f = [](double x) { return sin(4*Pi<double>*x); };
    function<double(double)> ze = [](double x) { return 0; };
    std::cout << Transpose(HeatDirichletFDMExplicit<double>(1,f,ze,ze,0,1, 64, 320, 0.03))
              << std::endl;
    system("pause");
    return 0;
}