#include <iostream>
#include "../lib/mathalgorithm.h"
using namespace mal;

/* minimal inclusion:
#include "../lib/plot.h"
using namespace plt; */

int main()
{
    function<double(double)> f = [](double x)
    { return sin(x); };
    Plot<double>(f, vector<double>({0., 2*Pi<double>}));
    system("pause");
    return 0;
}                  
                                    