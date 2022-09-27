//version 0.7.5 pre-release
#include <iostream>
#include <ctime>
#include <chrono>
#include "mathalgorithm.h"
using namespace mal;

int main()
{
    function<double(double)> f = [](double x)
    { return sin(x); };
    Plot<double>(f, vector<double>({0., 2*Pi<double>}));
    system("pause");
    return 0;
}
                                    
                    
                                    