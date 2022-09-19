//version 0.7.6
#include <iostream>
#include <ctime>
#include "plot.h"
using namespace std;
using namespace plt;
int main()
{
    Plot<double>([](double t)
                         { return sin(t); }, {0.,2*Pi<double>});
    system("pause");
    return 0;
}
                                    
                    
                                    