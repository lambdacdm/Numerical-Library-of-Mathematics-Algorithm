#include <iostream>
#include "../lib/mathalgorithm.h"
using namespace mal;

/* minimal inclusion:
#include "../lib/probability.h"
using namespace prb; */

int main()
{
    std::cout << CDF(NormalDistribution(1, 1), 1) << std::endl;
    system("pause");
    return 0;
} 
