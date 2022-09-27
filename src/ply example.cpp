#include <iostream>
#include "../lib/mathalgorithm.h"
using namespace mal;

/* minimal inclusion:
#include "../lib/polynomial.h"
using namespace ply; */

int main()
{
    Polynomial<int> f({1,1,1},2);
    Polynomial<int> g({5,-3,0,1},3);
    Polynomial<int> h=f*g;
    std::cout<<"h(x)="<<h<<std::endl;
    std::cout<<"h(5)="<<h(5)<<std::endl;
    system("pause");
    return 0;
}
