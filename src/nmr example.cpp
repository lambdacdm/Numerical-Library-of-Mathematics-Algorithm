#include <iostream>
#include "../lib/mathalgorithm.h"
using namespace mal;

/* minimal inclusion:
#include "../lib/numerical.h"
using namespace nmr; */

int main()
{
    Matrix<double> base(
    {{0.4,-0.916},{0.5,-0.693},{0.8,-0.223},{0.9,-0.105}});
    Matrix<double> data({{0.6},{0.7}},2,1);
    Matrix<double> result=Interpolation(base,data);
    std::cout<<result<<std::endl;

    auto p = [](double x) { return 4.0 / (1 + x * x); };
    auto q = Integrate<double>(p, 0, 1);
    N(q, 10);

    auto f=[](double x){ return x*x*exp(x); };
    auto g=D<double>(2,f);
    std::cout<<g(0)<<std::endl;

    system("pause");
    return 0;
}
