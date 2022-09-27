#include <iostream>
#include "../lib/mathalgorithm.h"
using namespace mal;

/* minimal inclusion:
#include "../lib/linearalgebra.h"
using namespace lag; */

int main()
{
    Matrix<double> A({{1,1,-1},{1,2,-2},{-2,1,1}},3,3);
    Matrix<double> b({{1}, {0}, {1}}, 3, 1);
    std::cout << LinearSolve(A, b) <<std::endl;
    system("pause");
    return 0;
}
