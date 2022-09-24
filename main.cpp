//version 0.7.5 pre-release
#include <iostream>
#include <ctime>
#include <chrono>
#include "numbertheory.h"
#include "highprecision.h"
#include "optimization.h"
using namespace std;
using namespace nbt;
using namespace hpc;
using namespace opt;

int main()
{
    /*BigInt a("986732000000000000000123567");
    BigInt m("5000000009");
    BigInt n("10000000019");
    BigInt x = Fibonacci(1000000);
    BigInt y = Fibonacci(100000);
        
    auto start = std::chrono::system_clock::now();
    std::cout << log2(a) << std::endl;
    auto end = std::chrono::system_clock::now();
    std::cout << "total time=" << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
    */
    function<double(Matrix<double>)> f = [](Matrix<double> input)
    {
        return input(0, 0) * input(0, 0) + 10. * input(1, 0) * input(1, 0);
    };
    double epsilon = 1e-2;
    Matrix<double> initial1(vector<vector<double>>({{10.}, {1.}}));
    Matrix<double> initial2(vector<vector<double>>({{-10.}, {-1.}}));
    double step = 0.085;
    std::cout << Transpose(GradientDescent(f, initial1, step, epsilon)) << std::endl;
    std::cout << Transpose(BarzilarBorwein(f, initial2, 0.8, 20, 0.5,0.5,epsilon)) << std::endl;
    system("pause");
    return 0;
}
                                    
                    
                                    