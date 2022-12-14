#include <iostream>
#include "../lib/mathalgorithm.h"
using namespace mal;

/* minimal inclusion:
#include "../lib/optimization.h"
using namespace opt; */

int main()
{
    double epsilon = 1e-6;

    std::cout << "约束优化问题：" << std::endl;
    function<double(double,double)> f_pre = [](double x, double y)
    {
        return x + sqrt(3.) * y;
    };
    auto f = MatrizeInputs(f_pre);
    function<Matrix<double>(double,double)> c_pre = [](double x, double y)
    {
        return Matrix<double>(x*x+y*y-1.,1,1);
    };
    auto c = MatrizeInputs(c_pre);
    Matrix<double> initial(1.,2,1);
    std::cout << "AugmentedLagrangian"<<std::endl<<
    Transpose(AugmentedLagrangian<double>(f, c, initial, epsilon, BFGS<double>))
    <<std::endl;

    std::cout << "无约束优化问题：" << std::endl; 
    function<double(double,double)> g_pre = [](double x,double y)
    {
        return x * x + 10. * y * y;
    };
    auto g = MatrizeInputs(g_pre);
    Matrix<double> initial2(vector<vector<double>>({{-10.}, {-1.}}));
    double step = 0.085;
    std::cout << "GradientDescentFixedStepsize"<<std::endl<<
    Transpose(GradientDescentFixedStepsize(g, initial2, epsilon, step)) << std::endl;
    std::cout << "GradientDescentWolfe"<<std::endl<<
    Transpose(GradientDescentWolfe(g, initial2, epsilon)) << std::endl;
    std::cout << "GradientDescentGoldstein"<<std::endl<<
    Transpose(GradientDescentGoldstein(g, initial2, epsilon)) << std::endl;
    std::cout << "BarzilarBorwein"<<std::endl<<
    Transpose(BarzilarBorwein(g, initial2, epsilon)) << std::endl;
    std::cout << "FletcherReeves" << std::endl<<
    Transpose(FletcherReeves(g, initial2, epsilon)) << std::endl;
    std::cout << "PolakRibiere" << std::endl<<
    Transpose(PolakRibiere(g, initial2, epsilon)) << std::endl;
    std::cout << "SR1"<<std::endl<<
    Transpose(SR1(g, initial2, epsilon)) << std::endl;
    std::cout << "DFP"<<std::endl<<
    Transpose(DFP(g, initial2, epsilon)) << std::endl;
    std::cout << "BFGS"<<std::endl<<
    Transpose(BFGS(g, initial2, epsilon)) << std::endl;

    Matrix<double> initial3(vector<vector<double>>({{-0.1}, {-0.1}}));
    std::cout << "ClassicNewton"<<std::endl<<
    Transpose(ClassicNewton(g, initial3, epsilon)) << std::endl;
    std::cout << "ModifiedNewtonGoldstein"<<std::endl<<
    Transpose(ModifiedNewtonGoldstein(g, initial2, epsilon)) << std::endl;
    system("pause");
    return 0;
}              
                                    