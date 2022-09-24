#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H
#include "numerical.h"
using namespace nmr;

namespace opt{

//-----声明部分-----

//-----定义部分-----

//梯度
template<class OP> Matrix<OP> Gradient(const function<OP(Matrix<OP>)> &f, const Matrix<OP> &x)
{
    const OP h=1e-7;
    const int n = RowSize(x);
    Matrix<OP> df(n, 1);
    for (int i = 0; i < n;++i)
    {
        Matrix<OP> x1 = x;
        Matrix<OP> x2 = x;
        x1(i, 0) = Get(x, i, 0) - h;
        x2(i, 0) = Get(x, i, 0) + h;
        df(i, 0) = (f(x2) - f(x1)) / (2*h);
    }
    return df;
}
template<class OP> Matrix<OP> GradientT(const function<OP(Matrix<OP>)> &f, const Matrix<OP> &x)
{
    const OP h=1e-7;
    const int n = RowSize(x);
    Matrix<OP> df(1, n);
    for (int i = 0; i < n;++i)
    {
        Matrix<OP> x1 = x;
        Matrix<OP> x2 = x;
        x1(i, 0) = Get(x, i, 0) - h;
        x2(i, 0) = Get(x, i, 0) + h;
        df(0, i) = (f(x2) - f(x1)) / (2*h);
    }
    return df;
}

//线搜索方法
template<class OP> Matrix<OP> OneStepLineSearch(const Matrix<OP> &x, 
const Matrix<OP> &direction, OP step)
{
    return x + (step * direction);
}
template<class OP> OP BacktrackingLineSearchStepsize(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &x, const Matrix<OP> &direction, OP gamma, OP c, OP initial_step)
{
    OP step = initial_step;
    while(f(x + (step * direction)) > f(x)+c*x*(GradientT(f,x)*direction))
    {
        step = gamma * step;
    }
    return step;
}

//梯度类算法
template<class OP> Matrix<OP> OneStepGradientDescent(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &x, OP step)
{
    return x - step * Gradient(f, x);
}
template<class OP> Matrix<OP> GradientDescent(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP step, OP epsilon)
{
    Matrix<OP> x = initial_x;
    while(Norm(Gradient(f,x), "Frobenius")>epsilon)
    {
        x = OneStepGradientDescent(f, x, step);
    }
    return x;
}
template<class OP> Matrix<OP> BarzilarBorwein(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP initial_step, uint32_t M, OP c, OP beta, OP epsilon)
{
    uint32_t cnt = 0;

    Matrix<OP> x = initial_x;
    vector<Matrix<OP>> x_list;
    x_list.push_back(x);

    Matrix<OP> gradient = Gradient(f, x);
    vector<Matrix<OP>> gradient_list;
    gradient_list.push_back(gradient);

    OP step = initial_step;
    vector<OP> step_list;
    step_list.push_back(step);
    
    OP norm_of_gradient = Norm(gradient, "Frobenius");

    while(norm_of_gradient>epsilon)
    {
        while(true)
        {
            bool flag = true;
            OP lhs = f(x - step * gradient);
            OP rhs_factor = c * norm_of_gradient * norm_of_gradient;
            for (uint32_t j = 0; j <= std::min(cnt, M);++j)
            {
                if(lhs < f(x_list[cnt-j])- rhs_factor * step)
                {
                    flag = false;
                    break;
                }
            }
            if(flag)
                step = beta * step;
            else
                break;
        }
        ++cnt;

        x = x - step * gradient;
        x_list.push_back(x);

        gradient = Gradient(f, x);
        gradient_list.push_back(gradient);

        Matrix<OP> y = gradient_list[cnt] - gradient_list[cnt - 1];
        step = InnerProductC<OP>(x_list[cnt] - x_list[cnt - 1], y) / InnerProductC<OP>(y, y);
        step_list.push_back(step);

        norm_of_gradient = Norm(gradient, "Frobenius");
    }
    return x;
}

}
#endif