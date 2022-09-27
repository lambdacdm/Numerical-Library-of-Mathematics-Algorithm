#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H
#include "linearalgebra.h"
using namespace lag;

// 参考资料：
// 理论：刘浩洋、户将、李勇锋、文再文《最优化：建模、算法与理论》
// 实践：https://zhuanlan.zhihu.com/p/454059326

namespace opt{

//-----声明部分-----

//无约束优化
template <class OP>
Matrix<OP> GradientDescentBacktracking(const function<OP(Matrix<OP>)> &f,
                                                 const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> GradientDescentGoldstein(const function<OP(Matrix<OP>)> &f,
                                    const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> GradientDescentWolfe(const function<OP(Matrix<OP>)> &f,
                                const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> BarzilarBorwein(const function<OP(Matrix<OP>)> &f,
                           const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> FletcherReeves(const function<OP(Matrix<OP>)> &f,
                  const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> PolakRibiere(const function<OP(Matrix<OP>)> &f,
                const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> ClassicNewton(const function<OP(Matrix<OP>)> &f,
                         const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> ModifiedNewtonBacktracking(const function<OP(Matrix<OP>)> &f,
                                                const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> ModifiedNewtonGoldstein(const function<OP(Matrix<OP>)> &f,
                                   const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> ModifiedNewtonWolfe(const function<OP(Matrix<OP>)> &f,
                               const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> SR1(const function<OP(Matrix<OP>)> &f,
               const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> DFP(const function<OP(Matrix<OP>)> &f,
               const Matrix<OP> &initial_x, OP epsilon);
template <class OP>
Matrix<OP> BFGS(const function<OP(Matrix<OP>)> &f,
                const Matrix<OP> &initial_x, OP epsilon);

//约束优化
template <class OP>
Matrix<OP> AugmentedLagrangian(const function<OP(Matrix<OP>)> &f,
                               const function<Matrix<OP>(Matrix<OP>)> &c, const Matrix<OP> &initial_x, OP epsilon,
                               const function<Matrix<OP>(const function<OP(Matrix<OP>)> &, const Matrix<OP> &, OP)> &argmin);

//-----定义部分-----

//求导算子
template<class OP> Matrix<OP> Gradient(const function<OP(Matrix<OP>)> &f, const Matrix<OP> &x)
{
    const OP h=1e-8;
    const int n = RowSize(x);
    Matrix<OP> df(n, 1);
    for (int i = 0; i < n;++i)
    {
        Matrix<OP> x1 = x;
        Matrix<OP> x2 = x;
        x1(i, 0) = x(i, 0) - h;
        x2(i, 0) = x(i, 0) + h;
        df(i, 0) = (f(x2) - f(x1)) / (2*h);
    }
    return df;
}
template<class OP> Matrix<OP> GradientT(const function<OP(Matrix<OP>)> &f, const Matrix<OP> &x)
{
    const OP h=1e-8;
    const OP twoh = 2 * h;
    const int n = RowSize(x);
    Matrix<OP> df(1, n);
    for (int i = 0; i < n;++i)
    {
        Matrix<OP> x1 = x;
        Matrix<OP> x2 = x;
        x1(i, 0) = x(i, 0) - h;
        x2(i, 0) = x(i, 0) + h;
        df(0, i) = (f(x2) - f(x1)) / twoh;
    }
    return df;
}
template<class OP> Matrix<OP> Hessian(const function<OP(Matrix<OP>)> &f, const Matrix<OP> &x)
{
    const OP h=1e-4;
    const OP twoh = 2 * h;
    const int n = RowSize(x);
    Matrix<OP> H(n, n);
    for (int i = 0; i < n;++i)
    {
        Matrix<OP> x1 = x;
        Matrix<OP> x2 = x;
        x1(i, 0) = x(i, 0) + twoh;
        x2(i, 0) = x(i, 0) - twoh;
        H(i, i) = ((f(x1) - f(x)) / twoh + (f(x2) - f(x)) / twoh) / twoh;
        for (int j = i + 1; j < n;++j)
        {
            Matrix<OP> x1 = x;
            Matrix<OP> x2 = x;
            Matrix<OP> x3 = x;
            Matrix<OP> x4 = x;
            x1(i, 0) = x(i, 0) + h; x1(j, 0) = x(j, 0) + h;
            x2(i, 0) = x1(i,0); x2(j, 0) = x(j, 0) - h;
            x3(i, 0) = x(i, 0) - h; x3(j, 0) = x1(j, 0);
            x4(i, 0) = x3(i,0); x4(j, 0) = x2(j, 0);
            H(i, j) = ((f(x1) - f(x2)) / twoh + (f(x4) - f(x3)) / twoh) / twoh;
        }
    }
    for (int i = 1; i < n;++i)
        for (int j = 0; j < i;++j)
            H(i, j) = H(j, i);
    return H;
}

//线搜索方法
template<class OP> Matrix<OP> OneStepLineSearch(const Matrix<OP> &x, 
const Matrix<OP> &direction, OP stepsize)
{
    return x + (stepsize * direction);
}
template<class OP> OP _BacktrackingStepsize(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &x, const Matrix<OP> &direction, OP stepsize, OP gamma, OP c)
{
    const uint8_t times = 10;
    uint8_t cnt = 0;
    while(f(x + (stepsize * direction)) > f(x)+c*stepsize*(GradientT<OP>(f,x)*direction)(0,0)
    && cnt<times)
    {
        stepsize = gamma * stepsize;
        ++cnt;
    }
    return stepsize;
}
template<class OP> OP BacktrackingStepsize(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &x, const Matrix<OP> &direction)
{
    return _BacktrackingStepsize<OP>(f, x, direction, OP(1), OP(0.5), OP(1e-3));
}
template<class OP> OP _GoldsteinStepsize(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &x, const Matrix<OP> &direction, OP stepsize, OP c)
{
    const uint8_t times = 10;
    uint8_t cnt = 0;
    OP fx = f(x);
    OP fd = (GradientT<OP>(f, x) * direction)(0, 0);
    OP lhs = f(x + stepsize * direction);
    OP afd = stepsize * fd;
    OP left_end = OP(0);
    OP right_end = OP(4) * stepsize;
    while (cnt<times)
    {
        if (lhs > fx + c * afd )
        {
            right_end = stepsize;
            stepsize = (left_end + right_end) / OP(2);
        }
        else if (lhs < fx + (1-c)*afd)
        {
            left_end = stepsize;
            stepsize = std::min(OP(2) * stepsize, (left_end + right_end) / OP(2));
        }
        else
            break;
        lhs = f(x + stepsize * direction);
        afd = stepsize * fd;
        ++cnt;
    }
    return stepsize;
}
template<class OP> OP GoldsteinStepsize(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &x, const Matrix<OP> &direction)
{
    return _GoldsteinStepsize<OP>(f, x, direction, OP(1), OP(0.2));
}
template<class OP> OP _GeneralizedWolfeStepsize(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &x, const Matrix<OP> &direction, OP stepsize, OP c1, OP c2, bool strong)
{
    const uint8_t times = 10;
    uint8_t cnt = 0;
    OP fd = (GradientT<OP>(f, x) * direction)(0, 0);
    Matrix<OP> xad = x + stepsize * direction;
    OP left_end = OP(0);
    OP right_end = OP(4) * stepsize;
    while (cnt<times)
    {
        if (f(xad) > f(x) + c1 * stepsize * fd)
        {
            right_end = stepsize;
            stepsize = (left_end + right_end) / OP(2);
        }
        else if (!strong && Abs((GradientT<OP>(f, xad)*direction)(0,0)) < c2 * fd)
        {
            left_end = stepsize;
            stepsize = std::min(OP(2) * stepsize, (left_end + right_end) / OP(2));
        }
        else if (strong && Abs((GradientT<OP>(f, xad)*direction)(0,0)) > -c2 * fd)
        {
            left_end = stepsize;
            stepsize = std::min(OP(2) * stepsize, (left_end + right_end) / OP(2));
        }
        else
            break;
        xad = x + stepsize * direction;
        ++cnt;
    }
    return stepsize;
}
template<class OP> OP WolfeStepsize(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &x, const Matrix<OP> &direction)
{
    return _GeneralizedWolfeStepsize<OP>(f, x, direction, OP(1), OP(1e-3), OP(0.9), false);
}
template<class OP> OP StrongWolfeStepsize(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &x, const Matrix<OP> &direction)
{
    return _GeneralizedWolfeStepsize<OP>(f, x, direction, OP(1), OP(1e-3), OP(0.9), true);
}

//梯度类算法
template<class OP> Matrix<OP> OneStepGradientDescent(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &x, OP stepsize)
{
    return x - stepsize * Gradient(f, x);
}
template<class OP> Matrix<OP> GradientDescentFixedStepsize(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon, OP stepsize)
{
    const uint8_t times = 200;
    uint8_t cnt = 0;
    Matrix<OP> x = initial_x;
    while(Norm(Gradient(f,x), "Frobenius")>epsilon && cnt<times)
    {
        x = OneStepGradientDescent<OP>(f, x, stepsize);
        ++cnt;
    }
    if(cnt==times)
    {
        std::cerr << "来自文件optimization.h中的函数opt::GradientDescentFixedStepsize："
        <<std::endl<<"警告：已到迭代次数上限。未达到指定精度。" << std::endl;
    }
    return x;
}
template<class OP> Matrix<OP> GradientDescent(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon,
const function<OP(const function<OP(Matrix<OP>)> &,const Matrix<OP> &,const Matrix<OP> &)>& update_stepsize)
{
    const uint8_t times = 200;
    uint8_t cnt = 0;
    Matrix<OP> x = initial_x;
    while(Norm(Gradient(f,x), "Frobenius")>epsilon && cnt<times)
    {
        Matrix<OP> direction = - Gradient(f, x);
        OP stepsize = update_stepsize(f, x, direction);
        x = OneStepLineSearch(x, direction, stepsize);
        ++cnt;
    }
    if(cnt==times)
    {
        std::cerr << "来自文件optimization.h中的函数opt::GradientDescent："
        <<std::endl<<"警告：已到迭代次数上限。未达到指定精度。" << std::endl;
    }
    return x;
}
template<class OP> Matrix<OP> GradientDescentBacktracking(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon)
{
    return GradientDescent<OP>(f, initial_x, epsilon, BacktrackingStepsize<OP>);
}
template<class OP> Matrix<OP> GradientDescentGoldstein(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon)
{
    return GradientDescent<OP>(f, initial_x, epsilon, GoldsteinStepsize<OP>);
}
template<class OP> Matrix<OP> GradientDescentWolfe(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon)
{
    return GradientDescent<OP>(f, initial_x, epsilon, WolfeStepsize<OP>);
}
template<class OP> Matrix<OP> BarzilarBorwein(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP initial_step, uint32_t M, OP beta, OP c, OP epsilon)
{
    const uint32_t times = 200;
    uint32_t cnt = 0;

    Matrix<OP> x = initial_x;
    vector<Matrix<OP>> x_list;
    x_list.push_back(x);

    Matrix<OP> gradient = Gradient(f, x);
    vector<Matrix<OP>> gradient_list;
    gradient_list.push_back(gradient);

    OP stepsize = initial_step;
    vector<OP> stepsize_list;
    stepsize_list.push_back(stepsize);
    
    OP norm_of_gradient = Norm<OP>(gradient, "Frobenius");

    while(norm_of_gradient>epsilon && cnt<times)
    {
        while(true)
        {
            bool flag = true;
            OP lhs = f(x - stepsize * gradient);
            OP rhs_factor = c * norm_of_gradient * norm_of_gradient;
            for (uint32_t j = 0; j <= std::min(cnt, M);++j)
            {
                if(lhs < f(x_list[cnt-j])- rhs_factor * stepsize)
                {
                    flag = false;
                    break;
                }
            }
            if(flag)
                stepsize = beta * stepsize;
            else
                break;
        }
        ++cnt;

        x = x - stepsize * gradient;
        x_list.push_back(x);

        gradient = Gradient(f, x);
        gradient_list.push_back(gradient);

        Matrix<OP> y = gradient_list[cnt] - gradient_list[cnt - 1];
        stepsize = InnerProductC<OP>(x_list[cnt] - x_list[cnt - 1], y) / InnerProductC<OP>(y, y);
        stepsize_list.push_back(stepsize);

        norm_of_gradient = Norm<OP>(gradient, "Frobenius");
    }
    if(cnt==times)
    {
        std::cerr << "来自文件optimization.h中的函数opt::BarzilarBorwein："
        <<std::endl<<"警告：已到迭代次数上限。未达到指定精度。" << std::endl;
    }
    return x;
}
template<class OP> Matrix<OP> BarzilarBorwein(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon)
{
    return BarzilarBorwein<OP>(f, initial_x, OP(1), 20u, OP(0.5), OP(1e-3), epsilon);
}

//共轭梯度法
template <class OP> Matrix<OP> ConjugateGradient(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon, OP restart_threshold,
const function<OP(const Matrix<OP> &, const Matrix<OP> &)> &update_beta)
{
    const uint8_t times = 200;
    uint8_t cnt = 0;
    Matrix<OP> x = initial_x;
    Matrix<OP> gradient = Gradient(f, x);
    Matrix<OP> direction = -gradient;
    Matrix<OP> x_pre = x;
    Matrix<OP> gradient_pre = gradient;
    while(Norm<OP>(gradient, "Frobenius")>epsilon && cnt<times)
    {
        x = OneStepLineSearch(x, direction, StrongWolfeStepsize<OP>(f, x, direction));
        gradient = Gradient(f, x);
        if(Abs(InnerProductC(gradient,gradient_pre))/InnerProductC(gradient,gradient)<restart_threshold)
            direction = -gradient;
        else
            direction = update_beta(gradient, gradient_pre) * direction - gradient;
        x_pre = x;
        gradient_pre = gradient;
        ++cnt;
    }
    if(cnt==times)
    {
        std::cerr << "来自文件optimization.h中的函数opt::ConjugateGradient："
        <<std::endl<<"警告：已到迭代次数上限。未达到指定精度。" << std::endl;
    }
    return x;
}
template<class OP> OP FletcherReevesUpdate(const Matrix<OP> &gradient,const Matrix<OP> &gradient_pre)
{
    return InnerProductC<OP>(gradient, gradient) / InnerProductC<OP>(gradient_pre, gradient_pre);
}
template<class OP> Matrix<OP> FletcherReeves(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &initial_x, OP epsilon)
{
    return ConjugateGradient(f, initial_x, epsilon, OP(0.1), FletcherReevesUpdate<OP>);
}
template<class OP> OP PolakRibiereUpdate(const Matrix<OP> &gradient,const Matrix<OP> &gradient_pre)
{
    return std::max(InnerProductC<OP>(gradient, gradient - gradient_pre) / 
        InnerProductC<OP>(gradient_pre, gradient_pre) , OP(0));
}
template<class OP> Matrix<OP> PolakRibiere(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &initial_x, OP epsilon)
{
    return ConjugateGradient(f, initial_x, epsilon, OP(0.1), PolakRibiereUpdate<OP>);
}

//牛顿类算法
template <class OP> Matrix<OP> OneStepClassicNewton(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &x)
{
    Matrix<OP> gradient = Gradient<OP>(f, x);
    Matrix<OP> hessian = gradient * Transpose<OP>(gradient);
    return x - Inverse<OP>(hessian) * gradient;
} 
template <class OP> Matrix<OP> ClassicNewton(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &initial_x, OP epsilon)
{
    const uint8_t times = 200;
    uint8_t cnt = 0;
    Matrix<OP> x = initial_x;
    Matrix<OP> gradient = Gradient<OP>(f, x);
    Matrix<OP> hessian = Hessian<OP>(f, x);
    while(Norm(gradient, "Frobenius")>epsilon && cnt<times)
    {
        x = x - Inverse<OP>(hessian) * gradient;
        gradient = Gradient<OP>(f, x);
        hessian = Hessian<OP>(f, x);
        ++cnt;
    }
    if(cnt==times)
    {
        std::cerr << "来自文件optimization.h中的函数opt::ClassicNewton："
        <<std::endl<<"警告：已到迭代次数上限。未达到指定精度。" << std::endl;
    }
    return x;
}
template <class OP> Matrix<OP> ModifiedNewton(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &initial_x, OP epsilon,
const function<OP(const function<OP(Matrix<OP>)> &,const Matrix<OP> &,const Matrix<OP> &)> &update_stepsize)
{
    const uint8_t times = 200;
    uint8_t cnt = 0;
    Matrix<OP> x = initial_x;
    Matrix<OP> gradient = Gradient<OP>(f, x);
    while(Norm(gradient, "Frobenius")>epsilon && cnt<times)
    {
        Matrix<OP> B = Hessian<OP>(f, x);
        vector<OP> eigenvalues = Eigenvalues<OP>(B);
        OP minimal_eigenvalue = *std::min_element(eigenvalues.begin(), eigenvalues.end());
        if(minimal_eigenvalue < OP(0.01))
            B = B + (OP(0.01) - minimal_eigenvalue) * Eye<OP>(RowSize(B));
        Matrix<OP> direction = LinearSolve<OP>(B, -gradient, "Cholesky");
        OP stepsize = update_stepsize(f, x, direction);
        x = OneStepLineSearch<OP>(x, direction, stepsize);
        gradient = Gradient<OP>(f, x);
        ++cnt;
    }
    if(cnt==times)
    {
        std::cerr << "来自文件optimization.h中的函数opt::ModifiedNewton："
        <<std::endl<<"警告：已到迭代次数上限。未达到指定精度。" << std::endl;
    }
    return x;
}
template<class OP> Matrix<OP> ModifiedNewtonBacktracking(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon)
{
    return ModifiedNewton<OP>(f, initial_x, epsilon, BacktrackingStepsize<OP>);
}
template<class OP> Matrix<OP> ModifiedNewtonGoldstein(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon)
{
    return ModifiedNewton<OP>(f, initial_x, epsilon, GoldsteinStepsize<OP>);
}
template<class OP> Matrix<OP> ModifiedNewtonWolfe(const function<OP(Matrix<OP>)> &f,
const Matrix<OP> &initial_x, OP epsilon)
{
    return ModifiedNewton<OP>(f, initial_x, epsilon, WolfeStepsize<OP>);
}

//拟牛顿类算法
template <class OP> Matrix<OP> QuasiNewton(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &initial_x, const Matrix<OP> &initial_H, OP epsilon,
const function<Matrix<OP>(const Matrix<OP> &,const Matrix<OP> &,const Matrix<OP> &,
const Matrix<OP> &,const Matrix<OP> &)> &update_H)
{
    const uint8_t times = 200;
    uint8_t cnt = 0;
    Matrix<OP> x = initial_x;
    Matrix<OP> H = initial_H;
    Matrix<OP> gradient = Gradient<OP>(f, x);
    while (Norm(gradient, "Frobenius")>epsilon && cnt<times)
    {
        Matrix<OP> x_pre = x;
        Matrix<OP> gradient_pre = gradient;

        Matrix<OP> direction = - (H * gradient);
        OP stepsize = WolfeStepsize<OP>(f, x, direction);

        x = OneStepLineSearch(x, direction, stepsize);
        gradient = Gradient<OP>(f, x);

        H = update_H(H, x, x_pre, gradient, gradient_pre);
        ++cnt;
    }
    if(cnt==times)
    {
        std::cerr << "来自文件optimization.h中的函数opt::QuasiNewton："
        <<std::endl<<"警告：已到迭代次数上限。未达到指定精度。" << std::endl;
    }
    return x;
}
template <class OP> Matrix<OP> QuasiNewton(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &initial_x, OP epsilon,
const function<Matrix<OP>(const Matrix<OP> &,const Matrix<OP> &,const Matrix<OP> &,
const Matrix<OP> &,const Matrix<OP> &)> &update)
{
    return QuasiNewton<OP>(f, initial_x, Eye<OP>(RowSize(initial_x)), epsilon, update);
}
template <class OP> Matrix<OP> SR1Update(const Matrix<OP> &H, const Matrix<OP> &x, 
const Matrix<OP> &x_pre, const Matrix<OP> &gradient, const Matrix<OP> &gradient_pre)
{
    Matrix<OP> s = x - x_pre;
    Matrix<OP> y = gradient - gradient_pre;
    Matrix<OP> M = s - H * y;
    Matrix<OP> MT = Transpose(M);
    return H + (OP(1) / (MT * y)(0, 0)) * (M * MT);
}
template<class OP> Matrix<OP> SR1(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &initial_x, OP epsilon)
{
    return QuasiNewton(f, initial_x, epsilon, SR1Update<OP>);
}
template <class OP> Matrix<OP> DFPUpdate(const Matrix<OP> &H, const Matrix<OP> &x, 
const Matrix<OP> &x_pre, const Matrix<OP> &gradient, const Matrix<OP> &gradient_pre)
{
    Matrix<OP> s = x - x_pre;
    Matrix<OP> y = gradient - gradient_pre;
    Matrix<OP> Hy = H * y;
    Matrix<OP> yT = Transpose(y);
    return H - (OP(1) / (yT * Hy)(0, 0)) * Hy * Transpose(Hy) + 
                (OP(1) / (yT * s)(0, 0)) * s * Transpose(s);
}
template<class OP> Matrix<OP> DFP(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &initial_x, OP epsilon)
{
    return QuasiNewton(f, initial_x, epsilon, DFPUpdate<OP>);
}
template <class OP> Matrix<OP> BFGSUpdate(const Matrix<OP> &H, const Matrix<OP> &x, 
const Matrix<OP> &x_pre, const Matrix<OP> &gradient, const Matrix<OP> &gradient_pre)
{
    Matrix<OP> s = x - x_pre;
    Matrix<OP> y = gradient - gradient_pre;
    OP rho = OP(1) / InnerProductC(s, y);
    Matrix<OP> rho_s = rho * s;
    Matrix<OP> M = Eye<OP>(RowSize(H)) - rho_s * Transpose(y);
    return Transpose(M) * H * M + rho_s * Transpose(s);
}
template<class OP> Matrix<OP> BFGS(const function<OP(Matrix<OP>)> &f, 
const Matrix<OP> &initial_x, OP epsilon)
{
    return QuasiNewton(f, initial_x, epsilon, BFGSUpdate<OP>);
}

//增广拉格朗日函数法
template<class OP> function<OP(Matrix<OP>)> AugmentedLagrangianFcuntion(
const function<OP(Matrix<OP>)> &f, const function<Matrix<OP>(Matrix<OP>)> &c,
Matrix<OP> multiplier, OP penalty_coeff)
{
    return [&](Matrix<OP> x)
    {
        Matrix<OP> cx = c(x);
        return f(x) + InnerProductC<OP>(multiplier, cx) +
               (penalty_coeff / OP(2)) * InnerProductC<OP>(cx, cx);
    };
}
template<class OP> Matrix<OP> AugmentedLagrangian(const function<OP(Matrix<OP>)> &f, 
const function<Matrix<OP>(Matrix<OP>)> &c, const Matrix<OP> &initial_x, OP epsilon,
const function<Matrix<OP>(const function<OP(Matrix<OP>)> &, const Matrix<OP> &, OP)> &argmin,
Matrix<OP> multiplier, OP penalty_coeff, OP penalty_update)
{
    const uint8_t times = 64;
    uint8_t cnt = 0;
    Matrix<OP> x = initial_x;
    Matrix<OP> cx = c(x);
    OP eta = OP(1);
    while(Norm<OP>(cx, "Frobenius")>epsilon && cnt<times)
    {
        x = argmin(AugmentedLagrangianFcuntion<OP>(f,c,multiplier,penalty_coeff),x,eta);
        eta = eta / OP(2);
        cx = c(x);
        multiplier = multiplier + penalty_coeff * cx;
        penalty_coeff = penalty_update * penalty_coeff;
        ++cnt;
    }
    if(cnt==times)
    {
        std::cerr << "来自文件optimization.h中的函数opt::AugmentedLagrangian："
        <<std::endl<<"警告：已到迭代次数上限。未达到指定精度。" << std::endl;
    }
    return x;
}  
template<class OP> Matrix<OP> AugmentedLagrangian(const function<OP(Matrix<OP>)> &f, 
const function<Matrix<OP>(Matrix<OP>)> &c, const Matrix<OP> &initial_x, OP epsilon,
const function<Matrix<OP>(const function<OP(Matrix<OP>)> &, const Matrix<OP> &, OP)> &argmin)
{
    return AugmentedLagrangian(f, c, initial_x, epsilon, argmin,
                               Matrix<OP>(OP(1), RowSize(c(initial_x)), 1), OP(1), OP(2));
}

}                                                                                              
#endif                                                                              