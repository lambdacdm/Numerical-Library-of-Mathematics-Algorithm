#ifndef PDE_H
#define PDE_H
#include "numerical.h"
using namespace nmr;

// 参考资料：
// 理论： https://zhuanlan.zhihu.com/p/65247362
// 实践： https://zhuanlan.zhihu.com/p/550912102

namespace pde{

//-----声明部分-----
template <class PS> Matrix<PS> HeatDirichletFDMExplicit(PS, std::function<PS(PS)>,
std::function<PS(PS)>, std::function<PS(PS)>, PS, PS, int, int, PS);

//-----定义部分-----

template<class PS> Matrix<PS> HeatDirichletFDMExplicit(PS para, std::function<PS(PS)> initial_cond,
    std::function<PS(PS)> left_bond_cond, std::function<PS(PS)> right_bond_cond, 
    PS left_end, PS right_end, int x_divide, int t_divide, PS time)
{
    PS x_step = (right_end - left_end) / PS(x_divide);
    PS t_step = time / PS(t_divide);
    PS sigma = para * t_step / x_step / x_step;
    PS one_minus_twosigma = PS(1) - PS(2) * sigma;
    Matrix<PS> w(PS(0), x_divide+1, t_divide+1);
    for (int i = 0; i <= x_divide;++i)
    {
        w(i, 0) = initial_cond(left_end + PS(i) * (right_end - left_end) / PS(x_divide));
    }
    for (int j = 0; j <= t_divide;++j)
    {
        w(0, j) = left_bond_cond(PS(j) * time / PS(t_divide));
        w(x_divide, j) = right_bond_cond(PS(j) * time / PS(t_divide));
    }
    for (int j = 1; j <= t_divide; ++j)
    {
        for (int i = 1; i < x_divide;++i)
        {
            w(i, j) = sigma * w(i + 1, j - 1) + one_minus_twosigma * w(i, j - 1) + sigma * w(i - 1, j - 1);
        }
    }
    return w;
}

}
#endif