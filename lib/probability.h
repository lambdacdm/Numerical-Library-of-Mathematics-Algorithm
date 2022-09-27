#ifndef PROBABILITY_H
#define PROBABILITY_H
#include "numerical.h"
using namespace nmr;

namespace prb{

//-----声明部分-----

//随机变量
class DiscreteRandomVariable
{
    public:
        DiscreteRandomVariable(const Matrix<double>&);
        friend Matrix<double> Distribution(const DiscreteRandomVariable &);
    private:
        Matrix<double> _distribution;
};
Matrix<double> Distribution(const DiscreteRandomVariable &);
double Expectation(Matrix<double> );
double Expectation(DiscreteRandomVariable X);

//分布
class NormalDistribution
{
    public:
        NormalDistribution(double, double);
        friend double Expectation(const NormalDistribution&);
        friend double StandardDeviation(const NormalDistribution&);

    private:
        double _mu;
        double _sigma;
};
double Expectation(const NormalDistribution &);
double StandardDeviation(const NormalDistribution &);
double PDF(const NormalDistribution &, double );
double CDF(const NormalDistribution &, double );

//-----定义部分-----

//随机变量
DiscreteRandomVariable::DiscreteRandomVariable(const Matrix<double> &M)
{
    _distribution = M;
}
Matrix<double> Distribution(const DiscreteRandomVariable &X)
{
    return X._distribution;
}
double Expectation(Matrix<double> A)
{
    if(RowSize(A)!=2)
    {
        cerr << "错误：输入的矩阵不是分布列" << '\n';
        return 0;
    }
    double s = 0;
    for (int j= 0; j< ColumnSize(A);++j)
    {
        s = s + A[0][j] * A[1][j];
    }
    return s;
}
double Expectation(DiscreteRandomVariable X)
{
    return Expectation(Distribution(X));
}

//分布
NormalDistribution::NormalDistribution(double mu,double sigma)
{
    _mu = mu;
    _sigma = sigma;
}
double Expectation(const NormalDistribution &N)
{
    return N._mu;
}
double StandardDeviation(const NormalDistribution &N)
{
    return N._sigma;
}
double PDF(const NormalDistribution &N,double x)
{
    const double mu = Expectation(N);
    const double sigma = StandardDeviation(N);
    return 1 / (sqrt(2 * Pi<double>)*sigma)*exp(-pow(x-mu,2)/(2*sigma*sigma));
}
double CDF(const NormalDistribution &N,double x)
{
    double leftend = Expectation(N) - 6 * StandardDeviation(N);
    if(x<=leftend)
        return 0;
    auto p = [&](double t) { return PDF(N, t); };
    return Integrate<double>(p, leftend, x);
}

}

#endif