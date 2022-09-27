#ifndef CORE_H
#define CORE_H
#include<iostream>
#include<iomanip>
#include<cmath>
#include<functional>
#include<vector>
using std::cerr;
using std::cout;
using std::function;
using std::vector;

namespace cor{

//-----常量部分-----
const double Machine_Epsilon = pow(2, -52);
template<class DC> const DC Pi = 3.1415926535897932;

//-----声明部分-----

//杂例I
template <class DM>
double Norm(const vector<DM> &, double);
template<class DJ> DJ Abs(DJ);
template <class DJ> DJ Sign(DJ);
template <class DJ> vector<DJ> Range(DJ, DJ,DJ);
template <class DJ> vector<DJ> Range(DJ, DJ);
template <class DJ> vector<DJ> Rand(int);
template <class DJ> DJ Limit(std::function<DJ(DJ)>, DJ);
template <class DJ> DJ N(DJ,int);
inline int Sig(){return rand()%2?1:-1;}//random Sign(+/-)

template <class DD> DD Iteration(std::function<DD(DD)>, DD);
template <class DD> DD Sqrt(DD);
template <> double Sqrt(double);


//-----定义部分-----


//杂例I
template<class DM> double Norm(const vector<DM> &a, double p)
{
    if(p==INFINITY)
    {
        double max = 0;
        for (long long unsigned int i = 0; i < a.size();++i)
        {
            if(Abs(a[i])>max)
                max = Abs(a[i]);
        }
        return max;
    }
    double s = 0;
    for (long long unsigned int i = 0; i < a.size();++i)
    {
        s += pow(Abs(a[i]), p);
    }
    s = pow(s, 1 / p);
    return s;
}
template<class DJ> DJ Abs(DJ a)
{
    if(a>0)
        return a;
    return -a;
}
template<class DJ> DJ Sign(DJ a)
{
    if(a>0)
        return 1;
    if(a<0)
        return -1;
    return 0;
}
template <class DJ> vector<DJ> Range(DJ a, DJ b,DJ c)
{
    vector<DJ> r;
    DJ s = a;
    while(s<=b)
    {
        r.emplace_back(s);
        s = s + DJ(c);
    }
    return r;
}
template <class DJ> vector<DJ> Range(DJ a, DJ b)
{
    return Range(a, b, DJ(1));
}
template <class DJ> vector<DJ> Rand(int a)
{
    vector<DJ> r(a);
    for (int i = 0;i<a;++i)
        r[i]=std::rand()*Sig();
    return r;
}
template <class DJ> DJ Limit(std::function<DJ(DJ)> f, DJ x0)
{
    const DJ h = 1e-7;
    return 0.5 * (f(x0+h)+f(x0-h));
}
template <class DJ> DJ N(DJ a,int n)
{
    cout << std::setprecision(n) << a << '\n';
    return a;
}

template <class DD> DD Iteration(std::function<DD(DD)> phi, DD a)
{
    int times = 100;//迭代次数
    int k = 0;
    DD pre=a;
    DD x = phi(a);
    while (Abs(x-pre)>1e-14 && k < times)
    {
        pre = x;
        x = phi(x);
        ++k;
    } 
    if(k>=times)
    {
        cerr<<"警告：未在指定次数内达到预期精度。"<<'\n';
        cerr << "可能是初值的选择使迭代法不收敛，或者迭代函数的选择使收敛速度很慢。" << '\n';
    }
    return x;
}
template<class DD> DD Sqrt(DD c)
{
    std::function<DD(DD)> phi = [&](DD x) { return (x + c / x) / 2; };
    return Iteration(phi,c);
}
template<> double Sqrt<double>(double c)
{
    if(c<0)
    {
        cerr << "错误：根号下不能是负数。" << '\n';
        return NAN;
    }
    double pre=0;
    double x = 0.5*c;
    while (fabs(pre-x)>1e-14)
    {
        pre = x;
        x = 0.5 * (x + c / x);
    }
    return x;
}

}
#endif