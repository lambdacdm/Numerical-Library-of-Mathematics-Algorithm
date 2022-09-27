#ifndef ODE_H
#define ODE_H
#include "numerical.h"
using namespace nmr;

namespace ode{

//-----声明部分-----

//常微分方程
template<class DS,class DO> DO DSolve(std::function<DO(DS,DO)>,const pair<DS,DO>&,DS,string);
template <class DS, class DO>DO DSolve(std::function<DO(DS, DO)>, const pair<DS, DO> &, DS, int);
template <class DS, class DO> DO DSolve(std::function<DO(DS, DO)>, const pair<DS, DO> &, DS);
template <class DS, class DO> 
std::function<DO(DS)> DSolve(std::function<DO(DS, DO)>, const pair<DS, DO> &, string);
template <class DS, class DO>
std::function<DO(DS)> DSolve(std::function<DO(DS, DO)>, const pair<DS, DO> &);

//-----定义部分-----

//常微分方程
template<class DS,class DO,class DP> DO Iteration(const DP &phi,const pair<DS,DO> &initial,DS x,int n)
{
    if(n==0)
        return get<1>(initial);
    DS h= (x - get<0>(initial))/n;//步长
    DS xk = get<0>(initial);
    DO yk= get<1>(initial);
    for (int i = 0; i < n;++i)
    {
        phi(xk, yk, h);
        xk = xk + h;
    }
    return yk;
}
template<class DS,class DO> DO RKSolve(std::function<DO(DS,DO)> f,const pair<DS,DO>& initial,DS x,int n)
{
    auto phi = [=](DS &xk, DO &yk, DS &h) {
            DS halfh = h / 2.0;
            DS xk_plus_halfh = xk + halfh;
            DO k1=f(xk,yk);
            DO k2 =f(xk_plus_halfh,yk+halfh*k1);
            DO k3 =f(xk_plus_halfh,yk+halfh*k2);
            DO k4=f(xk+h,yk+h*k3);
            yk = yk + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        };
    return Iteration(phi, initial, x, n);
}
template<class DS,class DO,class DP,class DF> 
DO LinearMultiStep(const DP &phi,const DF &f,const pair<DS,DO> &initial,DS x,int n,int startstep)
{
    if(n==0)
        return get<1>(initial);
    DS h = (x - get<0>(initial))/n;
    std::deque<DS> xk={get<0>(initial)};
    std::deque<DO> yk={get<1>(initial)};
    std::deque<DO> fk={f(xk[0],yk[0])};
    for (int i =1; i<=std::min(startstep-1,n);++i)
    {
        yk.push_back(RKSolve(f,{xk.back(),yk.back()},xk.back()+h,1));
        xk.push_back(xk.back()+h);
        fk.push_back(f(xk.back(), yk.back()));
    }
    DO y0adpt = yk.back();
    for (int i = startstep; i <=n;++i)
    {
        phi(xk, yk, fk,h,y0adpt);
        xk.push_back(xk.back() + h);
        fk.push_back(f(xk.back(), yk.back()));
        xk.pop_front();
        yk.pop_front();
        fk.pop_front();
    }
    return yk.back();
}
template<class DS,class DO> DO DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO>& initial,DS x,int n,string str)
{
    if(str=="Euler"|| str=="forward Euler" || str=="explicit Euler")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) 
        {   yk = yk + h * f(xk, yk);
        };
        return Iteration(phi, initial, x, n);
    }
    if(str=="backward Euler" || str=="implicit Euler")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) {
            DS x_plus_h = xk + h;
            DO y=yk+h*f(xk,yk);
            for (unsigned i = 0; i <3;++i)
                y = yk + h * f(x_plus_h, y);
            yk = yk + h * f(x_plus_h, y);
        };
        return Iteration(phi, initial, x, n);
    }
    if(str=="trapz" || str=="trapezoidal")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) {
            DO fk=f(xk,yk);
            DS halfh = h / 2.0;
            DS x_plus_h = xk + h;
            DO y0 = yk + h * fk;
            DO y1=yk+halfh*(fk+f(x_plus_h,y0));
            yk = yk + halfh * (fk + f(x_plus_h, y1));
        };
        return Iteration(phi, initial, x, n);
    }
    if(str=="Heun" || str=="improved Euler")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) {
            DO k1=f(xk,yk);
            DO k2=f(xk+h,yk+h*k1);
            yk=yk+h/2.0*(k1+k2);
        };
        return Iteration(phi, initial, x, n);
    }
    if(str=="midpoint")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) {
            DS halfh = h / 2.0;
            DO k1=f(xk,yk);
            DO k2=f(xk+halfh,yk+halfh*k1);
            yk=yk+h*k2;
        };
        return Iteration(phi, initial, x, n);
    }
    if(str=="RK4" || str=="Runge-Kutta" || str=="RK")
    {
        return RKSolve(f, initial, x, n);
    }
    if(str=="Adams-Bashforth")
    {
        auto phi = [=](deque<DS> &xk, deque<DO> &yk, deque<DO> &fk, DS &h, DO &y0adpt) {
            yk.push_back(yk[1] + 3.0 / 2.0 * h * fk[1] - 1.0 / 2.0 * h * fk[0]);
        };
        return LinearMultiStep(phi, f, initial, x, n, 2);
    }
    if(str=="Adams")
    {
        auto phi = [=](deque<DS> &xk, deque<DO> &yk, deque<DO> &fk, DS &h,DO &y0adpt) {
            yk.push_back(yk[3]+h/24.0*(55.0*fk[3]-59.0*fk[2]+37.0*fk[1]-9.0*fk[0]));
        };
        return LinearMultiStep(phi, f, initial, x, n, 4);
    }
    if(str=="Milne")
    {
        auto phi = [=](deque<DS> &xk, deque<DO> &yk, deque<DO> &fk, DS &h,DO &y0adpt) {
            yk.push_back(yk[0]+4.0*h/3.0*(2.0*fk[3]-fk[2]+2.0*fk[1]));
        };
        return LinearMultiStep(phi, f, initial, x, n, 4);
    }
    if(str=="adaptive Adams")
    {
        auto phi = [=](deque<DS> &xk, deque<DO> &yk, deque<DO> &fk, DS &h,DO &y0adpt) {
            DO y0=yk[3]+h/24.0*(55.0*fk[3]-59.0*fk[2]+37.0*fk[1]-9.0*fk[0]);
            DO y1 = y0 + 251.0 * (yk[3] - y0adpt)/270.0;
            DO y2 = yk[3] + h / 24.0 * (9.0 * f(xk[3] + h, y1) + 19.0 * fk[3] - 5.0* fk[2] + fk[1]);
            yk.push_back(y2-19.0*(y2-y0)/270.0);
            y0adpt = y0;
        };
        return LinearMultiStep(phi, f, initial, x, n, 4);
    }
    if(str=="adaptive Milne-Hamming")
    {
        auto phi = [=](deque<DS> &xk, deque<DO> &yk, deque<DO> &fk, DS &h,DO &y0adpt) {
            DO y0=yk[0]+4.0*h/3.0*(2.0*fk[3]-fk[2]+2.0*fk[1]);
            DO y1=y0 + 112.0*(yk[3]-y0adpt)/121.0;
            DO y2 = (9.0 * yk[3] - yk[1]) / 8.0 + 3.0 * h / 8.0 * (f(xk[3] + h, y1)+ 2.0 * fk[3] - fk[2]);
            yk.push_back(y2 - 9.0 * (y2 - y0) / 121.0);
            y0adpt = y0;
        };
        return LinearMultiStep(phi, f, initial, x, n, 4);
    }
    cerr<<"错误：没有定义该方法"<<'\n';
    return x;
}
template<class DS,class DO> DO DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO>& initial,DS x,string str)
{
    return DSolve(f, initial, x, 100, str);
}
template<class DS,class DO> DO DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO>& initial,DS x,int n)
{
    return DSolve(f, initial, x, n,"adaptive Milne-Hamming");
}
template<class DS,class DO> DO DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO>& initial,DS x)
{
    return DSolve(f, initial, x, 100);
}
template<class DS,class DO> std::function<DO(DS)> DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO> &initial,string str)
{
    return [=](DS x) { return DSolve(f, initial, x, str); };
}
template<class DS,class DO> std::function<DO(DS)> DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO> &initial)
{
    return [=](DS x) { return DSolve(f, initial, x); };
}

}
#endif