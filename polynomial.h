#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include "FFT.h"
using std::cerr;
using std::cout;
using std::ostream;
using namespace fft;

namespace ply{

//-----声明部分-----

//多项式类
template<class DC>
class Polynomial{
public:
    Polynomial();
    Polynomial(const vector<DC>&);
    Polynomial(const vector<DC>&, int);
    Polynomial(DC);

    Polynomial<DC> &operator=(const Polynomial<DC>&);
    DC operator()(DC);
    template <class DD> friend Polynomial<DD> operator+(const Polynomial<DD> &,const Polynomial<DD> &);
    template <class DD> friend Polynomial<DD> operator-(const Polynomial<DD> &,const Polynomial<DD> &);
    template <class DD> friend Polynomial<DD> operator+(DD, const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> operator*(DD, const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> operator*(const Polynomial<DD>&,const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> operator/(const Polynomial<DD>&,const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> operator^(const Polynomial<DD>&,int);
    template<class DD> friend ostream &operator << (ostream &,const Polynomial<DD> &);
    template <class DD> friend Polynomial<DD> operator>>(const Polynomial<DD>&, int);
    Polynomial<DC> &operator+=(const Polynomial<DC>&);
    Polynomial<DC> &operator-=(const Polynomial<DC>&);

    template <class DD> friend DD Get(const Polynomial<DD> &, DD);
    template <class DD> friend vector<DD> GetCoef(const Polynomial<DD> &);
    template <class DD> friend void Disp(const Polynomial<DD>&);
    template <class DD> friend int Deg(const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> SubPoly(const Polynomial<DD>&, int, int);
    template <class DD> friend Polynomial<DD> D(const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> Integrate(const Polynomial<DD> &);
    template <class DD> friend DD Integrate(const Polynomial<DD> &,DD,DD);
    template <class DD> friend Polynomial<DD> Schmidt(int,DD,DD);
    template <class DD> friend Polynomial<DD> Legendre(int n);

private:
    int poly_deg;
    vector<DC> poly_coeff;
};
Polynomial<double> operator*(const Polynomial<double>&,const Polynomial<double>&);
Polynomial<float> operator*(const Polynomial<float>&,const Polynomial<float>&);
Polynomial<complex<double>> operator*(const Polynomial<complex<double>> &, const Polynomial<complex<double>> &);
Polynomial<complex<float>> operator*(const Polynomial<complex<float>> &, const Polynomial<complex<float>> &);
/*Polynomial<unsigned> operator*(const Polynomial<unsigned> &, const Polynomial<unsigned> &);
Polynomial<unsigned long long> operator*(const Polynomial<unsigned long long> &, const Polynomial<unsigned long long> &);*/ /*test*/

//-----定义部分-----

//多项式类
const Polynomial<double> X({0, 1}, 1);
template<class DC> Polynomial<DC>::Polynomial()
{
    poly_deg = 0;
    poly_coeff.push_back(0);
}
template<class DC> Polynomial<DC>::Polynomial(const vector<DC> &coeff)
{
    poly_deg = coeff.size() - 1;
    poly_coeff = coeff;
}
template<class DC> Polynomial<DC>::Polynomial(const vector<DC> &s,int deg)
{
    poly_deg = deg;
    poly_coeff = s;
}
template<class DC> Polynomial<DC>::Polynomial(DC r)
{
    poly_deg = 0;
    poly_coeff.push_back(r);
}
template<class DC> Polynomial<DC> &Polynomial<DC>::operator=(const Polynomial<DC> &f)
{
    poly_deg = f.poly_deg;
    poly_coeff=f.poly_coeff;
    return *this;
}
template<class DC> DC Polynomial<DC>::operator()(DC t)
{
    DC s=poly_coeff[poly_deg];
    for (int i = poly_deg-1; i >= 0;--i)
        s = s * t + poly_coeff[i];
    return s;
}
template<class DC> Polynomial<DC> operator+(const Polynomial<DC> &a,const Polynomial<DC> &b)
{
    bool here = true;
    int degmax = a.poly_deg;
    if(b.poly_deg>a.poly_deg)
    {
        degmax=b.poly_deg;
        here = false;
    }
    vector<DC> r(degmax + 1);
    Polynomial<DC> c;
    c.poly_deg = degmax;
    c.poly_coeff = r;
    if(here)
    {
        for (int i = 0; i <=b.poly_deg;++i)
        {
            c.poly_coeff[i] = a.poly_coeff[i] + b.poly_coeff[i];
        }
        for (int i = b.poly_deg + 1;i<=c.poly_deg;++i)
        {
            c.poly_coeff[i] = a.poly_coeff[i];
        }
        return c;
    }
    for (int i = 0; i <=a.poly_deg;++i)
    {
        c.poly_coeff[i] = a.poly_coeff[i] + b.poly_coeff[i];
    }
    for (int i = a.poly_deg + 1;i<=c.poly_deg;++i)
    {
        c.poly_coeff[i] = b.poly_coeff[i];
    }
    return c;
}
template<class DC> Polynomial<DC> operator-(const Polynomial<DC> &a,const Polynomial<DC> &b)
{
    bool here = true;
    int degmax = a.poly_deg;
    if(b.poly_deg>a.poly_deg)
    {
        degmax=b.poly_deg;
        here = false;
    }
    vector<DC> r(degmax + 1);
    Polynomial<DC> c;
    c.poly_deg = degmax;
    c.poly_coeff = r;
    if(here)
    {
        for (int i = 0; i <=b.poly_deg;++i)
        {
            c.poly_coeff[i] = a.poly_coeff[i] - b.poly_coeff[i];
        }
        for (int i = b.poly_deg + 1;i<=c.poly_deg;++i)
        {
            c.poly_coeff[i] = a.poly_coeff[i];
        }
        return c;
    }
    for (int i = 0; i <=a.poly_deg;++i)
    {
        c.poly_coeff[i] = a.poly_coeff[i] - b.poly_coeff[i];
    }
    for (int i = a.poly_deg + 1;i<=c.poly_deg;++i)
    {
        c.poly_coeff[i] = DC(0)-b.poly_coeff[i];
    }
    return c;
}
template <class DD> Polynomial<DD> operator+(DD a,const Polynomial<DD> &f)
{
    Polynomial<DD> p = f;
    p.poly_coeff[0] += a;
    return p;
}
template <class DD> Polynomial<DD> operator*(DD a,const Polynomial<DD> &f)
{
    vector<DD> r;
    for (int i = 0; i <= f.poly_deg;++i)
    {
        r.push_back(a * f.poly_coeff[i]);
    }
    Polynomial<DD> p;
    p.poly_deg = f.poly_deg;
    p.poly_coeff = r;
    return p;
}
template <class DD> Polynomial<DD> operator*(const Polynomial<DD> &a, const Polynomial<DD> &b)
{
    int n = a.poly_deg > b.poly_deg ? a.poly_deg : b.poly_deg;
    if (a.poly_deg == 0)
        return a.poly_coeff[0] * b;
    if(b.poly_deg==0)
        return b.poly_coeff[0] * a;
    if(n==1)
    {
        vector<DD> r(3);
        r[0]=a.poly_coeff[0]*b.poly_coeff[0];
        r[2]=a.poly_coeff[1]*b.poly_coeff[1];
        r[1] = (a.poly_coeff[0] + a.poly_coeff[1]) * (b.poly_coeff[0] + b.poly_coeff[1]) - r[0] - r[2];
        Polynomial<DD> f;
        f.poly_deg=2;
        f.poly_coeff = r;
        return f;
    }
    Polynomial<DD> p=a;
    Polynomial<DD> q=b;
    if(a.poly_deg>b.poly_deg)
    {
        q.poly_coeff.insert(q.poly_coeff.end(), a.poly_deg - b.poly_deg, 0);
        q.poly_deg = n;
    }
    if(b.poly_deg>a.poly_deg)
    {
        p.poly_coeff.insert(p.poly_coeff.end(), b.poly_deg - a.poly_deg, 0);
        p.poly_deg = n;
    }
    Polynomial<DD> p0 = SubPoly(p, 0, n / 2);
    Polynomial<DD> p1 =SubPoly(p, n/2+1, n);
    Polynomial<DD> q0=SubPoly(q,0,n/2);
    Polynomial<DD> q1 = SubPoly(q, n/2+1, n);
    Polynomial<DD> r0 = p0 * q0;
    Polynomial<DD> r1=p1*q1;
    Polynomial<DD> r2=(p0+p1)*(q0+q1);
    return r0 + ((r2 - r0 - r1) >> (n/2+1)) + (r1 >> ((n/2+1)<<1));
}
template <class DD> Polynomial<DD> operator/(const Polynomial<DD>& f,const Polynomial<DD>& g)
{
    cerr<<"尚待开发"<<'\n';
    return f;
}
template<class DC> Polynomial<DC> operator^(const Polynomial<DC> &a,int n)
{
    if(n<0)
    {
        cerr << "错误：多项式不能有负的幂次。" << '\n';
        return a;
    }
    Polynomial<DC> m=a;
    Polynomial<DC> b({1},0);
    while(n>=1)
    {
        if(n&1)
        {
            b = m * b;
        }
        n=n>>1;
        m = m * m;
    }
    return b;
}
template<class DC> ostream & operator<<(ostream & os, const Polynomial<DC> &f)
{
    os << f.poly_coeff[0];
    for (int i =1;i<=f.poly_deg;++i)
    {
        if(f.poly_coeff[i]==0)
            continue;
        if(f.poly_coeff[i]>0)
        {
            os << '+';
        }
        os << f.poly_coeff[i] << "x" ;
        if(i!=1)
        {
            os << '^'<<i;
        }
    }
    return os;
}
template<class DC> Polynomial<DC> &Polynomial<DC>::operator+=(const Polynomial<DC> &b)
{
    *this = (*this) + b;
    return *this;
}
template<class DC> Polynomial<DC> &Polynomial<DC>::operator-=(const Polynomial<DC> &b)
{
    *this = (*this)- b;
    return *this;
}
template <class DD> Polynomial<DD> operator>>(const Polynomial<DD> &f, int a)
{
    vector<DD> r(f.poly_deg+1+a);
    for (int i = 0; i <= f.poly_deg;++i)
    {
        r[i + a] = f.poly_coeff[i];
    }
    Polynomial<DD> p;
    p.poly_deg = f.poly_deg + a;
    p.poly_coeff = r;
    return p;
}
template <class DD> DD Get(const Polynomial<DD> &f, DD t)
{
    DD s=f.poly_coeff[f.poly_deg];
    for (int i = f.poly_deg-1; i >= 0;--i)
        s = s * t + f.poly_coeff[i];
    return s;
}
template <class DD> vector<DD> GetCoef(const Polynomial<DD> &f)
{
    return f.poly_coeff;
}
template<class DC> void Disp(const Polynomial<DC> &f)
{
    cout << f<<'\n'<<'\n';
}
template<class DC> int Deg(const Polynomial<DC> &f)
{
    return f.poly_deg;
}
template <class DC> Polynomial<DC> SubPoly(const Polynomial<DC> &f, int a, int b)
{
    if(a<0)
        a = 0;
    if(b>f.poly_deg)
        b = f.poly_deg;
    vector<DC> r;
    for (int i = a; i <= b;++i)
    {
        r.push_back(f.poly_coeff[i]);
    }
    Polynomial<DC> p;
    p.poly_deg=b-a;
    p.poly_coeff = r;
    return p;
}
template <class DC> Polynomial<DC> D(const Polynomial<DC> &f)
{
    if(f.poly_deg==0)
    {
        Polynomial<DC> d({0},0);
        return d;
    }
    int n=f.poly_deg;
    vector<DC> r;
    for (int i = 1; i <= n;++i)
    {
        r.push_back(i * f.poly_coeff[i]);
    }
    Polynomial<DC> d;
    d.poly_deg = n - 1;
    d.poly_coeff = r;
    return d;
}
template <class DD> Polynomial<DD> Integrate(const Polynomial<DD> &f)
{
    if(f.poly_deg==0 && f.poly_coeff[0]==0)
    {
        return Polynomial<DD>({0}, 1);
    }
    int deg = f.poly_deg;
    vector<DD> r(deg+2);
    for (int i = 1; i <= deg + 1;++i)
    {
        r[i] = f.poly_coeff[i - 1] / i;
    }
    return Polynomial<DD>(r, deg + 1);
}
template <class DD> DD Integrate(const Polynomial<DD> &f,DD a,DD b)
{
    Polynomial<DD> g=Integrate(f);
    return g(b) - g(a);
}
template <class DD> Polynomial<DD> Schmidt(int n,DD a,DD b)
{
    if(n<0)
    {
        return Polynomial<DD>({0}, 0);
    }
    vector<Polynomial<DD>> phi(n+1);
    vector<DD> innerproduct(n);
    phi[0] = Polynomial<DD>({1}, 0);
    for (int k = 1; k <= n;++k)
    {
        vector<DD> r(k+1);
        r[k] = DD(1);
        Polynomial<DD> f(r, k);
        phi[k] = f;
        innerproduct[k - 1] = Integrate(phi[k - 1] * phi[k - 1], a, b);
        for (int j = 0; j < k;++j)
        {
            phi[k] = phi[k] - (Integrate(phi[j] * f, a, b) / innerproduct[j])*phi[j];
        }
    }
    return phi[n];
}
template<class DD> Polynomial<DD> Legendre(int n)
{
    if(n<0)
        return Polynomial<DD>({0}, 0);
    vector<Polynomial<DD>> p({
        Polynomial<DD>({1}, 0),
        Polynomial<DD>({0, 1}, 1),
        Polynomial<DD>({-0.5, 0, 1.5}, 2),
        Polynomial<DD>({0, -1.5, 0, 2.5}, 3),
        Polynomial<DD>({0.375, 0, -3.75, 0, 4.375}, 4)
    });
    if(n<=4)
        return p[n];
    Polynomial<DD> pre1 = p[4];
    Polynomial<DD> pre2 = p[3];
    Polynomial<DD> now;
    for (int k = 5; k <= n;++k)
    {
        now = (DD(2 *k- 1) /k) *(pre1>>1) - (DD(k - 1) /k) * pre2;
        pre2 = pre1;
        pre1 = now;
    }
    return now;
}

Polynomial<double> operator*(const Polynomial<double>& f,const Polynomial<double>& g)
{
    return Polynomial<double>(RealConvolution(GetCoef(f), GetCoef(g)));
}
Polynomial<float> operator*(const Polynomial<float>& f,const Polynomial<float>& g)
{
    return Polynomial<float> (RealConvolution(GetCoef(f), GetCoef(g)));
}
Polynomial<complex<double>> operator*(const Polynomial<complex<double>>& f,const Polynomial<complex<double>>& g)
{
    return Polynomial<complex<double>>(Convolution(GetCoef(f), GetCoef(g)));
}
Polynomial<complex<float>> operator*(const Polynomial<complex<float>> &f, const Polynomial<complex<float>> &g)
{
    return Polynomial<complex<float>> (Convolution(GetCoef(f), GetCoef(g)));
}
/*Polynomial<unsigned> operator*(const Polynomial<unsigned> &f, const Polynomial<unsigned> &g)
{
    return Polynomial<unsigned>(IntConvolution(GetCoef(f), GetCoef(g)));
}
Polynomial<unsigned long long> operator*(const Polynomial<unsigned long long> &f, const Polynomial<unsigned long long> &g)
{
    return Polynomial<unsigned long long>(IntConvolution(GetCoef(f), GetCoef(g)));
}*/ /*test*/

}
#endif