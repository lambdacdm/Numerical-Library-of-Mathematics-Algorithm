#ifndef NUMBERTHEORY_H
#define NUMBERTHEORY_H
#include<iostream>
#include<vector>
using std::ostream;
using std::vector;

namespace nbt{

//-----声明部分-----

//整数环
template <class ZZ> ZZ GCD(ZZ a, ZZ b); //最大公因数
template <class ZZ> ZZ LCM(ZZ a, ZZ b); //最小公倍数

//素数
template <class PP> PP EulerSieve(PP, vector<bool> &, vector<PP> &); //欧拉筛
template <class PP> vector<PP> PrimeTable(PP); //小于等于某个数的所有素数列表
template <class PP> PP PrimeCount(PP); //小于等于某个数的素数个数

//分式域
template <class FF> class Frac
{
    public:
        Frac();
        Frac(FF);
        Frac(FF, FF);

        template <class FG> friend bool operator==(const Frac<FG> &, const Frac<FG> &);
        template <class FG> friend Frac<FG> operator+(const Frac<FG> &, const Frac<FG> &);
        template <class FG> friend Frac<FG> operator-(const Frac<FG> &, const Frac<FG> &);
        template <class FG> friend Frac<FG> operator*(const Frac<FG> &, const Frac<FG> &);
        template <class FG> friend Frac<FG> operator/(const Frac<FG> &, const Frac<FG> &);
        template <class FG> friend ostream &operator<<(ostream &, const Frac<FG> &);

        template <class FG> friend Frac<FG> Reciprocal(const Frac<FG> &);
        template <class FG> friend Frac<FG> Simplify(const Frac<FG> &);
        template <class FG> friend Frac<FG> FullSimplify(const Frac<FG> &);
        template <class FG, class FH> friend FH Evaluate(const Frac<FG> &);

    private:
        FF numerator;
        FF denominator;
};

//-----定义部分-----

//整数环
template <class ZZ> ZZ GCD(ZZ a, ZZ b)
{
    while (b!=0)
    {
        ZZ t = a % b;
        a = b;
        b = t;
    }
    return a;
}
template <class ZZ> ZZ LCM(ZZ a, ZZ b)
{
    return a / GCD(a, b) * b;
}

//素数
template<class PP> PP EulerSieve(PP n,vector<bool> &ifprime,vector<PP> &primetable)
{
    PP count = 0;
    for (PP i = 2; i <= n;++i)
    {
        if(ifprime[i])
        {
            primetable.emplace_back(i);
            ++count;
        }
        for (auto j:primetable)
        {
            if(j*i>n)
                break;
            ifprime[j*i]=false;
            if(i%j==0)
                break;
        }
    }
    return count;
}
template<class PP> vector<PP> PrimeTable(PP n)
{
    vector<bool> ifprime(n + 1, true);
    vector<PP> primetable;
    EulerSieve(n, ifprime, primetable);
    return primetable;
}
template<class PP> PP PrimeCount(PP n)
{
    vector<bool> ifprime(n + 1, true);
    vector<PP> primetable;
    return EulerSieve(n, ifprime, primetable);
}

//分式域
template <class FF> Frac<FF>::Frac()
{
    numerator = 0;
    denominator = 1;
}
template <class FF> Frac<FF>::Frac(FF p)
{
    numerator = p;
    denominator = 1;
}
template <class FF> Frac<FF>::Frac(FF p, FF q)
{
    numerator = p;
    denominator = q;
}

template <class FG> const Frac<FG> FractionalZero = Frac<FG>(0);
template <class FG> const Frac<FG> FractionalOne = Frac<FG>(1);
template <class FG> const Frac<FG> FractionalUndefined = Frac<FG>(0,0);
template <class FG> const Frac<FG> FractionalInf = Frac<FG>(1,0);

template <class FG> bool operator==(const Frac<FG> &a, const Frac<FG> &b)
{
    if(a.numerator==0 && a.denominator==0)
        return (b.numerator == 0 && b.denominator == 0);
    if(b.numerator==0 && b.denominator==0)
        return (a.numerator == 0 && a.denominator == 0);    
    return (a.numerator * b.denominator == a.denominator * b.numerator);
}
template <class FG> Frac<FG> operator+(const Frac<FG> &a, const Frac<FG> &b)
{
    return Frac<FG>(a.numerator * b.denominator + a.denominator * b.numerator, 
        a.denominator * b.denominator);
}
template <class FG> Frac<FG> operator-(const Frac<FG> &a, const Frac<FG> &b)
{
    return Frac<FG>(a.numerator * b.denominator - a.denominator * b.numerator, 
        a.denominator * b.denominator);
}
template <class FG> Frac<FG> operator*(const Frac<FG> &a, const Frac<FG> &b)
{
    return Frac<FG>(a.numerator * b.numerator, a.denominator * b.denominator);
}
template <class FG> Frac<FG> operator/(const Frac<FG> &a, const Frac<FG> &b)
{
    return Frac<FG>(a.numerator * b.denominator, a.denominator * b.numerator);
}
template <class FG> ostream & operator<<(ostream & os, const Frac<FG> &a)
{
    if(a.numerator==0 && a.denominator==0)
    {
        os << "UNDEFINED";
        return os;
    }
    if(a.denominator==0)
    {
        os << "INF";
        return os;
    }
    if(a.numerator==0)
    {
        os << "0";
        return os;
    }
    if(a.denominator==1)
    {
        os << a.numerator;
        return os;
    }
    os << a.numerator << '/' << a.denominator;
    return os;
}
template <class FG> Frac<FG> Reciprocal(const Frac<FG> &a)
{
    return Frac<FG>(a.denominator, a.numerator);
}
template <class FG> Frac<FG> Simplify(const Frac<FG> &a)
{
    if(a.numerator==0 && a.denominator==0)
        return a;
    FG g = GCD(a.numerator, a.denominator);
    return Frac<FG>(a.numerator / g, a.denominator / g);
}
template <class FG> Frac<FG> FullSimplify(const Frac<FG> &a)
{
    if(a.numerator==0 && a.denominator==0)
        return a;
    FG g = GCD(a.numerator, a.denominator);
    FG n = a.numerator / g;
    FG d = a.denominator / g;
    if(d<0)
    {
        n = -n;
        d = -d;
    }
    return Frac<FG>(n,d);
}
template <class FG, class FH> FH Evaluate(const Frac<FG> &a)
{
    return FH(a.numerator) / FH(a.denominator);
}

}
#endif