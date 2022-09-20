#ifndef NUMBERTHEORY_H
#define NUMBERTHEORY_H
#include "FFT.h"
#include<iostream>
#include<vector>
#include<tuple>
#include<unordered_map>
using std::get;
using std::ostream;
using std::pair;
using std::tuple;
using std::unordered_map;
using std::vector;
using namespace fft;

namespace nbt{

//-----声明部分-----

//整数环
template <class ZZ> ZZ GCD(ZZ, ZZ); //最大公因数
template <class ZZ> ZZ LCM(ZZ, ZZ); //最小公倍数
template <class ZZ> tuple<ZZ, ZZ, ZZ> ExtendedGCD(ZZ a, ZZ b); //扩展欧几里得算法

//素数筛
template <class PP> PP EulerSieve(PP, vector<bool> &, vector<PP> &); //欧拉筛
template <class PP> vector<PP> PrimeTable(PP); //小于等于某个数的所有素数列表
template <class PP> PP PrimePi(PP); //小于等于某个数的素数个数

//素性测试
template <class PP> pair<PP, PP> Decomposition2Adic(PP); //分解 n=2^q*m
template <class PP> bool MillerRabin(PP); //Miller-Rabin素性测试
template <class PP> bool SimplePrimeQ(PP); //朴素素性测试
template <class PP> bool PrimeQ(PP, const string &); //素性测试，可调节方法
template <class PP> bool PrimeQ(PP); //素性测试
template <class PP> PP NextPrime(PP); //下一个素数

//素因子分解
template <class PP> PP MinimalPrimeFactor(PP); //最小素因子
template <class PP> PP PollandRho(PP); //Polland-Rho法找整数的一个因子
template <class PP> vector<PP> UnorderedFactorList(PP); //无序素因子列表
template <class PP> unordered_map<PP,PP> FactorInteger(PP); //整数的素因子分解

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

        template <class FG> friend FG Numerator(const Frac<FG> &);
        template <class FG> friend FG Denominator(const Frac<FG> &);
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
    if (b==0)
        return a;
    return GCD(b, a % b);
}
template <class ZZ> ZZ LCM(ZZ a, ZZ b)
{
    return a / GCD(a, b) * b;
}
template <class ZZ> tuple<ZZ, ZZ, ZZ> ExtendedGCD(ZZ a, ZZ b)
{
    if (b==0)
        return std::make_tuple(a,1,0);
    tuple<ZZ, ZZ, ZZ> t = ExtendedGCD(b, a % b);
    return std::make_tuple(get<0>(t), get<2>(t), get<1>(t) - (a / b) * get<2>(t));
}

//素数筛
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
template<class PP> PP PrimePi(PP n)
{
    vector<bool> ifprime(n + 1, true);
    vector<PP> primetable;
    return EulerSieve(n, ifprime, primetable);
}

//素性测试
template<class PP> pair<PP,PP> Decomposition2Adic(PP m)
{
    PP q=0;
    while(!(m&1))
    {
        ++q;
        m >>= 1;
    }
    return std::make_pair(q, m);
}
template<class PP> bool Witness(PP a,PP m,PP q,PP n)
{
    PP x= ModPower(a, m, n);
    for (PP i = 0; i < q; ++i)
    {
        PP x_pre = x;
        x = x * x % n;
        if (x==1 && x_pre!=1 && x_pre!=n-1)
            return true;
    }
    if (x!=1)
        return true;
    return false;
}
template<class PP> bool MillerRabin(PP n)
{ 
    if(n==2 || n==3 || n==5)
        return true;
    if(n%2==0 || n%3==0 || n%5==0 || n==1)
        return false;
    pair<PP, PP> p = Decomposition2Adic(n - 1);
    PP q = get<0>(p);
    PP m = get<1>(p);
    srand(time(0));
    for (PP i = 1; i < n; i<<=1) 
    {
        PP a =1+(rand()%(n-1));
        if(Witness(a,m,q,n))
            return false;   // definitely composite  
    }
    return true;    // almost surely prime
}
template<class PP> bool SimplePrimeQ(PP n)
{
    if(n==2 || n==3 || n==5)
        return true;
    if(n%2==0 || n%3==0 || n%5==0 || n==1)
        return false;
    const std::array<PP,8> step({4, 2, 4, 2, 4, 6, 2, 6}); //30n+k法
    PP c = 7;
    PP sqrtn = sqrt(n);
    while(c <= sqrtn)
        for (auto i:step)
        {
            if (n%c==0)
                return false;
            c += i;
        }
    return true; 
}
template<class PP> bool PrimeQ(PP n, const string &str)
{
    if (str=="Miller-Rabin" || str=="MillerRabin")
        return MillerRabin(n);
    if (str=="simple" || str=="trivial")
        return SimplePrimeQ(n);
    std::cerr << "错误：没有找到该方法！" << std::endl;
    return false;
}
template<class PP> bool PrimeQ(PP n)
{
    if (n<2000000000) //实测小于这个数据时simple更快
        return PrimeQ(n, "simple");
    return PrimeQ(n, "Miller-Rabin");
}
template<class PP> PP NextPrime(PP n)
{
    while(!MillerTest(++n));
    return n;
}

//素因子分解
template<class PP> PP MinimalPrimeFactor(PP n)
{
    if (n == 1) return 1;
    if (n%2 == 0) return 2;
    if (n%3 == 0) return 3;
    if (n%5 == 0) return 5;
    const std::array<PP,8> step({4, 2, 4, 2, 4, 6, 2, 6}); //30n+k法
    PP c = 7;
    PP sqrtn = sqrt(n);
    while(c <= sqrtn)
        for (auto i:step)
        {
            if (n%c==0)
                return c;
            c += i;
        }
    return n;
}
template<class PP> PP PollandRho(PP n)
{
    if (n == 1) return 1;
    if (n%2 == 0) return 2;
    if (n%3 == 0) return 3;
    if (n%5 == 0) return 5;
    if (PrimeQ(n)) return n;
    PP i = 1;
    srand(time(0));
    PP x = rand()%n;
    PP y = x; 
    PP k = 2;
    while(true)
    {
        ++i;
        x = (x * x - 1) % n;
        PP m = y > x ? y - x : x - y;
        PP d = GCD(m, n);
        if (d!=1 && d!=n)
            return d;
        if (i==k)
        {
            y = x;
            k <<= 1;
        }
    }
    return 2;
}
template<class PP> PP FindAFactor(PP n)
{
    if (n<2000000000) //实测小于这个数据时minimalprimefactor更快
        return MinimalPrimeFactor(n);
    return PollandRho(n);
}
template<class PP> void MakeFactorList(PP n, vector<PP> &v)
{
    if(n==0 || n==1)
        return;
    if(PrimeQ(n))
    {
        v.emplace_back(n);
        return;
    }
    PP m = FindAFactor(n);
    MakeFactorList(m, v);
    MakeFactorList(n / m, v);
}
template<class PP> vector<PP> UnorderedFactorList(PP n)
{
    vector<PP> v;
    MakeFactorList(n, v);
    return v;
}
template<class PP> unordered_map<PP,PP> FactorInteger(PP n)
{
    vector<PP> v = UnorderedFactorList(n);
    unordered_map<PP, PP> map;
    for(PP i:v)
    {
        if(map.find(i) != map.end()) // found
            ++map[i];
        else // not found
            map.insert({i, 1});
    }
    return map;
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

template <class FG> FG Numerator(const Frac<FG> &a)
{
    return a.numerator;
}
template <class FG> FG Denominator(const Frac<FG> &a)
{
    return a.denominator;
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