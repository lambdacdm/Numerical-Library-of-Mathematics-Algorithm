#ifndef FFT_H
#define FFT_H
#include "core.h"
#include<vector>
#include<cstdlib>
#include<complex>
#include<algorithm>
#include<iterator>
#include<cstdint>
using std::complex;
using std::string;
using namespace cor;

namespace fft{

//-----定义部分-----

//快速傅里叶变换
inline uint32_t Rev(uint32_t k,uint8_t m)
{
    uint32_t s = 0;
    while(m>0)
    {
        --m;
        s += ((k%2) << m);
        k = (k >> 1);
    }
    return s;
}
template<class DC> vector<DC> bit_reverse_copy(const vector<DC> &a)
{
    uint32_t n = a.size();
    uint8_t m = log2(n);
    vector<DC> A(n);
    for(uint32_t k=0;k<n;++k)
        A[Rev(k,m)] = a[k];
    return A;
}
template<class DC> vector<complex<DC>> iterative_FFT(const vector<complex<DC>>&a,bool _jud)
{
    vector<complex<DC>> A = bit_reverse_copy(a);
    int n = a.size();
    int m= 1;
    complex<DC> omega_m;
    complex<DC> omega;
    for (int s = 1; s <= log2(n);++s)
    {
        m= m<< 1;
        if(_jud)
            omega_m=complex<DC>(cos(2*Pi<DC>/m),sin(2*Pi<DC>/m));
        else
            omega_m=complex<DC>(cos(2*Pi<DC>/m),0-sin(2*Pi<DC>/m));
        for (int k = 0; k < n;k=k+m)
        {
            omega =complex<DC>(1,0);
            int halfm = m >> 1;
            for (int j = 0; j < halfm;++j)
            {
                complex<DC> t=omega*A[k+j+halfm];
                complex<DC> u = A[k + j];
                A[k+j]=u+t;
                A[k + j + halfm] = u - t;
                omega = omega * omega_m;
            }
        }
    }
    return A;
}
template<class DC> vector<complex<DC>> recursive_FFT(const vector<complex<DC>>& a,bool _jud)
{
    int n = a.size();
    if(n==1)
       return a;
    int halfn = (n >> 1);
    complex<DC> omega_n;
    if(_jud)
        omega_n=complex<DC>(cos(2 * Pi<DC>/n),sin(2*Pi<DC>/n));
    else
        omega_n=complex<DC>(cos(2 * Pi<DC>/n),sin(0-2*Pi<DC>/n));
    complex<DC> omega{1, 0};
    vector<complex<DC>> a_even(halfn);
    vector<complex<DC>> a_odd(halfn);
    for (int i = 0;i<halfn;++i)
    {
        int twicei = i << 1;
        a_even[i] = a[twicei];
        a_odd[i] = a[twicei+1];
    }
    vector<complex<DC>> y_even= recursive_FFT(a_even,_jud);
    vector<complex<DC>> y_odd = recursive_FFT(a_odd,_jud);
    vector<complex<DC>> y(n);
    for (int k = 0; k < halfn;++k)
    {
        complex<DC> temp = omega * y_odd[k];
        y[k] = y_even[k] + temp;
        y[k + halfn] = y_even[k] - temp;
        omega = omega * omega_n;
    }
    return y;
}
template<class DC> vector<complex<DC>> DFT(const vector<complex<DC>> &a)
{
    return iterative_FFT(a, true);
}
template<class DC> vector<DC> Product(const vector<DC>& a,const vector<DC> &b)
{
    int n = a.size()<b.size()?a.size():b.size();
    vector<DC> c(n);
    for (int i = 0; i < n;++i)
        c[i] = a[i] * b[i];
    return c;
}
template<class DC> vector<complex<DC>> IDFT(const vector<complex<DC>>& a)
{
    unsigned n = a.size();
    vector<complex<DC>> r = iterative_FFT(a, false);
    for (unsigned i=0;i<n;++i)
        r[i] = r[i] /DC(n);
    return r;
}
template<class DC> vector<complex<DC>> FFT(const vector<complex<DC>> &x,const vector<complex<DC>> &y)
{
    int n = x.size();
    vector<complex<DC>> a=x;
    vector<complex<DC>> b=y;
    a.resize(n << 1);
    b.resize(n << 1);
    auto r=IDFT(Product(DFT(a), DFT(b)));
    r.pop_back();
    return r;
}
template<class DC> vector<DC> RealFFT(const vector<DC> &x,const vector<DC> &y)
{
    int n = x.size();
    vector<complex<DC>> a(n);
    vector<complex<DC>> b(n);
    for (int i = 0; i < n;++i)
    {
        a[i] = complex<DC>(x[i], 0);
        b[i] = complex<DC>(y[i], 0);
    }
    auto c = FFT(a, b);
    int m = c.size();
    vector<DC> r(m);
    for (int i = 0;i<m;++i)
        r[i]=c[i].real();
    return r;
}
template<class DC> vector<vector<DC>> _AddZero(const vector<DC> &x,const vector<DC> &y)
{
    uint32_t n = x.size() > y.size() ? x.size() : y.size();
    uint32_t nearest = 1 << uint8_t(log2(n));
    if(n!=nearest)
        nearest = nearest << 1;
    vector<DC> a=x;
    vector<DC> b=y;
    a.resize(nearest);
    b.resize(nearest);
    return {a, b};
}
template<class DC> vector<complex<DC>> Convolution(const vector<complex<DC>>&x,const vector<complex<DC>>& y)
{
    unsigned n = x.size() > y.size() ? x.size() : y.size();
    unsigned nearest = 1 << unsigned(log2(n));
    if(n!=nearest)
        nearest = nearest << 1;
    vector<complex<DC>> a=x;
    vector<complex<DC>> b=y;
    a.resize(nearest);
    b.resize(nearest);
    auto r=FFT(a, b);
    r.resize(x.size() + y.size() - 1);
    return r;
}
template<class DC> vector<DC> RealConvolution(const vector<DC>&x,const vector<DC>& y)
{
    vector<vector<DC>> ab = _AddZero(x, y);
    auto r=RealFFT(ab[0], ab[1]);
    r.resize(x.size() + y.size() - 1);
    return r;
}
template<class DC> DC ModPower(DC a,DC n,DC mod)
{
    DC b = 1;
    while(n >= 1)
    {
        if(n & 1)
        {
            b = (a * b) % mod;
        }
        a = (a * a) % mod;
        n >>= 1; 
    }
    return b;
}
template<class DC> vector<DC> ModProduct(const vector<DC>& a,const vector<DC> &b,DC mod)
{
    uint32_t n = a.size()<b.size()?a.size():b.size();
    vector<DC> c(n);
    for (uint32_t i = 0; i < n;++i)
        c[i] = a[i]*b[i]%mod;
    return c;
}
template<class DC> vector<DC> iterative_NTT(const vector<DC>&a,bool _jud)
{
    const DC g = 3; /* g is the minimal primitive root of mod */
    const DC gi = 332748118; /* gi is the modular inverse of g modulo mod */ /* 332748118 */ /* 1398101334 */
    const DC mod = 998244353; /* mod must be a prime */ /* 998244353 */ /* 4194304001 */
    vector<DC> A= bit_reverse_copy(a);
    uint32_t n = a.size();
    DC omega, omega_m;
    uint32_t halfm;
    for (uint32_t m=2; m<=n; m<<=1)
    {
        if(_jud)
            omega_m=ModPower(g,(mod-1)/DC(m),mod);
        else
            omega_m=ModPower(gi,(mod-1)/DC(m),mod);
        for (uint32_t k = 0; k < n;k=k+m)
        {
            omega =1;
            halfm = m >> 1;
            for (uint32_t j = 0; j < halfm;++j)
            {
                DC t= omega*A[k+j+halfm]%mod;
                DC u = A[k + j];
                A[k+j]=(u+t)%mod;
                A[k + j + halfm] =(mod+u - t)%mod;
                omega = omega * omega_m%mod;
            }
        }
    }
    if(!_jud)
    {
        DC n_inverse=ModPower(DC(n),mod-2,mod);
        for (uint32_t i = 0; i < n;++i)
            A[i]=A[i]*n_inverse%mod;//这里要乘长度的逆元
    }
    return A;
}
template<class DC> vector<DC> NTT(const vector<DC> &x, const vector<DC> &y)
{
    uint32_t n=x.size();
    vector<DC> a = x;
    vector<DC> b = y;
    a.resize(n << 1);
    b.resize(n << 1);
    const DC mod=998244353; /* previous mod = 998244353 */ /* 4194304001 */
    vector<DC> r=iterative_NTT(ModProduct(iterative_NTT(a, true), iterative_NTT(b, true),mod), false);
    r.pop_back();
    return r;
}
template<class DC> vector<DC> IntConvolution(const vector<DC> &x,const vector<DC> &y)
{
    vector<vector<DC>> ab = _AddZero(x, y);
    vector<DC> r=NTT(ab[0], ab[1]);
    r.resize(x.size()+y.size()-1);
    return r;
}
template<class DC> vector<DC> CarryBit(const vector<DC>& a,uint32_t bit)
{
    vector<DC> r=a;
    uint32_t n = r.size();
    for (uint32_t i = 0; i <n-1;++i)
    {
        if(r[i]>=bit)
        {
            r[i + 1] += r[i] / bit;
            r[i] = r[i] % bit;
        }
    }
    while(r[r.size()-1]>=bit)
    {
        r.emplace_back(r[r.size()-1]/bit);
        r[r.size() - 2] = r[r.size()-2] % bit;
    }
    return r;
}
template<class DC> string BitToString(const vector<DC> &a)
{
    //auto start = std::chrono::system_clock::now(); 
    /*uint32_t n = a.size();
    string str,temp;
    for (uint32_t i = 0; i < a.size();++i)
    {
        temp = std::to_string(a[i]);
        reverse(temp.begin(), temp.end());
        if(temp.size()<d)
            temp.append(d - temp.size(), '0');
        str.append(temp);
    }*/
    string str(a.size(),' ');
    std::transform(a.begin(), a.end(), str.begin(), [](DC c)
                   { return '0'+c; });
    return str;
}

}
#endif