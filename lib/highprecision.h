#ifndef HIGHPRECISION_H
#define HIGHPRECISION_H
#include "fft.h"
using std::cerr;
using std::cout;
using std::ostream;
using std::pair;
using std::to_string;
using namespace fft;

namespace hpc{

//-----声明部分-----

//大整数类
class BigInt
{
    public:
        BigInt cutzero();
        BigInt();
        BigInt(const string&);
        BigInt(const char *);
        BigInt(int);
        BigInt(bool);
        BigInt(const string &, bool);

        friend BigInt operator-(const BigInt &);
        friend bool operator==(const BigInt&, const BigInt&);
        friend bool operator!=(const BigInt &, const BigInt &);
        friend bool operator<(const BigInt&, const BigInt&);
        friend bool operator>(const BigInt&, const BigInt&);
        friend bool operator<=(const BigInt&, const BigInt&);
        friend bool operator>=(const BigInt&, const BigInt&);
        friend ostream &operator << (ostream &, const BigInt&);
        friend BigInt operator+(const BigInt &,const BigInt&);
        friend BigInt operator-(const BigInt &,const BigInt&);
        friend BigInt operator*(const BigInt &,const BigInt&);
        friend BigInt operator/(const BigInt&, const BigInt&);
        friend BigInt operator%(const BigInt &, const BigInt &);
        friend BigInt operator&(const BigInt &, const BigInt &);
        friend BigInt operator^(const BigInt&, int);
        friend BigInt operator^(const BigInt&,unsigned long long);
        BigInt &operator++();
        BigInt &operator+=(const BigInt &);
        BigInt &operator-=(const BigInt &);
        friend BigInt operator<<(const BigInt &,int);
        friend BigInt operator>>(const BigInt &,int);
        BigInt &operator<<=(int);
        BigInt &operator>>=(int);

        explicit operator int();
        explicit operator bool();
        friend BigInt LeftShift(const BigInt&, int);
        friend BigInt DivideTwo(const BigInt &);
        friend BigInt Product_NTT(const BigInt &, const BigInt &);
        friend BigInt Product_DivideConquer(const BigInt &, const BigInt &);
        friend pair<BigInt, BigInt> DivisionAlgorithm(const BigInt &, const BigInt &);
        void Reassign(string );
        void _digitChange(string);
        void _signChange(bool);
        friend string Get_digit(const BigInt&);
        friend string GetString(const BigInt &);
        friend bool PositivityTest(const BigInt&);
        friend unsigned IntegerLength(const BigInt &);
        friend BigInt abs(const BigInt &);
        friend BigInt sqrt(const BigInt &);
        friend unsigned log2(const BigInt &);

    private:
        string _digit;
        bool _sign;
};

class HighPrecision: public BigInt
{
    public:
        HighPrecision CutTail();
        HighPrecision();
        HighPrecision(string);
        HighPrecision(const char[]);
        explicit HighPrecision(const BigInt&);
        HighPrecision(double);
        HighPrecision(string, int);
        HighPrecision(double, int);

        friend ostream &operator<<(ostream &, const HighPrecision &);
        friend HighPrecision operator+(const HighPrecision &, const HighPrecision &);
        friend HighPrecision operator- (const HighPrecision &);
        friend HighPrecision operator-(const HighPrecision &, const HighPrecision &);
        friend HighPrecision operator*(const HighPrecision &, const HighPrecision &);
        friend HighPrecision operator/(const HighPrecision &, const HighPrecision &);

        friend unsigned DecimalLength(const HighPrecision &);
        friend BigInt IntegerPart(const HighPrecision &);
        friend HighPrecision DecimalPart(const HighPrecision &);
        friend int Order(const HighPrecision &);
        friend string TopKDigit(const HighPrecision &,unsigned);
        friend unsigned TotalLength(const HighPrecision &);
        friend unsigned SignificantLength(const HighPrecision &);
        friend HighPrecision SetSignificantFigure(const HighPrecision &, unsigned, bool);
        friend HighPrecision SetSignificantFigure(const HighPrecision &, unsigned);
        friend HighPrecision LeftShift(const HighPrecision&, int);
        friend HighPrecision Divide(const HighPrecision &, const HighPrecision &, uint8_t, uint32_t);

    private:
        string _decimal;
};

//杂例II
BigInt Factorial(int);
BigInt Fibonacci_r(int);
BigInt Fibonacci_m(int);
BigInt Fibonacci(int);

//

//大整数类
BigInt BigInt::cutzero()
{
    int i = _digit.size()-1;
    while(i>=0 && _digit[i]=='0')
    {
        _digit.pop_back();
        --i;
    }
    if(i<0)
        _digit = "0";
    return *this;
}
BigInt::BigInt()
{
    _digit= "0";
    _sign = true;
}
BigInt::BigInt(const string &a)
{
    _digit = a;
    reverse(_digit.begin(), _digit.end());
    if(_digit.back()=='-')
    {
        _sign = false;
        _digit.pop_back();
    }
    else
        _sign = true;
    cutzero();
}
BigInt::BigInt(const char * a)
{
    _digit = a;
    reverse(_digit.begin(), _digit.end());
    if(_digit.back()=='-')
    {
        _sign = false;
        _digit.pop_back();
    }
    else
        _sign = true;
    cutzero();
}
BigInt::BigInt(int a)
{
    _digit=to_string(a);
    reverse(_digit.begin(), _digit.end());
    if(a<0)
    {
        _sign = false;
        _digit.pop_back();
    }
    else
        _sign = true;
}
BigInt::BigInt(bool flag)
{
    _digit = "1";
    if(flag)   
        _sign = true;
    else
        _sign = false;
}
BigInt::BigInt(const string &a, bool flag)
{
    _digit = a;
    _sign = flag;
    cutzero();
}
BigInt operator-(const BigInt &a)
{
    BigInt b = a;
    b._sign = !a._sign;
    return b;
}
bool operator==(const BigInt &a, const BigInt &b)
{
    return (a._digit == b._digit) && (a._sign==b._sign);
}
bool operator!=(const BigInt &a, const BigInt &b)
{
    return !(a == b);
}
bool operator<(const BigInt &a,const BigInt &b)
{
    if(a._sign!=b._sign)
        return b._sign;
    bool sign = a._sign;
    if(a._digit.size()<b._digit.size())
        return sign;
    if(a._digit.size()>b._digit.size())
        return !sign;
    int n = a._digit.size();
    for (int i = n - 1; i >= 0;--i)
    {
        if(a._digit[i]-'0'<b._digit[i ]-'0')
            return sign;
        if(a._digit[i]-'0'>b._digit[i]-'0')
            return !sign;
    }
    return !sign;
}
bool operator>(const BigInt &a,const BigInt &b)
{
    return b < a;
}
bool operator<=(const BigInt &a, const BigInt &b)
{
    return (a == b) || (a < b);
}
bool operator>=(const BigInt &a, const BigInt &b)
{
    return (a == b) || (a > b);
}
ostream & operator<<(ostream & os,const BigInt &a)
{
    os << GetString(a);
    return os;
}
BigInt LeftShift(const BigInt &a, int b)
{
    BigInt c = a;
    c._digit.insert(0,b, '0');
    return c;
}
BigInt operator+(const BigInt &a,const BigInt &b)
{
    if(a._sign==false && b._sign==true)
        return b - (-a);
    if(a._sign==true && b._sign==false)
        return a - (-b);
    int m = a._digit.size();
    int n = b._digit.size();
    if(n>m)
        return b + a;
    int flag = 0;
    int s = 0;
    int i = 0;
    BigInt r;
    r._digit = "";
    for (; i < n;++i)
    {
        s = a._digit[i]-'0' + b._digit[i]-'0' + flag;
        r._digit.push_back((s % 10)+'0');
        if(s>=10)
            flag = 1;
        else
            flag = 0;
    }
    for (; i < m;++i)
    {
        s =a._digit[i]-'0' + flag;
        r._digit.push_back((s % 10)+'0');
        if(s>=10)
            flag = 1;
        else
            flag = 0;
    }
    if (flag)
        r._digit.push_back(flag + '0');
    r.cutzero();
    r._sign = a._sign;
    return r;
}
BigInt operator-(const BigInt &a,const BigInt &b)
{
    if(a._sign!=b._sign)
        return a+(-b);
    if(a._sign==false)
        return -((-a)-(-b));
    if(a<b)
        return -(b - a);
    int m = a._digit.size();
    int n = b._digit.size();
    int flag = 0;
    int s = 0;
    int i = 0;
    BigInt r;
    r._digit = "";
    for (; i < n;++i)
    {
        s = a._digit[i]-'0' - (b._digit[i]-'0')+10-flag;
        r._digit.push_back((s % 10)+'0');
        if(s<10)
            flag = 1;
        else
            flag = 0;
    }
    for (; i < m;++i)
    {
        s =a._digit[i]-'0'+10- flag;
        r._digit.push_back((s % 10)+'0');
        if(s<10)
            flag = 1;
        else
            flag = 0;
    }
    if (flag)
        r._digit.push_back(flag + '0');
    r.cutzero();
    return r;
}
BigInt DivideTwo(const BigInt &x)
{
    BigInt r;
    r._sign = x._sign;
    size_t m = x._digit.size();
    string new_digit(m,'0');
    int8_t temp=0;
    for (size_t i = m; i>0;--i)
    {
        int8_t current = int8_t(x._digit[i - 1] - '0') + temp*10;
        new_digit[i - 1] = (current>>1) + '0';
        temp = current & 1;
    }
    r._digit = new_digit;
    r.cutzero();
    if(x._sign==false && (int(x._digit[0] - '0') & 1))
        return r - BigInt("1");
    return r;
}
BigInt Product_DivideConquer(const BigInt &x,const BigInt &y)
{
    BigInt p = x;
    BigInt q = y;
    if(p._digit.size()<q._digit.size())
        p._digit.append(q._digit.size()-p._digit.size(),'0');
    else if(p._digit.size()>q._digit.size())
        q._digit.append(p._digit.size() - q._digit.size(), '0');
    int m=p._digit.size();
    int n =q._digit.size();
    int k = (m >> 1);
    int u = (n >> 1);
    if(m<=8 && n<=8)
    {
        BigInt r;
        string a = p._digit;
        string b = q._digit;
        reverse(a.begin(),a.end());
        reverse(b.begin(), b.end());
        b=to_string(stoull(a) * stoull(b));
        reverse(b.begin(),b.end());
        r._digit = b;
        r._sign = (x._sign == y._sign);
        return r;
    }
    BigInt p0, p1, q0, q1;
    p0._digit=p._digit.substr(0,k);
    p1._digit=p._digit.substr(k, m - k);
    q0._digit=q._digit.substr(0,u);
    q1._digit=q._digit.substr(u,n-u);
    BigInt r0= Product_DivideConquer(p0,q0);
    BigInt r1 = Product_DivideConquer(p1,q1);
    BigInt r2 = Product_DivideConquer((p0 + p1), (q0 + q1));
    BigInt r=r0 + LeftShift((r2 - r0 - r1), k) + LeftShift(r1, (k << 1));
    r._sign = (x._sign == y._sign);
    r.cutzero();
    return r;
}
vector<unsigned long long> CompressBit(const string &a_str)
{
    vector<unsigned long long> a(a_str.size());
    std::transform(a_str.begin(), a_str.end(), a.begin(), [](char c)
                   { return c - '0'; });
    return a;
}
BigInt Product_NTT(const BigInt &x,const BigInt &y)
{
    auto rr =IntConvolution(CompressBit(x._digit), CompressBit(y._digit));
    auto r = CarryBit<unsigned long long>(rr, 10);
    return BigInt(BitToString(r), (x._sign == y._sign));
}
BigInt operator*(const BigInt &x,const BigInt &y)
{
    uint32_t n=std::max((x._digit).size(),(y._digit).size());
    if (n<=15)
        return Product_DivideConquer(x, y);
    return Product_NTT(x, y);
}
pair<BigInt,BigInt> DivisionAlgorithm(const BigInt &a,const BigInt &b)
{
    const uint32_t snf = a._digit.size()>b._digit.size()? 
        (a._digit.size() -b._digit.size() + 1):1;
    const uint8_t times = std::max(0., std::ceil(std::log2(snf)) - 2);
    BigInt q = Divide(HighPrecision(a), HighPrecision(b), times, snf);
    BigInt r = a - q * b;
    if(r<0)
    {
        q -= BigInt("1");
        r += b;
    }
    if(r>=b && PositivityTest(b))
    {
        q += BigInt("1");
        r -= b;
    }
    return {q,r};
}
BigInt operator/(const BigInt &a,const BigInt &b)
{
    if(b==BigInt(0))
    {
        cerr << "错误：不能除以0." << '\n';
        return BigInt(0);
    }
    if(a==BigInt(0))
        return BigInt(0);
    if(b==BigInt("1"))
        return a;
    if(b==BigInt("2"))
        return DivideTwo(a);
    return std::get<0>(DivisionAlgorithm(a,b));
}
BigInt operator%(const BigInt &a, const BigInt &b)
{
    if(b==BigInt(0))
    {
        cerr << "错误：不能除以0." << '\n';
        return BigInt(0);
    }
    if(a==BigInt(0))
        return BigInt(0);
    if(b==BigInt("1"))
        return BigInt(0);
    if(b==BigInt("2"))
        return BigInt(int(a._digit[0] - '0') & 1);
    if(b==BigInt("10"))
        return BigInt(int(a._digit[0] - '0'));
    if(b==BigInt("5") && PositivityTest(a))
        return BigInt(int(a._digit[0] - '0') % 5);
    if(b==BigInt("3") && PositivityTest(a))
    {
        int r = 0;
        for(auto i:a._digit)
            r = (r + (i - '0')) % 3;
        return BigInt(r);
    }
    if(b==BigInt("9") && PositivityTest(a))
    {
        int r = 0;
        for(auto i:a._digit)
            r = (r + (i - '0')) % 9;
        return BigInt(r);
    }
    return std::get<1>(DivisionAlgorithm(a,b));
}
BigInt operator&(const BigInt &a, const BigInt &b)
{
    if(b==BigInt(0) || a==BigInt(0))
        return BigInt(0);
    if(b==BigInt("1"))
        return BigInt(int(a._digit[0] - '0') & 1);
    if(a==BigInt("1"))
        return BigInt(int(b._digit[0] - '0') & 1);
    BigInt r;
    size_t t = std::min(a._digit.size(), b._digit.size());
    for (size_t i = 0; i < t;++i)
    {
        r._digit[i] = ((a._digit[i] - '0') & (b._digit[i] - '0')) + '0';
    }
    r._sign = a._sign & b._sign;
    return r;
}
BigInt operator^(const BigInt &a, int n)
{
    if(n<0)
    {
        cerr << "错误：尚不支持大整数的负数次幂。" << '\n';
        return a;
    }
    BigInt m=a;
    BigInt b("1");
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
BigInt operator^(const BigInt &a,unsigned long long n)
{
    if(n<0)
    {
        cerr << "错误：尚不支持大整数的负数次幂。" << '\n';
        return a;
    }
    BigInt m=a;
    BigInt b("1");
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
BigInt &BigInt::operator++()
{
    *this = *this + BigInt("1");
    return *this;
}
BigInt &BigInt::operator+=(const BigInt &b)
{
    *this = *this + b;
    return *this;
}
BigInt &BigInt::operator-=(const BigInt &b)
{
    *this = *this -b;
    return *this;
}
BigInt operator<<(const BigInt &a,int n)
{
    if (n==0)
        return a;
    if (n==1)
        return (a * BigInt("2"));
    if (n>1)
        return a *(BigInt("2") ^ n);
    return a >> (-n);
}
BigInt operator>>(const BigInt &a, int n)
{
    if (n==0)
        return a;
    if (n==1)
        return DivideTwo(a);
    if (n>0)
    {
        BigInt b = a;
        for (int i = 0; i < n;++i)
            b = DivideTwo(b);
        return b;
    }
    return a << (-n);
}
BigInt &BigInt::operator<<=(int n)
{
    *this = *this << n;
    return *this;
}
BigInt &BigInt::operator>>=(int n)
{
    *this = *this >> n;
    return *this;
}
BigInt::operator int()
{
    string _str = _digit;
    reverse(_str.begin(),_str.end());
    if(_sign)
        return stoi(_str);
    return -stoi(_str);
}
BigInt::operator bool()
{
    if (_digit == "0")
        return false;
    return true;
}
void BigInt::Reassign(string str)
{
    *this = BigInt(str);
}
void BigInt::_digitChange(string str)
{
    _digit = str;
}
void BigInt::_signChange(bool sign)
{
    _sign = sign;
}
string Get_digit(const BigInt&a)
{
    return a._digit;
}
string GetString(const BigInt&a)
{
    string _string=a._digit;
    if(a._sign==false)
        _string.push_back('-');
    reverse(_string.begin(), _string.end());
    return _string;
}
bool PositivityTest(const BigInt &a)
{
    return a._sign;
}
unsigned IntegerLength(const BigInt &a)
{
    if(a==BigInt(0))
        return 0;
    return a._digit.size();
}
BigInt abs(const BigInt &c)
{
    return BigInt(c._digit, true);
}
BigInt sqrt(const BigInt &c)
{
    if(!PositivityTest(c))
    {
        std::cerr << "错误：根号下不能是负数。" << std::endl;
        return 0;
    }
    if(c< BigInt("10"))
    {
        BigInt c_copy = c;
        return BigInt(int(std::sqrt(int(c_copy))));
    }
    BigInt pre("0");
    BigInt prepre("0");
    const size_t length = c._digit.size();
    string str((length + 1) / 2, '0');
    if(length & 1)
    {
        str[0] = '0' + int(std::sqrt(c._digit[length - 1]-'0'));
    }
    else
    {
        str[0] = '0' + int(std::sqrt((c._digit[length - 1] - '0') * 10 + (c._digit[length - 2] - '0')));
    }
    BigInt x(str);
    while (x!=pre && x!=prepre)
    {
        prepre = pre;
        pre = x;
        x = DivideTwo(x + (c / x));
    }
    if(x!=pre) // this case happens only if (c+1) is a square
        return pre > x ? x : pre;
    return pre;
}
unsigned log2(const BigInt &c)
{
    if(!PositivityTest(c) || c==BigInt("0"))
    {
        std::cerr << "错误：只有正数才能取对数。" << std::endl;
        return 0;
    }
    unsigned i = 0;
    BigInt a = c;
    while(a)
    {
        a = DivideTwo(a);
        ++i;
    }
    return (i - 1);
}

HighPrecision HighPrecision::CutTail()
{
    const string &str = _decimal;
    int d =-1;
    int n = str.size();
    for (int i = 0; i < n;++i)
    {
        if(str[i]!='0')
        {
            d=i;
            break;
        }
    }
    if(d>0)
        _decimal.erase(0,d);
    else if(d<0)
        _decimal = "";
    return *this;
}
HighPrecision::HighPrecision()
{
    _decimal = "";
}
HighPrecision::HighPrecision(string str)
{
    unsigned d = str.size();
    for (unsigned i = 0; i <str.size();++i)
        if(str[i]=='.')
        {
            d = i;
            break;
        }
    Reassign(str.substr(0,d));
    if(d>=str.size()-1)   
        _decimal = "";
    else
    {
        _decimal = str.substr(d + 1);
        reverse(_decimal.begin(), _decimal.end());
    }
    CutTail();
}
HighPrecision::HighPrecision(const char a[])
{
    new (this) HighPrecision(string(a));
}
HighPrecision::HighPrecision(const BigInt &a)
{
    _digitChange(Get_digit(a));
    _signChange(PositivityTest(a));
    _decimal = "";
}
HighPrecision::HighPrecision(double a)
{
    new (this) HighPrecision(to_string(a));
}
HighPrecision::HighPrecision(string a,int n)
{
   *this=LeftShift((HighPrecision(a)),n);
}
HighPrecision::HighPrecision(double a,int n)
{
    new (this) HighPrecision(to_string(a), n);
}
ostream &operator<<(ostream &os, const HighPrecision &a)
{ 
    string result = GetString(a);
    if(a._decimal.size()>0)
    {
        string _string = a._decimal;
        reverse(_string.begin(), _string.end());
        result.append(".");
        result.append(_string);
    }
    os << result;
    return os;
}
HighPrecision LeftShift(const HighPrecision &a,int n)
{
    if(n==0)
        return a;
    HighPrecision t = a;
    if(n>0)
    {
        const string &cat = t._decimal + Get_digit(t);
        int len = t._decimal.size();
        if(n>len)
        {
            string zeros;
            zeros.insert(zeros.end(), n - len, '0');
            t._decimal = "";
            t._digitChange(zeros + cat);
        }
        else if(n==len)
        {
            t._decimal ="";
            t._digitChange(cat);
        }
        else
        {
            t._decimal = cat.substr(0, len - n);
            t._digitChange(cat.substr(len - n));
        }
        t.cutzero();
        return t;
    }
    int len= Get_digit(t).size();
    string cat=t._decimal+Get_digit(t);
    cat.insert(cat.end(),-n,'0');
    t._digitChange(cat.substr(cat.size() - len));
    t._decimal=cat.substr(0, cat.size() - len);
    t.cutzero();
    t.CutTail();
    return t;
}
HighPrecision operator+(const HighPrecision &a, const HighPrecision &b)
{
    string str1=a._decimal;
    string str2=b._decimal;
    if(str1.size()<str2.size())
        str1.insert(str1.begin(),str2.size()-str1.size(),'0');
    else if(str2.size()<str1.size())
        str2.insert(str2.begin(),str1.size()-str2.size(),'0');
    unsigned dsize= str1.size();
    str1.append(Get_digit(a));
    str2.append(Get_digit(b));
    BigInt c, d;
    c._digitChange(str1);
    d._digitChange(str2);
    c._signChange(PositivityTest(a));
    d._signChange(PositivityTest(b));
    HighPrecision r(c+d);
    const string &str = Get_digit(r);
    if(dsize<str.size())
    {
        r._decimal = str.substr(0,dsize);
        r._digitChange(str.substr(dsize));
    }
    else
    {
        r._decimal = str;
        r._decimal.insert(r._decimal.end(), dsize - str.size(), '0');    
        r._digitChange("0");
    }
    r.CutTail();
    return r;
}
HighPrecision operator-(const HighPrecision &a)
{
    HighPrecision b = a;
    b._signChange(!PositivityTest(a));
    return b;
}
HighPrecision operator-(const HighPrecision &a,const HighPrecision &b)
{
    return a +(-b);
}
HighPrecision operator*(const HighPrecision &a, const HighPrecision &b)
{
    const string &str1 = a._decimal + Get_digit(a);
    const string &str2 = b._decimal + Get_digit(b);
    BigInt c, d;
    c._digitChange(str1);
    d._digitChange(str2);
    c._signChange(PositivityTest(a));
    d._signChange(PositivityTest(b));
    HighPrecision r(c*d);
    const string &str=Get_digit(r);
    unsigned dsize=a._decimal.size() + b._decimal.size();
    if(dsize<str.size())
    {
        r._decimal = str.substr(0,dsize);
        r._digitChange(str.substr(dsize));
    }
    else
    {
        r._decimal = str;
        r._decimal.insert(r._decimal.end(), dsize - str.size(), '0');    
        r._digitChange("0");
    }
    r.CutTail();
    return r;
}
HighPrecision Divide(const HighPrecision &a, const HighPrecision &b, uint8_t times, uint32_t snf)
{
    if(BigInt(b)==0 && b._decimal.size()==0)
    {
        cerr << "错误：不能除以0" << '\n';
        return HighPrecision("0");
    }  
    const HighPrecision two("2");
    int order_b = Order(b); // b = (...something lying in [1,10)...) * 10 ^ (order_b)
    // shift b to B = b * 10 ^ (order_b), where B lies in [1,10)
    const HighPrecision B = LeftShift(b,(-order_b)); 
    // Next we apply Newton's method to compute 1/B.
    // Note that TopKDigit(b,5) lies in [10^4,10^5).
    /* So 10^4/TopKDigit(b,5) is a good apporixmate of 1/B, 
        with accuracy of at least 4 significant figures. */
    // We thus set it as the inital value.
    uint32_t accuracy = 4;
    HighPrecision x(10000. / stod(TopKDigit(b, 5)));
    // Newton's method
    for (uint8_t i = 0; i < times; ++i)
    {
        x = SetSignificantFigure(x, accuracy + 1, false); // cut off unnecessary decimal digits
        x=x*(two-x*B);
        accuracy <<= 1; // under each iteartion, the accuracy of siginificant figures is doubled.
    }
    // After the above iteration, x = 1/B, with accuracy of 2 ^ (times + 2) significant figures
    x = SetSignificantFigure(x, std::min(snf, accuracy) + 1, false); 
    x = LeftShift((a*x),(-order_b)); // As a result, a / b = A * (1/B) * 10^ (-order_b) 
    x._signChange(!(PositivityTest(a) ^ PositivityTest(b)));
    x.CutTail();
    return x;
}
HighPrecision NDivide(const HighPrecision &a, const HighPrecision &b, uint32_t snf)
{
    if(snf == 0)
        return HighPrecision("0");
    const uint8_t times = std::max(0., std::ceil(std::log2(snf)) - 2);
    return SetSignificantFigure(Divide(a, b, times, snf), snf, true);
}
HighPrecision operator/(const HighPrecision &a, const HighPrecision &b)
{
    const uint32_t snf = std::max(8u,
        std::max(SignificantLength(a), SignificantLength(b)));
    return NDivide(a, b, snf);
}
unsigned DecimalLength(const HighPrecision &a)
{
    return a._decimal.size();
}
BigInt IntegerPart(const HighPrecision &a)
{
    BigInt r;
    r._digitChange(Get_digit(a));
    r._signChange(PositivityTest(a));
    return r;
}
HighPrecision DecimalPart(const HighPrecision &a)
{
    HighPrecision r=a;
    r._digitChange("0");
    return r;
}
int Order(const HighPrecision &a)
{
    const BigInt &intpart=BigInt(a);
    if(intpart!=BigInt("0"))
        return IntegerLength(a)-1;
    return -(a._decimal.size());
}
string TopKDigit(const HighPrecision &a,unsigned k)
{
    string result;
    if(BigInt(a)==BigInt(0) && a._decimal.size()==0)
    {
        result.insert(result.end(), k, '0');
        return result;
    }
    if(!PositivityTest(a))
        result.push_back('-');
    const string &cat = a._decimal + Get_digit(a);
    unsigned n = cat.size();
    unsigned m = std::min(k, n);
    for (unsigned i = 1; i <= m; ++i)
        result.push_back(cat[n-i]);
    if(n<k)
        result.insert(result.end(), k - n, '0');
    return result;
}
unsigned TotalLength(const HighPrecision &a)
{
    if(BigInt(a)==BigInt(0) && a._decimal.size()==0)
        return 0;
    return IntegerLength(a) + a._decimal.size();
}
unsigned SignificantLength(const HighPrecision & a)
{
    const string &str = a._decimal + Get_digit(a);
    int n = str.size();
    int lft=-1, rght=-1;
    for (int i = 0; i < n;++i)
        if(str[i]!='0')
        {
            lft = i;
            break;
        }
    if(lft<0)
        return 0;
    for (int i = n-1; i>=0;--i)
        if(str[i]!='0')
        {
            rght = i;
            break;
        }
    return rght - lft + 1;
}
HighPrecision SetSignificantFigure(const HighPrecision &a, unsigned snf, bool round)
{
    if(snf == 0)
        return HighPrecision("0");
    string str=a._decimal+Get_digit(a);
    const unsigned len= str.size();
    unsigned d=len;
    for (int i = len - 1; i >= 0;--i) // find the highest non-zero position
        if(str[i]!='0')
        {
            d=i;
            break;
        }
    if(d==len)
        return HighPrecision("0");
    if(d<snf)
        return a;
    unsigned d_minus_snf_plus_one = d - snf + 1;
    if(round && str[d-snf]>='5') // 四舍五入
    {
        unsigned i = d_minus_snf_plus_one;
        while (i <= d)
        {
            if(str[i]<'9')
            {
                str[i] = str[i] + 1;
                break;
            }
            str[i] = '0';
            ++i;
        }
        if(i == d+1)
        {
            if(i<len) 
                str[i] = '1';
            else
                str.push_back('1'); 
        }
    }
    unsigned decilen = a._decimal.size();
    HighPrecision r;
    if(d_minus_snf_plus_one<decilen)
        r._decimal=str.substr(d_minus_snf_plus_one,decilen-d_minus_snf_plus_one);
    for (unsigned i = decilen; i < d_minus_snf_plus_one; ++i)
        str[i] = '0';
    r._digitChange(str.substr(decilen));
    r._signChange(PositivityTest(a));
    r.CutTail();
    return r;
}
HighPrecision SetSignificantFigure(const HighPrecision &a, unsigned snf)
{
    return SetSignificantFigure(a, snf, true);
}

//杂例II
BigInt Factorial(int n)
{
    if(n<0)
    {
        cerr<<"错误：负数没有阶乘。"<<'\n';
        return 0;
    }
    BigInt s("1");
    for (int i = 1; i <=n;++i)
    {
        s =s* BigInt(i);
    }
    return s;
}
BigInt Fibonacci_r(int n)
{
    if(n==0)
        return 0;
    if(n==1)
        return 1;
    BigInt pre2(0);
    BigInt pre1(1);
    BigInt s(1);
    for (int i = 2; i <=n;++i)
    {
        s = pre2 + pre1;
        pre2 = pre1;
        pre1 = s;
    }
    return s;
}
BigInt Fibonacci_m(int n,vector<BigInt> &g)
{
    if(g[n]>0)
        return g[n];
    if(n%2)
        return g[n]=Fibonacci_m(n/2,g)*Fibonacci_m(n/2,g)+Fibonacci_m(n/2+1,g)*Fibonacci_m(n/2+1,g);
    return g[n] = (Fibonacci_m(n / 2 - 1, g) + Fibonacci_m(n / 2 + 1, g)) * Fibonacci_m(n / 2, g);
}
BigInt Fibonacci(int n)
{
    if(n<255)
        return Fibonacci_r(n);
    vector<BigInt> g({0, 1, 1,2,3});
    g.resize(n+1);
    return Fibonacci_m(n, g);
}

}

#endif