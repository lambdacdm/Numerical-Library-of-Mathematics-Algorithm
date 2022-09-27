#ifndef NUMERICAL_H
#define NUMERICAL_H
#include "polynomial.h"
#include "linearalgebra.h"
#include<deque>
#include<tuple>
using std::deque;
using std::tuple;
using namespace ply;
using namespace lag;

namespace nmr{

//-----声明部分-----

//插值
template <class DB> vector<std::function<DB(DB)>> LagrangeBasis(const Matrix<DB> &);
template <class DB> DB LBasisToInterp(const vector<std::function<DB(DB)>>&,const Matrix<DB> &,DB);
template <class DB> Matrix<DB> DifferenceQuotientTable(const Matrix<DB> &);
template <class DB> Matrix<DB> DividedDifferenceTable(const Matrix<DB> &);
template <class DB> DB DifferenceQuotient(const Matrix<DB> &);
template <class DB> DB DifferenceQuotient(std::function<DB(DB)>, const Matrix<DB> &);
template <class DB> DB DQTableToInterp(const Matrix<DB> &, const Matrix<DB> &,DB);
template <class DB> DB DDTableToInterp(const Matrix<DB> &, const Matrix<DB> &,DB);
template <class DB> DB _2P3DHermiteInterp(const Matrix<DB> &,int i,int j,DB);
template <class DB> DB _3P3DHermiteInterp(const Matrix<DB> &,DB);
template <class DE> Matrix<DE> SplineSlope(const Matrix<DE> &);
template <class DE> DE Interpolation(const Matrix<DE> &,DE,string);
template <class DE> DE Interpolation(const Matrix<DE> &,DE);
template <class DE> Matrix<DE> Interpolation(const Matrix<DE> &A, const Matrix<DE> &x, string);
template <class DE> std::function<DE(DE)> Interpolation(const Matrix<DE> &,string);

//拟合与逼近
template <class DB> Matrix<DB> PseudoInverse(const Matrix<DB>&);
template <class DB> Matrix<DB> LeastSquares(const Matrix<DB>&, const Matrix<DB>&);
template <class DB> Polynomial<DB> PolyFit(const Matrix<DB> &, const vector<Polynomial<DB>> &);
template <class DB> Polynomial<DB> PolyFit(const Matrix<DB>&,int);
template <class DB> Polynomial<DB> PolyFit(std::function<DB(DB)>,const vector<Polynomial<DB>>&,DB,DB);

//积分
template <class DI> DI Integrate(std::function<DI(DI)>,DI,DI,string);
template <class DI> DI Integrate(std::function<DI(DI)>,DI,DI);
template <class FS> FS FourierSinCoefficient(const std::function<FS(FS)> &, int32_t, FS);
template <class FS> FS FourierCosCoefficient(const std::function<FS(FS)> &, int32_t, FS);

//微分
template <class DD> DD D(std::function<DD(DD)>, DD);

//方程(组)求根
template<class DD> DD FindRoot(std::function<DD(DD)>, DD,DD);
template <class DD> DD FindRoot(std::function<DD(DD)>, DD,string);
template <class DD> DD FindRoot(std::function<DD(DD)>, DD);
template <class DD> DD FindRoot(std::function<DD(DD)>, DD,string,int);


//-----定义部分-----

//插值
template <class DB> vector<std::function<DB(DB)>> LagrangeBasis(const Matrix<DB> &A)
{
    int n = RowSize(A);
    vector<std::function<DB(DB)>> r;
    std::function<DB(DB)> f;
    for (int i = 0; i <n;++i)
    {
        f = [=](DB x) {
            DB s1=DB(1);
            DB s2=DB(1);
            for (int j = 0; j < i;++j)
            {
                s1=s1*(x-A(j,0));
                s2=s2*(A(i,0) - A(j,0));
            }
            for (int j = i + 1;j<n;++j)
            {
                s1=s1*(x-A(j,0));
                s2=s2*(A(i,0) - A(j,0));
            }
            return s1 / s2;
        };
        r.push_back(f);
    }
    return r;
}
template <class DB> DB LBasisToInterp(const vector<std::function<DB(DB)>> &lb,const Matrix<DB> &A,DB x)
{
    int n = RowSize(A);
    DB s = lb[0](x) * A(0,1);
    for (int i = 1; i < n;++i)
    {
        s += lb[i](x) * A(i,1);
    }
    return s;
}
template <class DB> Matrix<DB> DifferenceQuotientTable(const Matrix<DB> & A)
{
    int n = RowSize(A);
    Matrix<DB> table(n,n);
    for (int i = 0; i < n;++i)
        table(i,0) = A(i,1);
    for (int j = 1; j <n;++j)
    {
        for (int i = j; i < n;++i)
        {
            table(i,j) =(table(i,j-1)-table(i-1,j-1))/(A(i,0)-A(i-j,0));
        }
    }
    return table;
}
template <class DB> Matrix<DB> DividedDifferenceTable(const Matrix<DB> &A)
{
    int data_num= RowSize(A);
    int der_deg = ColumnSize(A)-1;
    int n = data_num * der_deg;
    int factorial = 1;
    vector<DB> z(n);
    for (int i = 0; i <n;++i)
        z[i] = A(i/der_deg,0);
    Matrix<DB> table(n, n);
    for (int j = 0; j <der_deg;++j)
    {
        for (int i = j; i < n;++i)
        {
            if((i%der_deg)<j)
                table(i,j)=(table(i,j-1)-table(i-1,j-1))/(z[i]-z[i-j]);
            else
                table(i,j) = A(i/der_deg, j+1)/factorial;
        }
        factorial = factorial * (j + 1);
    }
    for(int j=der_deg;j<n;++j)
        for (int i = j; i < n;++i)
            table(i,j)=(table(i,j-1)-table(i-1,j-1))/(z[i]-z[i-j]);
    return table;
}
template <class DB> DB DifferenceQuotient(const Matrix<DB> &A)
{
    int n = RowSize(A);
    return (DifferenceQuotientTable(A))(n-1,n-1);
}
template <class DB> DB DifferenceQuotient(std::function<DB(DB)> f, const Matrix<DB> &x)
{
    int n = RowSize(x);
    Matrix<DB> A(n,2);
    for (int i = 0; i < n;++i)
    {
        A(i,0) = x(i,0);
        A(i,1) = f(x(i,0));
    }
    return DifferenceQuotient(A);
}
template <class DB> DB DQTableToInterp(const Matrix<DB> &A, const Matrix<DB> &table,DB x)
{
    int n = RowSize(A);
    DB s = table(n-1,n-1);
    for (int i = n - 2; i >= 0;--i)
    {
        s = table(i,i) + (x -A(i,0)) * s;
    }
    return s;
}
template <class DB> DB DDTableToInterp(const Matrix<DB> &A, const Matrix<DB> &table,DB x)
{
    int data_num = RowSize(A);
    int der_deg=ColumnSize(A)-1;
    int n = data_num * der_deg;
    vector<DB> nest(n);
    nest[0] = DB(1);
    for(int i=0;i<n-1;++i)
        nest[i+1] = nest[i] * (x - A(i/der_deg,0));
    DB s=DB(0);
    for (int i = 0;i<n;++i)
        s += nest[i] * table(i,i);
    return s;
}
template <class DE> DE _2P3DHermiteInterp(const Matrix<DE> &M,int i,int j,DE x)
{
    DE l0x = (x-M(j,0)) / (M(i,0) - M(j,0));
    DE l1x = (x - M(i,0)) / (M(j,0) - M(i,0));
    DE l0x2 = l0x * l0x;
    DE l1x2 = l1x * l1x;
    return M(i,1) * (1 + 2 * l1x) * l0x2 + M(j,1) * (1 + 2 * l0x) * l1x2 +
           M(i,2) * (x - M(i,0)) * l0x2 + M(j,2) * (x - M(j,0)) *l1x2;
}
template <class DE> DE _3P3DHermiteInterp(const Matrix<DE> &M,DE x)
{
    std::function<DE(DE)> L2 = [=](DE t) {
        return DQTableToInterp(M, DifferenceQuotientTable(M),t);
    };
    DE k = (M(1,2)-D(L2,M(1,0)))/
           ((M(1,0)-M(0,0))*(M(1,0)-M(2,0)));
    return L2(x) + k * (x - M(0,0)) * (x - M(1,0)) * (x - M(2,0));
}
template <class DE> Matrix<DE> SplineSlope(const Matrix<DE> &A)
{
    int n = RowSize(A);
    vector<DE> h(n);
    vector<DE> lambda(n);
    vector<DE> mu(n);
    Matrix<DE> g(n, 1);
    h[0] = A(1,0) - A(0,0);
    lambda[0]=DE(0);
    mu[0] = DE(1);
    g(0,0) = 3 * (A(1,1)-A(0,1)) / h[0];
    for (int k = 1; k<= n - 2;++k)
    {
        h[k]=A(k+1,0)-A(k,0);
        lambda[k] = h[k] / (h[k] + h[k - 1]);
        mu[k]=1-lambda[k];
        g(k,0)=3*(lambda[k]*(A(k,1)-A(k-1,1))/h[k-1]
        +mu[k]*(A(k+1,1)-A(k,1))/h[k]);
    }
    lambda[n-1]=DE(1);
    mu[n-1]=DE(0);
    g(n-1,0) = 3 * (A(n-1,1) - A(n-2,1)) / h[n - 2];
    Matrix<DE> m=TridiagonalSolve({lambda, vector<DE>(n, DE(2)), mu}, g);
    return m;
}
template <class DE> DE Interpolation(const Matrix<DE> &A,DE x,string str)
{
    if(RowSize(A)<2)
    {
        cerr << "错误：至少需要两个点才能插值" << '\n';
        return x;
    }
    if(ColumnSize(A)<2)
    {
        cerr<<"错误：传入的数据矩阵至少两列"<<'\n';
        return x;
    }
    if(str=="Lagrange")
    {
        return LBasisToInterp(LagrangeBasis(A), A, x);
    }
    if(str=="Newton")
    {
        return DQTableToInterp(A, DifferenceQuotientTable(A),x);
    }
    if(str=="Hermite")
    {
        return DDTableToInterp(A, DividedDifferenceTable(A), x);
    }
    if(str=="non-standard Hermite")
    {
        return _3P3DHermiteInterp(A, x);
    }
    if(str=="linear" || str=="piecewise Lagrange")
    {
        int n = RowSize(A);
        int k= 1;
        while(k<n)
        {
            if(x<A(k,0))
            {
                return A(k-1,1) * (x - A(k,0)) / (A(k-1,0) - A(k,0)) + 
                       A(k,1) * (x - A(k-1,0)) / (A(k,0) - A(k-1,0));
            }
            ++k;
        }
        return A(n-2,1)*(x-A(n-2,0))/(A(n-2,0)-A(n-1,0))+
               A(n-1,1)*(x-A(n-2,0))/(A(n-1,0)-A(n-2,0));
    }
    if(str=="pchip" || str=="piecewise Hermite")
    {
        int n = RowSize(A);
        int i = 1;
        while(i<n)
        {
            if(x<A(i,0))
                return _2P3DHermiteInterp(A,i-1,i, x);
            ++i;
        }
        return _2P3DHermiteInterp(A,n-2,n-1, x);
    }
    if(str=="spline")
    {
        int n = RowSize(A);
        int k= 1;
        const Matrix<DE> &m = SplineSlope(A);
        Matrix<DE> B(2, 3);
        while(k<n)
        {
            if(x<A(k,0))
            {
                for (int i = 0; i < 2;++i)
                {
                    B(i,2) = m(k-1+i,0);
                    for (int j = 0; j < 2;++j)
                        B(i,j) = A(k-1+i,j);
                }
                return _2P3DHermiteInterp(B, 0, 1, x);
            }               
            ++k;
        }
        for (int i = 0; i < 2;++i)
        {
            B(i,2) = m(n-2+i,0);
            for (int j = 0; j < 2;++j)
                B(i,j) = A(n-2+i,j);
        }
        return _2P3DHermiteInterp(B, 0, 1, x);
    }
    cerr<<"错误：没有定义该方法"<<'\n';
    return x;
}
template <class DE> DE Interpolation(const Matrix<DE> &A,DE x)
{
    if(ColumnSize(A)<=2)
        return Interpolation(A,x,"linear");
    return Interpolation(A, x, "piecewise Hermite");
}
template <class DE> Matrix<DE> Interpolation(const Matrix<DE> &A,const Matrix<DE> &x,string str)
{
    if(RowSize(A)<2)
    {
        cerr << "错误：至少需要两个点才能插值" << '\n';
        return x;
    }
    if(ColumnSize(A)<2)
    {
        cerr<<"错误：传入的数据矩阵至少两列"<<'\n';
        return x;
    }
    if(str=="Lagrange")
    {
        vector<std::function<DE(DE)>> L = LagrangeBasis(A);
        Matrix<DE> s(RowSize(x), 1);
        for (int i = 0; i < RowSize(x);++i)
            s(i,0) = LBasisToInterp(L, A, x(i,0));
        return s;
    }
    if(str=="Newton")
    {
        Matrix<DE> D = DifferenceQuotientTable(A);
        Matrix<DE> s(RowSize(x), 1);
        for (int i = 0; i < RowSize(x);++i)
            s(i,0) = DQTableToInterp(A, D, x(i,0));
        return s;
    }
    if(str=="Hermite")
    {
        Matrix<DE> D = DividedDifferenceTable(A);
        Matrix<DE> s(RowSize(x), 1);
        for (int i = 0; i < RowSize(x);++i)
            s(i,0) = DDTableToInterp(A, D, x(i,0));
        return s;
    }
    if(str=="non-standard Hermite")
    {
        Matrix<DE> s(RowSize(x),1);
        for (int i = 0; i < RowSize(x);++i)
            s(i,0) = _3P3DHermiteInterp(A, x(i,0));
        return s;
    }
    if(str=="linear" || str=="piecewise Lagrange")
    {
        int n = RowSize(A);
        Matrix<DE> s(RowSize(x), 1);
        int k= 1;
        int j = 0;
        while(j<RowSize(x) && k<n)
        {
            if(x(j,0)<A(k,0))
            {
                s(j,0)=A(k-1,1) * (x(j,0)- A(k,0)) / (A(k-1,0) - A(k,0)) + 
                       A(k,1) * (x(j,0) - A(k-1,0)) / (A(k,0) - A(k-1,0));
                ++j;
                continue;
            }
            ++k;
        }
        for (int t = j; t < RowSize(x);++t)
        {
            s(t,0) =A(n-2,1)*(x(t,0)-A(n-2,0))/(A(n-2,0)-A(n-1,0))+
               A(n-1,1)*(x(t,0)-A(n-2,0))/(A(n-1,0)-A(n-2,0));
        }
        return s;
    }
    if(str=="pchip" || str=="piecewise Hermite")
    {
        int n = RowSize(A);
        Matrix<DE> s(RowSize(x), 1);
        int i=1;
        int j = 0;
        while(j<RowSize(x) && i<n)
        {
            if(x(j,0)<A(i,0))
            {
                s(j,0) = _2P3DHermiteInterp(A,i-1,i, x(j,0));
                ++j;
                continue;
            }
            ++i;
        }
        for (int t = j; t < RowSize(x);++t)
            s(t,0) = _2P3DHermiteInterp(A,n-2,n-1, x(t,0));
        return s;
    }
    if(str=="spline")
    {
        int n = RowSize(A);
        Matrix<DE> s(RowSize(x), 1);
        int k= 1;
        int j = 0;
        const Matrix<DE> &m = SplineSlope(A);
        Matrix<DE> B(2, 3);
        while(j<RowSize(x) && k<n)
        {
            if(x(j,0)<A(k,0))
            {
                for (int i = 0; i < 2;++i)
                {
                    B(i,2) = m(k-1+i,0);
                    for (int j = 0; j < 2;++j)
                        B(i,j) = A(k-1+i,j);
                }
                s(j,0) = _2P3DHermiteInterp(B, 0, 1, x(j,0));
                ++j;
                continue;
            }
            ++k;
        }
        for (int t = j; t < RowSize(x);++t)
        {
            for (int i = 0; i < 2;++i)
            {
                B(i,2) = m(n-2+i,0);
                for (int j = 0; j < 2;++j)
                    B(i,j) = A(n-2+i,j);
            }
            s(t,0)=_2P3DHermiteInterp(B, 0, 1, x(t,0));
        }
        return s;
    }
    cerr<<"错误：没有定义该方法"<<'\n';
    return x;
}
template <class DE> Matrix<DE> Interpolation(const Matrix<DE> &A, const Matrix<DE> &x)
{
    if(ColumnSize(A)<=2)
        return Interpolation(A,x,"linear");
    return Interpolation(A, x, "piecewise Hermite");
}
template <class DE> std::function<DE(DE)> Interpolation(const Matrix<DE> &A,string str)
{
    return [=](DE x) { return Interpolation(A, x, str); };
}

//积分
template<class DI> DI Integrate(std::function<DI(DI)> f,DI a,DI b,string str,int n)
{
    if(str=="Newton-Cotes"|| str=="N-C")
    {
        if(n<=0)
        {
            cerr<<"错误：Newton-Cotes法要求最后一个参数大于0"<<'\n';
            return (b - a) * f(a);
        }
        if(n>6)
        {
            cerr << "警告：插入过多的点可能会造成Newton-Cotes法数值不稳定。为了确保稳定性，已减少插入的点。" << '\n';
            return Integrate(f, a, b, "Newton-Cotes", 6);
        }
        const vector<vector<DI>> CotesTable({
            {0.5,0.5},
            {1.0/6,4.0/6,1.0/6},
            {0.125,0.375,0.375,0.125},
            {7.0/90,32.0/90,12.0/90,32.0/90,7.0/90},
            {19.0/288,75.0/288,50.0/288,50.0/288,75.0/288,19.0/288},
            {41.0/840,216.0/840,27.0/840,272.0/840,27.0/840,216.0/840,41.0/840}
        }); 
        DI s=0;
        DI divpoint = a;
        DI interval = b - a;
        DI h = interval / n;
        for (int k = 0; k <= n;++k)
        {
            s=s+CotesTable[n-1][k]*f(divpoint);
            divpoint = divpoint + h;
        }
        return interval * s;
    }
    if(str=="Gauss-Legendre" || str=="G-L" || str=="Gauss")
    {
        if(n<0)
        {
            cerr << "错误：Gauss-Legendre法要求至少插入一个点" << '\n';
            return (b - a) * f((a + b) / 2);
        }
        if(n>4)
        {
            cerr<<"警告：不支持插入过多点的Gauss-Legendre法。已自动减少插入的点。"<<'\n';
            return Integrate(f, a, b, "Gauss", 4);
        }
        const vector<vector<DI>> GaussPoint({
            {0},
            {0.5773502692,-0.5773502692},
            {0.7745966692,-0.7745966692,0},
            {0.8611363116,-0.8611363116,0.3399810436,-0.3399810436},
            {0.9061798459,-0.9061798459,0.5384693101,-0.5384693101,0}
        });
        const vector<vector<DI>> GaussCoef({
            {2},
            {1,1},
            {5.0/9,5.0/9,8.0/9},
            {0.3478548451,0.3478548451,0.6521451549,0.6521451549},
            {0.2369268851,0.2369268851,0.4786286705,0.4786286705,0.5688888889}
        });
        DI interval=b-a;
        DI halfinterval=interval/2;
        DI midpoint = (a + b) / 2;
        DI s = 0;
        for (int k = 0; k <=n;++k)
        {
            s = s + GaussCoef[n][k] * f(halfinterval*GaussPoint[n][k]+midpoint);
        }
        return halfinterval * s;
    }
    cerr << "错误：没有定义这种方法" << '\n';
    return f(a) * (b - a);
}
template<class DI> DI CompoundIntegrate(std::function<DI(DI)> f,DI a,DI b,string str,int n,int m)
{
    if(m<0)
    {
        cerr<<"错误：需要划分至少一个区间"<<'\n';
        return f(a) * (b - a);
    }
    DI s = 0;
    DI interval = b - a;
    DI h = interval / m;
    DI leftdiv = a;
    DI rightdiv = a + h;
    for (int i = 0; i < m;++i)
    {
        s = s + Integrate(f, leftdiv, rightdiv, str, n);
        leftdiv = rightdiv;
        rightdiv = rightdiv + h;
    }
    return s;
}
template <class DI> DI Integrate(std::function<DI(DI)> f,DI a,DI b,string str)
{
    if(str=="Romberg")
    {
        const DI epsilon = 1e-14;
        vector<DI> pre(1);
        vector<DI> now(1, (b - a) / 2 * (f(a) + f(b)));
        int k= 0;
        int times_max=20;
        DI accu = 0;
        int powof2 = 1;
        int powof4 = 1;
        do{
            ++k;
            pre = now;
            accu = 0;
            powof2 =(powof2<<1);
            for (int j = 0; j <(powof2>>1);++j)
            {
                accu = accu + f(a + (2 * j + 1) * (b - a) / powof2);
            }
            now[0] = 0.5 * pre[0] + (b - a) / powof2 * accu;
            powof4 = 4;
            for (int m= 1; m< k;++m)
            {
                now[m] = 1 / DI(powof4 - 1) * (powof4 * now[m - 1] - pre[m - 1]);
                powof4 = powof4 * 4;
            }
            now.emplace_back(1 / DI(powof4 - 1)*(powof4 * now[k-1] - pre[k-1]));
        }while(k<times_max && Abs(now[k]-pre[k-1])>=epsilon);
        if(k>=times_max)
            cerr << "警告：未在指定次数内达到预期精度"<<'\n';
        return now[k];
    }
    if(str=="trapz")
    {
        return Integrate(f, a, b, "Newton-Cotes", 1);
    }
    if(str=="compound trapz")
    {
        return CompoundIntegrate(f, a, b, "Newton-Cotes", 1, 100);
    }
    if(str=="Simpson")
    {
        return Integrate(f, a, b, "Newton-Cotes", 2);
    }
    if(str=="compound Simpson")
    {
        return CompoundIntegrate(f, a, b, "Newton-Cotes", 2, 100);
    }
    if(str=="Simpson3/8")
    {
        return Integrate(f, a, b, "Newton-Cotes", 3);
    }
    if(str=="compound Simpson3/8")
    {
        return CompoundIntegrate(f, a, b, "Newton-Cotes", 3, 100);
    }
    if(str=="Cotes"||str=="Milne")
    {
        return Integrate(f, a, b, "Newton-Cotes", 4);
    }
    if(str=="compound Cotes" || str=="compound Milne")
    {
        return CompoundIntegrate(f, a, b, "Newton-Cotes", 4, 100);
    }
    cerr << "错误：没有定义这种方法" << '\n';
    return f(a) * (b - a);
}
template <class DI> DI Integrate(std::function<DI(DI)> f,DI a,DI b)
{
    return Integrate(f, a,b,"Romberg");
}
template<class DI> DI CompoundIntegrate(std::function<DI(DI)>f ,DI a,DI b,string str,int m)
{
    if(m<0)
    {
        cerr<<"错误：需要划分至少一个区间"<<'\n';
        return f(a) * (b - a);
    }
    DI s = 0;
    DI interval = b - a;
    DI h = interval / m;
    DI leftdiv = a;
    DI rightdiv = a + h;
    for (int i = 0; i < m;++i)
    {
        s = s + Integrate(f, leftdiv, rightdiv, str);
        leftdiv = rightdiv;
        rightdiv = rightdiv + h;
    }
    return s;
}
template <class FS> FS FourierSinCoefficient(const std::function<FS(FS)> &f, int32_t n, FS w)
{
    if (Abs(w) < Machine_Epsilon)
    {
        std::cerr << "错误：频率不能为0" << std::endl;
        return 0;
    }
    return 2 * w / Pi<FS> * Integrate<FS>([&](FS t)
                           { return f(t) * sin(w * n * t); }, 0, Pi<FS> / w);
}
template <class FS> FS FourierCosCoefficient(const std::function<FS(FS)> &f, int32_t n, FS w)
{
    if (Abs(w) < Machine_Epsilon)
    {
        std::cerr << "错误：频率不能为0" << std::endl;
        return 0;
    }
    return 2 * w / Pi<FS> * Integrate<FS>([&](FS t)
                           { return f(t) * cos(w * n * t); }, 0, Pi<FS> / w);
}

//微分
template <class DD> DD D(std::function<DD(DD)> f,DD x,string str)
{
    if(str=="center")
    {
        const DD h = 1e-7;
        return (f(x + h) - f(x - h)) / (2 * h);
    }
    if(str=="forward")
    {
        const DD h = 1e-7;
        return (f(x + h) - f(x))/h;
    }
    if(str=="backward")
    {
        const DD h = 1e-7;
        return (f(x) - f(x-h))/h;
    }
    cerr << "错误：没有定义该方法。" << '\n';
    return x;
}
template <class DD> DD D(std::function<DD(DD)> f, DD x)
{
    return D(f, x, "center");
}
template<class DD> std::function<DD(DD)> D(std::function<DD(DD)> f)
{
    return [=](DD x) { return D(f, x); };
}
template<class DD> std::function<DD(DD)> D(int n,std::function<DD(DD)> f)
{
    if(n==0)
        return f;
    return D(n - 1, D(f));
}
template<class DD> DD D(int n,std::function<DD(DD)>f, DD x)
{
    return D(n, f)(x);
}
template<class DD> DD D(std::function<DD(Matrix<DD>)>f,int n,const Matrix<DD> &x)
{
    if(n>=RowSize(x))
    {
        cerr << "错误：求偏导的分量越界" << '\n';
        return f(x);
    }
    const DD h=1e-7;
    auto y=x;
    y(n,0)=Get(x,n,0)+h;
    return (f(y) - f(x)) / h;
}
template<class DD> std::function<DD(Matrix<DD>)> D(std::function<DD(Matrix<DD>)>f,int n)
{
    return [=](Matrix<DD> x){return D(f,n,x);};
}
template<class DD> Matrix<DD> D(const vector<std::function<DD(Matrix<DD>)>> &f,const Matrix<DD> &x)
{
    int n=f.size();
    Matrix<DD> J(n,n);
    for(int i=0;i<n;++i)
        for(int j=0;j<n;++j)
            J(i, j) = D(f[i], j, x);
    return J;
}

//方程(组)求根
template<class DD,class DM,class DF> 
Matrix<DD> Iteration(const DF &ini,const DM &phi,const Matrix<DD>&x0)
{
    const int times = 20;//迭代次数
    const DD eps = 1e-10;//精度
    int count = 0;
    auto x = x0;
    int n = RowSize(x0);
    Matrix<DD> fx(n, 1);
    Matrix<DD> B;
    ini(x,fx, B);
    do
    {
        phi(x, fx,B);
        ++count;
    } while (VectorNorm(fx,2) >= eps && count<times);
    if(count>=times)
        cerr << "警告：未在指定次数内达到预期精度" << '\n';
    return x;
}
template<class DD> 
Matrix<DD> FindRoot(const vector<std::function<DD(Matrix<DD>)>> &f,const Matrix<DD> &x0,string str)
{
    if(str=="Newton")
    {
        int n = RowSize(x0);
        auto ini=[](Matrix<DD> &,Matrix<DD>&,Matrix<DD>&){};
        auto phi = [=](Matrix<DD> &x, Matrix<DD> &fx, Matrix<DD> &B) 
        {
            for(int j=0;j<n;++j)
                fx(j, 0) = f[j](x);
            B = Inverse(D(f, x));
            x=x-B*fx;
        };
        return Iteration(ini, phi, x0);
    }
    if(str=="Broyden")
    {
        int n = RowSize(x0);
        auto ini=[=](Matrix<DD> &x,Matrix<DD> &fx,Matrix<DD> &B)
        {
            for(int j=0;j<n;++j)
                fx(j, 0) = f[j](x);
            B = Inverse(D(f, x));
        };
        auto phi = [=](Matrix<DD> &x, Matrix<DD> &fx, Matrix<DD> &B) 
        {
            auto fxpre = fx;
            auto xpre = x;
            x = x - B * fx;
            auto s = x - xpre;
            for(int j=0;j<n;++j)
                fx(j, 0) = f[j](x);
            auto y = fx - fxpre;
            auto stb = Transpose(s) * B;
            B = B + (s - B * y) * stb / (stb * y);
        };
        return Iteration(ini, phi, x0);
    }
    cerr << "错误：没有定义该方法" << '\n';
    return x0;
}
template<class DD> 
Matrix<DD> FindRoot(const vector<std::function<DD(Matrix<DD>)>> &f,const Matrix<DD> &x0)
{
    return FindRoot(f, x0, "Newton");
}
template<class DD> DD FindRoot(std::function<DD(DD)> f, DD a,DD b,string str)
{
    if(str=="bisection")
    {
        if(a>=b)
        {
            cerr << "错误：区间左端值需小于右端。" << '\n';
            return a;
        }
        if(f(a)*f(b)>0)
        {
            cerr << "错误：区间两端函数值同号，无法使用二分法。" << '\n';
            return a;
        }
        int times = 30;//迭代次数
        DD _left = a;
        DD _right=b;
        int n = 1;
        DD x = (_left + _right) / 2;
        while(n<=times)
        {
            if(f(x)==0) return x;
            if(f(_left)*f(x)<0)
            {
                _right=x;
            }
            else
            {
                _left = x;
            }
            x = (_left + _right) / 2;
            ++n;
        }
        return x;
    }
    if(str=="secant")
    {
        DD epsilon = 1e-14;
        if(Abs(a-b)<=epsilon)
        {
            cerr << "错误：弦截法的前两个初值选取的过于接近。" << '\n';
            return a;
        }
        int times = 100;//迭代次数
        int k = 0;
        DD pre = a;
        DD now = b;
        DD x =0;
        while (Abs(now-pre)>epsilon && k < times)
        {
            x = now - f(now) / (f(now) - f(pre)) * (now - pre);
            pre = now;
            now = x;
            ++k;
        } 
        if(k>=times)
        {
            cerr<<"警告：未在指定次数内达到预期精度。"<<'\n';
            cerr << "可能是初值的选择使迭代法不收敛，或者迭代函数的选择使收敛速度很慢。" << '\n';
        }
        return x;
    }
    cerr << "错误：没有定义该方法。" << '\n';
    return a;
}
template <class DD> DD FindRoot(std::function<DD(DD)> f, DD a,string str)
{
    if(str=="Newton")
    {
        std::function<DD(DD)> phi = [&](DD x) { return x-f(x) / D(f, x); };
        return Iteration(phi, a);
    }
    if(str=="multiplicity")
    {
        std::function<DD(DD)> F = [&](DD x) { return f(x) / D(f, x); };
        std::function<DD(DD)> phi = [&](DD x) { return x - F(x) / D(F, x); };
        return Iteration(phi, a);
    }
    if(str=="simplified Newton")
    {
        DD g = D(f, a);
        std::function<DD(DD)> phi = [&](DD x) { return x - f(x) / g; };
        return Iteration(phi, a);
    }
    if(str=="downhill")
    {
        DD lambda = 1;//下山因子
        DD b = 0;
        int k = 0;
        int times = 30;//下山次数
        for (k = 0; k <times;++k)
        {
            b = a - lambda * f(a) / D(f, a);
            if(Abs(f(b))<Abs(f(a)))
                break;
            lambda = lambda / 2;
        }
        if(k==times)
        {
            cerr << "警告：超出下山次数，未找到合适的下山因子。" << '\n';
        } 
        std::function<DD(DD)> phi=[&](DD x) { return x - lambda * f(x) / D(f, x); };
        return Iteration(phi, a);
    }
    cerr << "错误：没有定义该方法" << '\n';
    return a;
}
template <class DD> DD FindRoot(std::function<DD(DD)> f, DD a)
{
    return FindRoot(f, a, "Newton");
}
template <class DD> DD FindRoot(std::function<DD(DD)> f, DD a,string str,int m)
{
    if(str=="multiplicity")
    {
        std::function<DD(DD)> phi=[&](DD x) { return x - m* f(x) / D(f, x); };
        return Iteration(phi, a);
    }
    cerr << "错误：没有定义该方法" << '\n';
    return a;
}

//拟合与逼近
template<class DI,class DF,class DG> DI InnerProduct(DF f,DG g,DI a,DI b)
{
    std::function<DI(DI)> h = [&](DI x) { return f(x) * g(x); };
    return Integrate(h, a, b);
}
template<class DF> Matrix<DF> FindFit(const Matrix<DF> &data,const vector<std::function<DF(DF)>> &phi)
{
    int m=RowSize(data);
    int n = phi.size();
    Matrix<DF> f(1,m);
    vector<Matrix<DF>> Phi(n,Matrix<DF>(1,m));
    Matrix<DF> Gram(n,n);
    Matrix<DF> B(n, 1);
    for(int j=0;j<m;++j)
    {
        for (int i = 0; i < n;++i)
            Phi[i](0,j)= phi[i](Get(data,j,0));   
        f(0, j) = Get(data, j, 1);  
    }    
    for (int i = 0; i < n;++i)
    {
        for (int j =i; j < n;++j)
            Gram(i, j) = InnerProductR(Phi[i], Phi[j]);
        B(i, 0) = InnerProductR(Phi[i], f);
    }
    for (int i = 0; i < n;++i)
        for (int j = 0; j < i;++j)
            Gram(i, j) = Gram(j, i);
    return Transpose(LinearSolve(Gram,B));
}
template<class DF> Matrix<DF> FindFit(std::function<DF(DF)> f,const vector<std::function<DF(DF)>> &phi,DF a,DF b)
{
    int n = phi.size();
    Matrix<DF> Gram(n,n);
    Matrix<DF> B(n, 1);
    for(int i=0;i<n;++i)
    {
        for (int j = i;j<n;++j)
            Gram(i, j) = InnerProduct(phi[i],phi[j],a,b);
        B(i, 0) = InnerProduct(phi[i], f, a, b);
    }
    for (int i = 0; i < n;++i)
        for (int j = 0; j < i;++j)
            Gram(i, j) = Gram(j, i);
    return Transpose(LinearSolve(Gram,B));
}
template<class DF> DF Fit(const Matrix<DF> &data,const vector<std::function<DF(DF)>> &phi,DF x)
{
    int n = phi.size();
    const Matrix<DF> &Coef=FindFit(data,phi);
    DF s = 0;
    for (int i = 0;i<n;++i)
        s = s + Get(Coef,0,i)* phi[i](x);
    return s;
}
template<class DF> DF Fit(std::function<DF(DF)> f,const vector<std::function<DF(DF)>> &phi,DF a,DF b,DF x)
{
    int n = phi.size();
    const Matrix<DF> &Coef=FindFit(f,phi,a,b);
    DF s = 0;
    for(int i=0;i<n;++i)
        s=s+Get(Coef,0,i)*phi[i](x);
    return s;
}
template<class DF> std::function<DF(DF)> Fit(const Matrix<DF> &data,const vector<std::function<DF(DF)>> &phi)
{
    return [=](double x) { return Fit(data, phi, x); };
} 
template<class DF> std::function<DF(DF)> Fit(std::function<DF(DF)> f,const vector<std::function<DF(DF)>> &phi,DF a,DF b)
{
    return [=](double x) { return Fit(f, phi, a, b, x); };
} 
template <class DB> Matrix<DB> PseudoInverse(const Matrix<DB> &A)
{
    Matrix<DB> T = Transpose(A);
    return Inverse(T * A) * T;
}
template <class DB> Matrix<DB> LeastSquares(const Matrix<DB> &A, const Matrix<DB> &y)
{
    Matrix<DB> T = Transpose(A);
    return LinearSolve(T * A, T * y);
}
template <class DB> Polynomial<DB> PolyFit(const Matrix<DB> &data, const vector<Polynomial<DB>> &phi)
{
    int m = data.row_num;
    int n = phi.size();
    Matrix<DB> A(m,n);
    Matrix<DB> y(m, 1);
    for (int i = 0;i<m;++i)
    {
        for (int j = 0;j<n;++j)
        {
            A(i,j) =Get(phi[j],data(i,0));
        }
        y(i,0) = data(i,1);
    }
    const Matrix<DB> &c = LeastSquares(A, y);
    vector<DB> r(n);
    for (int i = 0; i <n;++i)
        r[i]=c(i,0);
    return Polynomial<DB>(r);
}
template <class DB> Polynomial<DB> PolyFit(const Matrix<DB> &xy,int deg)
{
    if(xy.column_num!=2)
    {
        cerr << "错误：传入数据矩阵的列数必须为2" << '\n';
        return X;
    }
    int n = xy.row_num;
    Matrix<DB> A(n, deg + 1);
    Matrix<DB> b(n, 1);
    for (int i = 0; i < n;++i)
    {
        b(i,0) = xy(i,1);
        A(i,0) = 1;
        for(int j=1;j<=deg;++j)
        {
            A(i,j) = A(i,j-1) * xy(i,0);
        }
    }
    const Matrix<DB> &c=LeastSquares(A, b);
    vector<DB> r;
    for (int i = 0; i <= deg;++i)
        r.push_back(c(i,0));
    Polynomial<DB> f(r);
    return f;
}
template <class DB> Polynomial<DB> PolyFit(std::function<DB(DB)> f,const vector<Polynomial<DB>>&phi,DB a,DB b)
{
    int n = phi.size();
    Matrix<DB> Gram(n, n);
    Matrix<DB> B(n, 1);
    for(int i=0;i<n;++i)
    {
        for (int j = i;j<n;++j)
            Gram(i, j) = InnerProduct(phi[i],phi[j],a,b);
        B(i, 0) = InnerProduct(phi[i], f, a, b);
    }
    for (int i = 0; i < n;++i)
        for (int j = 0; j < i;++j)
            Gram(i, j) = Gram(j, i);
    const Matrix<DB> &c = LinearSolve(Gram, B);
    vector<DB> r(n);
    for (int i = 0; i <n;++i)
        r[i]=c(i,0);
    return Polynomial<DB>(r);
}

}
#endif