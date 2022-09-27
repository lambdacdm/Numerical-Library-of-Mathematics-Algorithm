#ifndef PLOT_H
#define PLOT_H
#include "linearalgebra.h"
using namespace lag;

namespace plt{

//-----声明部分-----

//图像绘制
template <class DG>
class Graphics
{
    public:
        Graphics();
        Graphics(const Matrix<bool> &, const std::pair<DG, DG> &, const std::pair<DG, DG> &);
        void show();

    private:
        Matrix<bool> pixel;
        std::pair<DG,DG> xrange;
        std::pair<DG,DG> yrange;
};

//-----定义部分-----

//图像绘制
template<class DG> Graphics<DG>::Graphics()
{
    pixel = Matrix<bool>(0, 0);
    xrange = std::pair<DG, DG>{0, 0};
    yrange = std::pair<DG, DG>{0, 0};
}
template<class DG> Graphics<DG>::Graphics
(const Matrix<bool> &d, const std::pair<DG, DG> &xr, const std::pair<DG, DG> &yr)
{ 
    pixel=d;
    xrange=xr;
    yrange = yr;
}
template<class DG> void Graphics<DG>::show()
{
    int n=RowSize(pixel);
    int m=ColumnSize(pixel);
    std::cout << '\n';
    for (int i = 0;i<n;++i)
    {
        for(int j=0;j<m;++j)
        {
            if(Get(pixel,i,j))
                std::cout<<'*';
            else
                std::cout << ' ';
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}
template <class DG> Graphics<DG> Plot
(const Matrix<DG> &data,std::pair<DG,DG> xr,std::pair<DG,DG> yr,unsigned xdots, unsigned ydots)
{
    DG xmin=std::get<0>(xr);
    DG xmax=std::get<1>(xr);
    DG ymin=std::get<0>(yr);
    DG ymax=std::get<1>(yr);
    Matrix<bool> pixelmatrix(ydots+1, xdots+1);
    DG xh,yh;
    if(xdots==0)
        xh = 1;
    else
        xh = (xmax - xmin) / xdots;
    if(ydots==0)
        yh = 1;
    else
        yh = (ymax - ymin) / ydots;
    const int datarows=RowSize(data);
    for (int i = 0; i < datarows;++i)
    {
        int dotcolumn=int((Get(data,i,0)-xmin)/xh);
        int dotrow = int((ymax-Get(data, i, 1)) / yh);
        Make(pixelmatrix,dotrow,dotcolumn,true);
    }
    Graphics<DG> gr(pixelmatrix, xr, yr);
    gr.show();
    return gr;
}
template <class DG> Graphics<DG> Plot(const Matrix<DG> &data, unsigned xdots, unsigned ydots)
{
    DG xmin = ColumnMin(data, 0);
    DG xmax = ColumnMax(data, 0);
    DG ymin = ColumnMin(data, 1);
    DG ymax = ColumnMax(data, 1);
    std::pair<DG,DG> xr({xmin,xmax});
    std::pair<DG,DG> yr({ymin, ymax});
    return Plot(data, xr, yr, xdots, ydots);
}
template <class DG> Graphics<DG> Plot(const Matrix<DG> &data, unsigned xdots, string str)
{
    if(str=="unitscale")
    {
        DG xmin = ColumnMin(data, 0);
        DG xmax = ColumnMax(data, 0);
        DG ymin = ColumnMin(data, 1);
        DG ymax = ColumnMax(data, 1);
        std::pair<DG,DG> xr({xmin,xmax});
        std::pair<DG,DG> yr({ymin, ymax});
        unsigned ydots=1;
        DG ratio=(ymax - ymin) / (xmax - xmin);
        const DG maxratio = DG(3);
        if(ratio>maxratio)
        {
            cerr<<"警告：y方向上函数值变动过大，已自动调整图像大小\n";
            ratio=maxratio;
        }
        ydots=unsigned(xdots*ratio);
        return Plot(data, xr, yr, xdots, ydots);
    }
    if(str=="square")
    {
        return Plot(data, xdots, xdots);
    }
    else
    {
        cerr<<"错误：未定义该方法\n";
        return Plot(data,xdots,"unitscale");
    }
}
template <class DG> Graphics<DG> Plot(const Matrix<DG> &data, unsigned xdots)
{
    return Plot(data, xdots, "unitscale");
}
template <class DG> Graphics<DG> Plot(const Matrix<DG> &data,string str)
{
    return Plot(data, RowSize(data), str);
}
template <class DG> Graphics<DG> Plot(const Matrix<DG> &data)
{
    return Plot(data, "unitscale");
}
template <class DG> Graphics<DG> Plot(const std::function<DG(DG)> &f,const vector<DG> &range)
{
    DG xmin=range[0];
    DG xmax=range[1];
    int xdots = 1;
    if(range.size()>=3)
        xdots=range[2];
    else
        xdots = 100;
    Matrix<DG> data(xdots+1, 2);
    DG xh = (xmax - xmin) / xdots;
    DG pre = xmin;
    for (int i = 0; i<=xdots;++i)
    {
        data(i, 0) = pre;
        data(i, 1) = f(Get(data,i,0));
        pre = pre + xh;
    }
    return Plot(data);
}

}
#endif