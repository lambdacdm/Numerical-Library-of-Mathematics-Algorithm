//version 0.7.5
#include <iostream>
#include <ctime>
#include <chrono>
#include "numbertheory.h"
#include "highprecision.h"
using namespace std;
using namespace nbt;
using namespace hpc;

template<class DC> DC testmodpower(DC a,DC n,DC mod)
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

int main()
{
    BigInt a("100000000000000000000000567");
    BigInt m("5000000009");
    BigInt n("10000000019");
    //BigInt x = Fibonacci(1000000);
    //BigInt y = Fibonacci(100000);
        
    auto start = std::chrono::system_clock::now();
    /*for (unsigned long long i = 0; i < 10000000;++i)
    {
        if(MillerRabin(i)!=SimplePrimeQ(i))
        {
            cout << "wrong at " << i << endl;
            break;
        }
    }*/
    std::cout << MillerRabin<BigInt>(a) << std::endl;
    auto end = std::chrono::system_clock::now();
    std::cout << "total time=" << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
    system("pause");
    return 0;
}
                                    
                    
                                    