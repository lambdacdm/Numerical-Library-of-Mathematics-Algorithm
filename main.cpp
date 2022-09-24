//version 0.7.5
#include <iostream>
#include <ctime>
#include <chrono>
#include "numbertheory.h"
#include "highprecision.h"
using namespace std;
using namespace nbt;
using namespace hpc;

int main()
{
    BigInt a("986732000000000000000123567");
    BigInt m("5000000009");
    BigInt n("10000000019");
    //BigInt x = Fibonacci(1000000);
    BigInt y = Fibonacci(100000);
        
    auto start = std::chrono::system_clock::now();
    std::cout << log2(a) << std::endl;
    auto end = std::chrono::system_clock::now();
    std::cout << "total time=" << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
    system("pause");
    return 0;
}
                                    
                    
                                    