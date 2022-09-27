#include <iostream>
#include <ctime>
#include <chrono>
#include "../lib/mathalgorithm.h"
using namespace mal;

/* minimal inclusion:
#include "../lib/highprecision.h"
#include "../lib/numbertheory.h"
using namespace hpc;
using namespace nbt; */

int main()
{
    auto start = std::chrono::system_clock::now();
    BigInt a("23948576235");
    BigInt b("823761298");
    BigInt c("2");

    BigInt m("986732000000000000000123761");
    BigInt n("5000000009");
    BigInt x = Fibonacci(1000000);
    BigInt y = Fibonacci(100000);

    std::cout << "a * b = " << a*b << std::endl;
    std::cout << "c ^ 100 = " << (c^100) << std::endl;
    std::cout << "Prime test for m: " << PrimeQ(m) << std::endl;
    std::cout << "Find a factor of n: " << FindAFactor(n) << std::endl;
    std::cout << "integer length of x is " << IntegerLength(x) << std::endl;
    std::cout << "x mod y is " << x % y << std::endl;
    std::cout << "log2 of a is " << log2(a) << std::endl;
    auto end = std::chrono::system_clock::now();
    std::cout << "total time=" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
    << "ms"<<std::endl;
    system("pause");
    return 0;
}                
                                    