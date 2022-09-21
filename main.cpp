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
    BigInt a = ("10000000019");
    auto start = std::chrono::system_clock::now();
    cout<< MillerRabin<BigInt>(a) <<endl;
    auto end = std::chrono::system_clock::now();
    std::cout << "total time=" << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
    system("pause");
    return 0;
}
                                    
                    
                                    