//version 0.7.5
#include <iostream>
#include <ctime>
#include <chrono>
#include "numbertheory.h"
using namespace std;
using namespace nbt;
int main()
{
    for (unsigned long long i = 0; i < 1000000;++i)
    {
        if (SimplePrimeQ(i)!=MillerRabin(i))
        {
            cout << "wrong at " << i << endl;
            break;
        }   
    }
    system("pause");
    return 0;
}
                                    
                    
                                    