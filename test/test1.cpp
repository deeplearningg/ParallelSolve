/******************************************************************************
*******************************************************************************/

#include <iostream>
#include <iostream>
#include <vector>
#include <ctime>
#include <assert.h>
#include <cstring>
using namespace std;

int main()
{
    double *x1,*x2,*x3;
    size_t N = 100000;
    cout << N << endl;
    x1 = new double[N];
    x2 = new double[N];
    x3 = new double[N];
    for(int i = 0;i < N;i++){
        x1[i] = i;
    }
    clock_t t1 = clock();
    for(int i = 0;i < N;i++)
        x2[i] = x1[i];
    clock_t t2 = clock();
    memcpy(x3,x1,sizeof(double)*N);
    clock_t t3 = clock();
    for(int i = 0;i < N;i++)
        assert(x2[i] == x3[i]);

    cout << " for loop takes " << (t2 - t1)/CLOCKS_PER_SEC << endl;
    cout << " memcpy takes " << (t3 - t2)/CLOCKS_PER_SEC << endl;

    return 0;
}

