#include<iostream>
#include<thread>
#include<cmath>
#include<vector>
#include<chrono>
#include<omp.h>
// #include "/usr/local/opt/libomp/include/omp.h"


using namespace std;
using namespace std::chrono;

double square(int x) {
    return x*x;
}

int main () {
    vector<int> a = {};
    for(int i=0; i<160000000; i++) {
        a.push_back(i);
    }
    double sum = 0;
    double start = omp_get_wtime();
    for(int i=0; i<160000000; i++) {
        sum += a[i]*a[i];
    }
    double end = omp_get_wtime();
    cout<<"sequential sum = "<<sum<<", time taken ="<<end-start<<"\n";
    
    sum = 0;
    start = omp_get_wtime();
    #pragma omp parallel for reduction(+:sum)
    for(int i=0; i<160000000; i++) {
        sum += a[i]*a[i];
    }
   end=omp_get_wtime();
    cout<<"parallel sum = "<<sum<<", time taken ="<<end-start<<"\n";
}