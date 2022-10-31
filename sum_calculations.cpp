#include<iostream>
#include<thread>
#include<cmath>
#include<chrono>
#include<vector>
#include<omp.h>

using namespace std;
using namespace std::chrono;

const int d = 5;

//TODO: parallelize this
vector<vector<int>>  ktuples(int r, int d) {
    vector<vector<int>> sol(pow(r, d));
    if (d == 0) {
        vector<int> l;
        sol[0] =l;
    } else {
        vector<vector<int>> res = ktuples(r, d-1);
        //#pragma omp parallel for collapse(2)
        for(int i=0; i<r; i++) {
            for(int j=0; j<res.size(); j++) {
                vector<int> copy = res[j];
                copy.push_back(i);
                sol[i*r+j] = copy;
            }
        }
    }
    return sol;
}


double term(vector<int> w, vector<int> k, int n) { 
    double num = 1;
    double denom = d;

    for (int i = 0; i < d; i++) { 
        num *= pow(cos(M_PI*k[i]*(2*w[i]-1)/(2*n)), 2);
        denom -= cos(M_PI*k[i]/n);
    }
    return 0.5*num/denom;
}

double oldterm(vector<int> k, int n) {
    double denom = d;
    for (int i = 0; i < d; i++) { 
        denom -= cos(M_PI*k[i]/n);
    }
    return 0.5/denom;
}

double eval_periodic_sum(int n) {
    double start = omp_get_wtime();
    vector<vector<int>> k_vals = ktuples(n, d);
    double end = omp_get_wtime();
    cout<<"tuple time taken ="<<end-start<<"\n";
    double sum = 0;
    start = omp_get_wtime();
    for(int i=1; i<k_vals.size(); i++) {
        sum += oldterm(k_vals[i], n);
    }
    end = omp_get_wtime();
    cout<<"sequential sum = "<<sum<<", time taken ="<<end-start<<"\n";
    sum = 0;
    start = omp_get_wtime();
    #pragma omp parallel for reduction(+:sum)
    for(int i=1; i<k_vals.size(); i++) {
        sum += oldterm(k_vals[i], n);
    }
    end = omp_get_wtime();
    cout<<"parallel sum = "<<sum<<", time taken ="<<end-start<<"\n";
    return pow((double) 1/n, d)*sum;
}

double eval_aperiodic_sum(vector<int> w, int n) {
    vector<vector<int>> k_vals = ktuples(n, d);
    double sum = 0;
    //TODO: Optimize this by threading
    for(int i=0; i<k_vals.size(); i++) {
        sum += term(w, k_vals[i], n);
    }
    return pow(((double) 1/n), d)*sum;
}



int main () {
    int n = 35;
    double sum = eval_periodic_sum(n);
    cout<<"n = "<<n<<", periodic sum = "<<sum;
}