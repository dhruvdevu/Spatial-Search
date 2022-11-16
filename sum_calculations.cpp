#include<iostream>
#include<thread>
#include<cmath>
#include<chrono>
#include<vector>
#include<omp.h>
#include<cstdlib>
#include<time.h>

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
        #pragma omp parallel for collapse(2)
        for(int i=0; i<r; i++) {
            for(int j=0; j<res.size(); j++) {
                vector<int> copy = res[j];
                copy.push_back(i);
                sol[i*(pow(r, d-1))+j] = copy;
            }
        }
    }
    return sol;
}


double term(vector<int> w, vector<int> k, int n) { 
    double num = 1;
    double denom = d;

    for (int i = 0; i < d; i++) { 
        num *= pow(cos((M_PI*k[i]*(2*w[i]-1))/(2*n)), 2);
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

double errorterm(vector<int> w, vector<int> k, int n) {
    double denom = d;
    double num = (1-k[0]%2)*(cos((2*M_PI*(k[0]-1)*w[0])/n))*(M_PI/n)*(sin(M_PI*k[0]/n));
    for (int i = 0; i < d; i++) { 
        denom -= cos(M_PI*k[i]/n);
    }
    return 0.5*num/denom;
}

double eval_periodic_sum(vector<vector<int>>& k_vals, int n) {
    double sum = 0;
    double start = omp_get_wtime();
    #pragma omp parallel for reduction(+:sum)
    for(int i=1; i<k_vals.size(); i++) {
        sum += oldterm(k_vals[i], n);
    }
    double end = omp_get_wtime();
    //cout<<"parallel sum = "<<sum<<", time taken ="<<end-start<<"\n";
    return pow((double) 1/n, d)*sum;
}

double eval_aperiodic_sum(vector<int> w, vector<vector<int>> &k_vals, int n) {
    
    double sum = 0;
    double start = omp_get_wtime();
    #pragma omp parallel for shared(w) reduction(+:sum)
    for(int i=1; i<k_vals.size(); i++) {
        sum += term(w, k_vals[i], n);
    }
    double end = omp_get_wtime();
    //cout<<"parallel sum = "<<sum<<", time taken ="<<end-start<<"\n";
    return pow((double) 2/n, d)*sum;
}

double eval_first_order_error_sum(vector<int> w, vector<vector<int>> &k_vals, int n) {
    
    double sum = 0;
    #pragma omp parallel for shared(w) reduction(+:sum)
    for(int i=1; i<k_vals.size(); i++) {
        sum += errorterm(w, k_vals[i], n);
    }
    //cout<<"parallel sum = "<<sum<<", time taken ="<<end-start<<"\n";
    return pow((double) 2/n, d)*sum;
}

void printvec(vector<int> v) {
    cout<<"[";
    for(auto i : v) {
        cout<<i<<",";
    }
    cout<<"]";
}

void eval_first_order_error() {
    //Eval first order error
    vector<double> errorvec;
    for(int n = 20; n < 50; n+=3) {
        vector<int> w(d, floor(n/2));
        vector<vector<int>> k_vals = ktuples(n, d);
        double sum = eval_first_order_error_sum(w, k_vals, n);
        cout<<"n = "<<n<<", error sum = "<<sum<<"\n";
        errorvec.push_back(sum);
    }
    //printvec(errorvec);
}

void tail_sums(int n) {
    vector<int> w(d, floor(n/2));
    vector<vector<int>> k_vals = ktuples(n, d);
    for(int i=0; i < 5; i++) {
        w[i] = n;
        double sum = eval_aperiodic_sum(w, k_vals, n);
        printvec(w);
        cout<<"\naperiodic sum = "<<sum<<"\n";
    }
}

void random_points(int n, int num_points) {
        //Eval at random points
    srand(time(0));
    //int num_points = 100;
    //int n= 31;

    vector<vector<int>> points;
    vector<double> sums;
    vector<vector<int>> k_vals = ktuples(n, d);

    
    cout<<"n = "<<n<<"\n";
    double psum = eval_periodic_sum(k_vals, n);
    cout<<"periodic sum = "<<psum<<"\n";
    
    for(int i=0; i< num_points; i++){ 
        vector<int> w(d);
        for(int j = 0; j < d; j++) {
            w[j] = 1 + rand() % n;
        }
        points.push_back(w);
        
        double sum = eval_aperiodic_sum(w, k_vals, n);
        //printvec(w);
        sums.push_back(sum);
        //cout<<"\naperiodic sum = "<<sum<<"\n";
    }
    cout<<"[";
    for(auto p: points) {
        printvec(p);
        cout<<",";
    }
    cout<<"]\n[";

    for(auto s: sums) {
        cout<<s<<",";
    }
    cout<<"]\n";
    
}

void eval_along_line(int n) {
    int center = floor((n+1)/2);
    vector<int> w(d, center);
    double start = omp_get_wtime();
    vector<vector<int>> k_vals = ktuples(n, d);
    double end = omp_get_wtime();
    cout<<"tuple time taken ="<<end-start<<"\n";
    
    start = omp_get_wtime();
    double sum = eval_periodic_sum(k_vals, n);
    end = omp_get_wtime();
    cout<<"n = "<<n<<", periodic sum = "<<sum<<", time taken = "<<end-start<<"\n";
    
    start = omp_get_wtime();
    sum = eval_aperiodic_sum(w, k_vals, n);
    end = omp_get_wtime();
    cout<<"n = "<<n<<", aperiodic sum = "<<sum<<", time taken = "<<end-start<<"\n";
    for (int i=center-n/4; i<center+n/4; i++) {
        vector<int> w(d, i);
        sum = eval_aperiodic_sum(w, k_vals, n);
        cout<<"n = "<<n<<", w_0 = "<<i<<", aperiodic sum = "<<sum<<"\n";
    }
}

void eval_over_n(int start, int stop, int inc) {
    // Eval at center from n=start to stop
    for(int n = start; n < stop; n+=inc) {
        vector<int> w(d, floor((n+1)/2));
        vector<vector<int>> k_vals = ktuples(n, d);
        double psum = eval_periodic_sum(k_vals, n);
        double sum = eval_aperiodic_sum(w, k_vals, n);
        cout<<"n = "<<n<<", aperiodic sum = "<<sum<<", periodic sum = "<<psum<<"\n";
    }
}

int main () {
    random_points(32, 100);
    //eval_over_n(20, 50, 3);
    //tail_sums(31);
    //eval_along_line(31);
}