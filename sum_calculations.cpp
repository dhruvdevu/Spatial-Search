#include<iostream>
#include<thread>
#include<cmath>
#include<chrono>
#include<vector>
#include<omp.h>
#include<cstdlib>
#include<time.h>
#include<iomanip>  


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

vector<vector<int>>  ktuples_short(int r, int d) {
    vector<vector<int>> sol(pow(r, d-1)*(r-1.0)/2.0);
    vector<vector<int>> res = ktuples(r, d-1);
    //#pragma omp parallel for collapse(2)
    for(int i=0; i<(r-1.0)/2.0; i++) {
        for(int j=0; j<res.size(); j++) {
            vector<int> copy = res[j];
            copy.push_back(i);
            sol[i*(pow(r, d-1))+j] = copy;
        }
    }
    
    return sol;
}


double term(vector<int> w, vector<int> k, int n) { 
    double num = 1;
    double denom = d;

    for (int i = 0; i < d; i++) { 
        num *= pow(cos((M_PI*k[i]*(2*w[i]-1))/(2*n)), 2);
        denom -= cos(2*M_PI*k[i]/n);
    }
    return 0.5*num/denom;
}

double term_single_axis(vector<int> w, vector<int> k, int n) { 
    double num = pow(cos((M_PI*k[0]*(2*w[0]+n))/(2*n)), 2);
    //double num = pow(cos((M_PI*k[0]*(2*w[0]-1))/(2*n)), 2);
    /*
    if (k[0] % 2 != 0) {
        return 0;
    }
    */
    //double num = 1.0;
    double denom = d;
    for (int i = 0; i < d; i++) { 
        denom -= cos(2*M_PI*k[i]/n);
    }
    return 0.5*num/denom;
}

double oldterm(vector<int> k, int n) {
    double denom = d;
    for (int i = 0; i < d; i++) { 
        denom -= cos(2*M_PI*k[i]/n);
    }
    return 0.5/denom;
}

double oldterm_short(vector<int> k, int n) {
    double denom = d - cos(4*M_PI*k[0]/(n+1.0));
    for (int i = 1; i < d; i++) { 
        denom -= cos(2*M_PI*k[i]/n);
    }
    return 0.5/denom;
}

double error1helperterm(vector<int> k, int n) {
    double denom = d;
    for (int i = 1; i < d; i++) { 
        denom -= cos(2*M_PI*k[i]/n);
    }
    return 0.5/denom;
}


double error1term(vector<int> w, vector<int> k, int n) {
    double denom = d - cos(2*M_PI*k[0]/n);
    //double num = (1-k[0]%2)*(cos((2*M_PI*(k[0]-1)*w[0])/n))*(M_PI/n)*(sin(M_PI*k[0]/n));
    if (k[0] % 2 == 1) {
        return 0;
    }
    double num = 2*sin(M_PI*w[0]/n)*sin(M_PI*w[0]*(2*k[0]+1)/n);
    for (int i = 1; i < d; i++) { 
        denom -= cos(2*M_PI*k[i]/n);
    }
    if (k[0] == n - 1) {
        //return error1helperterm(k, n) + 0.5*num/denom;
    }
    
    return 0.5*num/denom;
}

// double error2term() {
//     if (k[0] % 2 == 1) {
//         return 0;
//     }
//     double num = cos(2*M_PI*w[0]*(k[0]+1)/n);
//     double denom_homo = d - cos(M_PI*(k[0]+1)/n;
//     double denom_inhomo = d - cos(2*M_PI*k[0]/n);
//     return 0.0;

// }
double error2helperterm(vector<int> w, vector<int> k, int n) {
    double denom = d-cos(2*M_PI*(n-1)/n);
    for (int i = 1; i < d; i++) { 
        denom -= cos(2*M_PI*k[i]/n);
    }
    return 0.5/denom;
}

double error2homoterm(vector<int> w, vector<int> k, int n) {
    double denom = d - cos(2*M_PI*k[0]/n);
    //double num = (1-k[0]%2)*(cos((2*M_PI*(k[0]-1)*w[0])/n))*(M_PI/n)*(sin(M_PI*k[0]/n));
    if (k[0] == n - 1) {
        return error2helperterm(w, k, n);
    }
    if (k[0] % 2 == 0) {
        return 0;
    }
    double num = cos(2*M_PI*w[0]*k[0]/n);
    for (int i = 1; i < d; i++) { 
        denom -= cos(2*M_PI*k[i]/n);
    }
    return -0.5*num/denom;
}

double error2inhomoterm(vector<int> w, vector<int> k, int n) {
    double denom = d;
    //double num = (1-k[0]%2)*(cos((2*M_PI*(k[0]-1)*w[0])/n))*(M_PI/n)*(sin(M_PI*k[0]/n));
    if (k[0] % 2 == 0) {
        return 0;
    }
    double num = cos(2*M_PI*w[0]*k[0]/n);
    for (int i = 1; i < d; i++) { 
        denom -= cos(2*M_PI*k[i]/n);
    }
    denom -= cos(2*M_PI*(k[0]-1)/n);
    if (denom == 0) {
        return 0;
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

double eval_aperiodic_single_axis_sum(vector<int> w, vector<vector<int>> &k_vals, int n) {
    
    double sum = 0;
    double start = omp_get_wtime();
    //#pragma omp parallel for shared(w) reduction(+:sum)
    for(int i=1; i<k_vals.size(); i++) {
        sum += term_single_axis(w, k_vals[i], n);
    }
    double end = omp_get_wtime();
    //cout<<"parallel sum = "<<sum<<", time taken ="<<end-start<<"\n";
    return 2*pow((double) 1/n, d)*sum;
}

double eval_error1_sum(vector<int> w, vector<vector<int>> &k_vals, int n) {
    
    double sum = 0;
    #pragma omp parallel for shared(w) reduction(+:sum)
    for(int i=1; i<k_vals.size(); i++) {
        sum += error1term(w, k_vals[i], n);
    }
    //cout<<"parallel sum = "<<sum<<", time taken ="<<end-start<<"\n";
    return pow((double) 1/n, d)*sum;
}

double eval_error2_sum(vector<int> w, vector<vector<int>> &k_vals, int n) {
    
    double homo_sum = 0;
    double inhomo_sum = 0;
    #pragma omp parallel for shared(w) reduction(+:homo_sum)
    for(int i=1; i<k_vals.size(); i++) {
        homo_sum += error2homoterm(w, k_vals[i], n);
    }
    #pragma omp parallel for shared(w) reduction(+:inhomo_sum)
    for(int i=1; i<k_vals.size(); i++) {
        inhomo_sum += error2inhomoterm(w, k_vals[i], n);
    }
    //cout<<"hsum"<<homo_sum<<" ihsum"<<inhomo_sum<<"\n";
    //cout<<"parallel sum = "<<sum<<", time taken ="<<end-start<<"\n";
    return pow((double) 1/n, d)*(homo_sum + inhomo_sum);
}

void printvec(vector<int> v) {
    cout<<"[";
    for(auto i : v) {
        cout<<i<<",";
    }
    cout<<"]";
}

void printvec(vector<double> v) {
    cout<<"[";
    for(auto i : v) {
        cout<<setprecision(10)<<i<<",";
    }
    cout<<"]";
}


// void eval_first_order_error(int t) {
//     //Eval first order error
//     vector<double> errorvec;
//     //cout<<"n="<<n<<"\n";
    
//     for(int n = 21; n <= 60; n+=2) {
//         vector<int> w(d, floor(n/2));
//         vector<vector<int>> k_vals = ktuples(n, d);
//         //w[0] = t;
//         w[0] = sqrt(n);
//         double sum = eval_first_order_error_sum(w, k_vals, n);
//         cout<<"n = "<<n<<", error sum = "<<sum<<"\n";
//         errorvec.push_back(sum);
//     }
//     printvec(errorvec);
// }

void eval_error1() {
    //Eval first order error
    vector<double> errorvec;
    //cout<<"n="<<n<<"\n";
    
    for(int n = 21; n <= 50; n+=2) {
        //vector<int> w(d, floor(n/2));
        vector<int> w(d, 0);
        vector<vector<int>> k_vals = ktuples(n, d);
        //w[0] = t;
        w[0] = 1;
        double sum = eval_error1_sum(w, k_vals, n);
        cout<<"n = "<<n<<", error1 sum = "<<setprecision(10)<<sum<<"\n";
        errorvec.push_back(sum);
    }
    printvec(errorvec);
}

void eval_error2() {
    //Eval first order error
    vector<double> errorvec;
    //cout<<"n="<<n<<"\n";
    
    for(int n = 21; n <= 50; n+=2) {
        //vector<int> w(d, floor(n/2));
        vector<int> w(d, 0);
        vector<vector<int>> k_vals = ktuples(n, d);
        //w[0] = t;
        w[0] = 1;
        double sum = eval_error2_sum(w, k_vals, n);
        cout<<"n = "<<n<<", error2 sum = "<<setprecision(10)<<sum<<"\n";
        errorvec.push_back(sum);
    }
    printvec(errorvec);
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

void eval_over_n(int start, int stop, int inc, int p) {
    // Eval at center from n=start to stop
    vector<double> apvals;
    bool pflag = false;
    for(int n = start; n < stop; n+=inc) {
        //vector<int> w(d, floor((n+1.0)/2.0));
        vector<int> w(d, 0);
        vector<vector<int>> k_vals = ktuples(n, d);
        //vector<vector<int>> k_vals = ktuples_short(n, d);
        //w[0] = p;
        //w[0] = ((double) n)/2.0;
        if (pflag) {
        double psum = 2*eval_periodic_sum(k_vals, n);
        cout<<"n = "<<n<<", periodic sum = "<<setprecision(10)<<psum<<"\n";
        apvals.push_back(psum);
        } else {
        double sum = eval_aperiodic_single_axis_sum(w, k_vals, n);
        //double sum = eval_aperiodic_sum(w, k_vals, n);
        cout<<"n = "<<n<<", aperiodic sum = "<<setprecision(10)<<sum<<"\n";//", periodic sum = "<<psum<<"\n";
        apvals.push_back(sum);
        }
    }
    printvec(apvals);
}

void test_sums() {
    int n = 21;
    vector<int> w(d, 0);
    vector<vector<int>> k_vals = ktuples(n, d);
    double true_val = eval_periodic_sum(k_vals, n);
    double s = eval_aperiodic_single_axis_sum(w, k_vals, n);
    double e1 = eval_error1_sum(w, k_vals, n);
    double e2 = eval_error2_sum(w, k_vals, n);
    cout<<s-true_val<<" "<<e1<<" "<<e2<<" "<<s-true_val-e2<<"\n";
    if (d == 1) { 
        double calc_periodic_sum = 0.0;
        for (int k = 1; k < n; k++) {
            calc_periodic_sum += (0.5/n)/(1.0-cos(2*M_PI*k/n));
        }
        
        double calc_aperiodic_sum = 0.0;
        for (int k = 1; k < n; k++) {
            calc_aperiodic_sum += (1.0/n)*pow(cos(M_PI*k/2),2)/(1.0-cos(2*M_PI*k/n));
        }
        //double calc_error = 0.5*(-1.0/(1.0-cos(2*M_PI/3)) + 1.0/(1.0-cos(4*M_PI/3)))/n;
        double calc_e1 = 0; //0.5/(d-cos(2*M_PI*(n-1)/n));
        //double calc_e2 = -0.5/(d-cos(2*M_PI*(1)/n)) + 0.5/(d-cos(2*M_PI*(n-1)/n));
        double calc_e2 = -(1.0/n)*1.0/(1.0-cos(2*M_PI/n))+(1.0/n)*1.0/(1.0-cos(2*M_PI*(n-1)/n));
        for (int k = 3; k < n; k+=2) {
            calc_e2 += (1.0/n)*(1.0/(1.0-cos(2*M_PI*(k-1)/n)) - 1.0/(1.0-cos(2*M_PI*(k)/n)));
        }
        cout<<true_val<<" "<<calc_periodic_sum<<" "<<s<<" "<<calc_aperiodic_sum<<"\n";
        cout<<true_val-s<<" "<<calc_e1<<" "<<calc_e2<<" "<<e1<<" "<<e2;
    }
    if (d == 2) {
        double calc_periodic_sum = 0.2222222222222222;
        double calc_aperiodic_sum =  0.2962962962962962;
        double calc_error = -0.07407407407407401;
        double calc_e1 = 0; //0.5/(d-cos(2*M_PI*(n-1)/n));
        double calc_e2 = 0;
        cout<<true_val<<" "<<calc_periodic_sum<<" "<<s<<" "<<calc_aperiodic_sum<<"\n";
    }
    //(cos(2*M_PI*k*w/n) - cos(2*M_PI*(k+1)*w/n))/(1.0-cos(2*M_PI*k/n));
    //double calc_e2 = 0.0; 
    
}

int main () {
    //random_points(32, 100);
    eval_over_n(21, 36, 2, 0);
    //eval_over_n(20, 32, 2, 0);
    //test_sums();
    //eval_over_n(21, 50, 2, 1);
    //Need to make this a torus on all but one axis
    //tail_sums(31);
    //eval_along_line(31);
    //eval_error1();
    //eval_error2();
    // vector<int> w(d, 0);
    // w[0] = 1;
    // int n = 63;
    // vector<vector<int>> k_vals = ktuples(n, d);
    // double s = eval_periodic_sum(k_vals, n);
    // cout<<"n = "<<n<<", periodic sum = "<<setprecision(14)<<s<<"\n";
    


    
    //cout<<s-true_val<<","<<e1-e2<<","<<e2;
}