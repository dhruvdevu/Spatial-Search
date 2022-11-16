#include<iostream>
#include<vector>

using namespace std;
int main() {
    vector< vector<int> > a(10);
    for(int i=0;i<10;i++) {
        a[i].push_back(i);
    }
    for(int i=0;i<10;i++) {
        cout<<a[i][0];
    }
}