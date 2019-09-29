#include <iostream>
using namespace std;

#include <string>
#include "Long.h"

int binpoisk(int a, int b, int M, int n){
    if (a == b){
        return a;
    }
    int c = (a + b) / 2 + 1;
    if ( n * c > M)
        return binpoisk(a, c - 1, M, n);
    else
        return binpoisk(c, b, M, n);
}

int main() {
    const Long a("12345");
    Long b("12345");
    cout << (a / 70);
    //cout << a << "\n" << b << "\n";
    //cout << a - b << "\n";
    a.toom3(b);
    return 0;
}