#include <iostream>
using namespace std;

#include <string>
#include "Long.h"

int main() {
    const Long a("12345");
    Long b("12345");
    cout << (a / 70);
    //cout << a << "\n" << b << "\n";
    //cout << a - b << "\n";
    a.toom3(b);
    return 0;
}