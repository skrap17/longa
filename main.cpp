#include <iostream>
using namespace std;

#include <string>
#include "Long.h"

int main() {
    Long a("12345");
    Long b("1234");
    cout << string(a) << "\n" << string(b) << "\n";
    cout << string(a.karatsuba(b)) << "\n";
    return 0;
}