#include <iostream>
using namespace std;

#include "Long.h"

int main() {
    Long a("198660");
    Long b("101101");
    cout << string(a) << "\n" << string(b) << "\n";
    cout << string(a - b) << "\n";
    return 0;
}