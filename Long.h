//
// Created by Slava Karpii on 9/13/19.
//

#ifndef LONGA_LONG_H
#define LONGA_LONG_H

using namespace std;

#include <string>


class Long {
    static const int base = 10;
    string digits = "0";
    char sign = '+';
    char negative_sign() const;
    int Size = 1;
    int digit(int i) const;
    void set(int i, int val);
    void trim();
    bool comp_abs(const Long& other) const;
    Long add_abs(const Long& other);
    Long sub_abs(const Long& other);
public:
    Long(string s = "0");
    bool operator==(const Long& other);
    bool operator>(const Long& other);
    bool operator>=(const Long& other);
    bool operator<(const Long& other);
    bool operator<=(const Long& other);
    bool operator!=(const Long& other);
    Long operator+(const Long& other);
    Long operator-(const Long& other);
    Long operator*(const Long& other);
    Long operator<<(int n);
    Long karatsuba(const Long& other);
    explicit operator string();
};


#endif //LONGA_LONG_H
