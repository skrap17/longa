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
    Long add_abs(const Long& other) const ;
    Long sub_abs(const Long& other) const;
public:
    Long(string s = "0");
    Long(int n);
    bool operator==(const Long& other);
    bool operator>(const Long& other);
    bool operator>=(const Long& other);
    bool operator<(const Long& other);
    bool operator<=(const Long& other);
    bool operator!=(const Long& other);
    Long operator+(const Long& other) const;
    Long operator-(const Long& other) const;
    Long operator*(const Long& other) const;
    Long operator*(int n) const;
    Long operator/(int n) const;
    Long& operator<<(int n);
    Long operator<<(int n) const;

    Long karatsuba(const Long& other) const;
    Long toom3(const Long& other) const;
    explicit operator string();
};

ostream &operator<<(ostream &out, Long i);

ostream& operator<<(ostream &out, Long i);


#endif //LONGA_LONG_H
