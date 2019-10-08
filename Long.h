//
// Created by Slava Karpii on 9/13/19.
//

#ifndef LONGA_LONG_H
#define LONGA_LONG_H

using namespace std;

#include <string>
#include <iostream>
#include <complex>
#include <vector>
class Mult;

class Long {
    friend class Mult;
    friend class Karatsuba;
    friend class Toom3;
    friend class Sen_strass;
protected:
    static Mult *t;
    static int base;

    string digits = "0";
    void binary();
    char sign = '+';
    char negative_sign() const;
    int Size = 1;
    int digit(int i) const;
    void set(int i, int val);
    void trim();
    bool comp_abs(const Long& other) const;
    Long add_abs(const Long& other) const ;
    Long sub_abs(const Long& other) const;
    int bin_search(int a, int b, const Long& M, const Long& n) const;
    int fl = 0;
public:
    static void set(Mult *p){
        t = p;
    }

    Long(string s = "0");
    Long(int n);
    string bin = "";
    void back();
    static void bswitch(int b);
    bool operator==(const Long& other) const;
    bool operator>(const Long& other) const;
    bool operator>=(const Long& other) const;
    bool operator<(const Long& other) const;
    bool operator<=(const Long& other) const;
    bool operator!=(const Long& other) const;
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
istream& operator>> (istream& in, Long& a);
Long GCD(const Long& first, const Long& other);
Long Jacobi( Long n,  Long k);


class Mult{
public:
    Long naive(const Long& a, const Long& other) const;;

    virtual Long mult(const Long& a, const Long& other) const {
        return naive(a, other);
    };
};

class Karatsuba: public Mult {
    Long karatsuba(const Long& a, const Long& other) const;
    virtual Long mult(const Long& a, const Long& other) const {
        return karatsuba(a, other);
    };
};

class Toom3: public Mult {
    Long toom3(const Long& a, const Long& other) const;
    virtual Long mult(const Long& a, const Long& other) const {
        return toom3(a, other);
    };
};

class Sen_strass: public Mult {
    Long sen_strass(const Long& a, const Long& other) const;
    virtual Long mult(const Long& a, const Long& other) const {
        return sen_strass(a, other);
    };
private:
    typedef complex<double> b;
    void fft(vector<b> &a, bool invert) const;
};

ostream& operator<<(ostream &out, Long i);


#endif //LONGA_LONG_H
