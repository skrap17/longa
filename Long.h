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
#include <tuple>
class Mult;
class Long {
    friend Long abs(const Long& a);
    friend class Mult;
    friend class Karatsuba;
    friend class Toom3;
    friend class Sen_strass;
    friend class Modular;
protected:
    static Mult *t;
    static int base;
    string bin = "";
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
    static void bswitch(int b);
    Long naive(const Long &other) const;
public:
    static void set(Mult *p){
        t = p;
    }
    Long(string s = "0");
    Long(int n);
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
    Long pow(const Long& other) const;
    Long pow_mod(const Long& other, const Long& modulus) const;
    bool fermat() const;
    bool rab_mil() const;
    bool sol_strass() ;
    string cook(bool decimal = true) const;
    Long div_cook(const Long& other);
    operator string();
};




ostream &operator<<(ostream &out, Long i);
istream& operator>> (istream& in, Long& a);
Long GCD(const Long& first, const Long& other);
Long Jacobi( Long n,  Long k);
Long abs(const Long& a);

class Mult{
public:
    Long naive(const Long& a, const Long& other) const;

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

class Modular: public Mult {
    Long modular(const Long& a, const Long& other) const;
    virtual Long mult(const Long& a, const Long& other) const {
        return modular(a, other);
    };
private:
    int q(int k) const;
    vector<Long> m(int k) const;
    Long gcdex (Long a, Long b, Long & x, Long & y) const;
    Long inv_mod(Long a, Long m) const;
};



ostream& operator<<(ostream &out, Long i);


#endif //LONGA_LONG_H
