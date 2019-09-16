//
// Created by Slava Karpii on 9/13/19.
//

#include <iostream>
#include "Long.h"

Long::Long(string s) {
    if (s.find_first_not_of("+-0123456789") == string::npos) {
        digits = s;
        if (s[0] == '-' || s[0] == '+') {
            sign = s[0];
            digits = digits.substr(1);
        } else
            sign = '+';
        trim();
    }
    else {
        digits = "0";
        Size = 1;
        sign = '+';
    }
}

bool Long::comp_abs(const Long& other) const {
    if (Size > other.Size)
        return true;
    else if (Size < other.Size)
        return false;
    else {
        return digits > other.digits;
    }
}

bool Long::operator==(const Long& other) {
    return digits == other.digits;
}

bool Long::operator>(const Long& other) {
    if (sign == '+'){
        if (other.sign == '+')
            return comp_abs(other);
        else
            return true;
    }
    else {
        if (other.sign == '-')
            return !comp_abs(other);
        else
            return false;
    }
}

bool Long::operator>=(const Long& other) {
    return *this > other || *this == other;
}

bool Long::operator<(const Long& other) {
    return !(*this >= other);
}

bool Long::operator<=(const Long& other) {
    return !(*this > other);
}

bool Long::operator!=(const Long& other) {
    return !(*this == other);
}

Long Long::add_abs(const Long &other) {
    Long res, tmp;
    if (comp_abs(other)) {
        res = *this;
        tmp = other;
    }
    else {
        res = other;
        tmp = *this;
    }
    int n = tmp.Size;
    int N = res.Size;
    int carry = 0, s = 0;
    int i = 1;
    for (; i <= N; i++){
        s = res.digit(N - i) + tmp.digit(n - i) + carry;
        carry = s / base;
        s %= base;
        //cout << s << ' ' << carry << '\n';
        res.set(N - i, s);
    }
    if (carry != 0){
        if (n == N)
            res.set(N - i, carry);
        else
            res.set(N - i, res.digit(N - i) + carry);
    }
    return res;
}

int Long::digit(int i) const {
    if (i >= 0 && i < Size)
        return int(digits[i]) - 48;
    else
        return 0;
}

void Long::set(int i, int val) {
    if (i >= 0)
        digits[i] = val + 48;
    else{
        digits = to_string(val) + digits;
        Size++;
    }
}

Long::operator string() {
    string s;
    if (sign == '-')
        s = "-";
    return s + digits;
}

Long Long::sub_abs(const Long &other) {
    Long res, tmp;
    if (comp_abs(other)) {
        res = *this;
        tmp = other;
    }
    else {
        res = other;
        tmp = *this;
    }
    int n = tmp.Size;
    int N = res.Size;
    int i;
    int s = 0, carry = 0;
    for (i = 1; i <= n; i++){
        s = res.digit(N - i) - tmp.digit(n - i) - carry;
        carry = 0;
        //cout << s << " " << carry << "   ";
        if (s <  0){
            s += base;
            carry = 1;
        }
        //cout << s << " " << carry << "\n";
        res.set(N - i, s);
    }
    if (carry != 0){
        res.set(N - i, res.digit(N - i) - 1);
    }
    res.trim();
    return res;
}

void Long::trim() {
    int n = digits.find_first_not_of('0');
    if (n != string::npos){
        digits = digits.substr(n);
        Size = digits.size();
    } else{
        digits = "0";
        Size = 1;
    }
}

Long Long::operator+(const Long &other) {
    Long res;
    if (sign == other.sign){
        res = add_abs(other);
        res.sign = sign;
    }
    else {
        res = sub_abs(other);
        if (comp_abs(other))
            res.sign = sign;
        else
            res.sign = "+-"[string("+-").find_first_not_of(sign)];
    }
    return res;
}

Long Long::operator-(const Long &other) {
    Long res;
    if (sign != other.sign) {
        res = add_abs(other);
        res.sign = sign;
    }
    else {
        res = sub_abs(other);
        if (comp_abs(other))
            res.sign = sign;
        else
            res.sign = negative_sign();
    }
    return res;
}

Long Long::operator*(const Long &other) {
    Long res;
    if (digits == "0" || other.digits == "0"){
        return res;
    }
    if (sign != other.sign)
        res.sign = '-';
    int N = Size + other.Size;
    res.digits.resize(N, '0');
    res.Size = N;
    int n = Size, m = other.Size;
    int mult = 0;
    int carry = 0;
    for (int i = 1; i <= n; i++){
        carry = 0;
        for (int j = 1; j <= m; j++){
            //cout << res.digit(N - i - j + 1) << " " << digit(n - i) << " " << other.digit(m - j) << " " << carry << "\n";
            mult = res.digit(N - i - j + 1) + digit(n - i) * other.digit(m - j) + carry;
            carry = mult / base;
            mult %= base;
            //cout << mult << " " << carry << "\n";
            res.set(N - i - j + 1, mult);
        }
        res.set(N - i - m, carry);
    }
    res.trim();
    return res;
}

char Long::negative_sign() const {
    if (sign == '+')
        return '-';
    return '-';
}

Long Long::operator<<(int n) {
    for (int i = 0; i < n; i++){
        digits = digits + '0';
    }
    Size += n;
    return *this;
}

Long Long::karatsuba(const Long& other) {
    Long res;
    res.Size = (Size + other.Size);
    res.digits.resize(res.Size, '0');

    if (Size < 4 || other.Size < 4) {
        return *this * other;
    }
    else if (digits == "0" || other.digits == "0") {
        return Long();
    }
    else {
        int m = (max(Size, other.Size) / 2);
        Long a1(digits.substr(0, Size - m));
        Long a0(digits.substr(Size - m));
        Long b1(other.digits.substr(0, other.Size - m));
        Long b0(other.digits.substr(other.Size - m));
        //m = 6;

        Long a0b0 = a0.karatsuba(b0);
        Long a1b1 = a1.karatsuba(b1);
        Long a0a1 = a0 + a1;
        Long b0b1 = b0 + b1;

        res = a0b0 + ((a0a1.karatsuba(b0b1) - a0b0 - a1b1) << m) + (a1b1 << (2 * m));

        //res =  a0.karatsuba(b0) + (((a0 + a1).karatsuba(b0 + b1) - a0.karatsuba(b0) - a1.karatsuba(b1)) << m ) + (a1.karatsuba(b1) << ( 2 * m));
        res.trim();
        return res;
    }
}


