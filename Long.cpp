//
// Created by Slava Karpii on 9/13/19.
//

#include <iostream>
#include <cmath>
#include<ctime>
#include "Long.h"
#include "Rand.h"

test* Long::t = new class test();

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

bool Long::operator==(const Long& other) const {
    return digits == other.digits && sign == other.sign;
}

bool Long::operator>(const Long& other) const {
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

bool Long::operator>=(const Long& other) const {
    return *this > other || *this == other;
}

bool Long::operator<(const Long& other) const {
    return !(*this >= other);
}

bool Long::operator<=(const Long& other) const {
    return !(*this > other);
}

bool Long::operator!=(const Long& other) const {
    return !(*this == other);
}

Long Long::add_abs(const Long &other) const {
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

Long Long::sub_abs(const Long &other) const{
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
    for (i = 1; (i <= n || carry != 0); i++){
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

Long Long::operator+(const Long &other) const{
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
    if (res.digits == "0")
        res.sign = '+';
    return res;
}

Long Long::operator-(const Long &other) const{
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
    if (res.digits == "0")
        res.sign = '+';
    return res;
}

Long Long::operator*(const Long &other) const{
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
    return '+';
}

Long Long::operator<<(int n) {
    if (digits != "0") {
        if (n > 0) {
            for (int i = 0; i < n; i++) {
                digits = digits + '0';
            }
            Size += n;
        } else if (Size + n > 0) {
            digits = digits.substr(0, Size + n);
            Size = digits.size();
        } else {
            Size = 1;
            digits = "0";
        }
    }
    return *this;
}

Long Long::karatsuba(const Long& other) const{
    Long res;
    res.Size = (Size + other.Size);
    res.digits.resize(res.Size, '0');

    if (Size < 4 || other.Size < 4 || Size + other.Size < 8) {
        return *this * other;
    }
    else {
        int m = (max(Size, other.Size) / 2);
        const Long a1 = (*this << -m);
        Long a0 = *this - (a1 << m);

        const Long b1 = (other << -m);
        Long b0 = other - (b1 << m);
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

Long::Long(int n) {
    sign = '+';
    if (n < 0) {
        sign = '-';
        n *= -1;
    }
    digits = to_string(n);
    Size = digits.size();
}

Long Long::toom3(const Long &other) const{
    Long res;
    res.Size = (Size + other.Size);
    res.digits.resize(res.Size, '0');

    if (Size < 4 || other.Size < 4 || Size + other.Size < 8) {
        return *this * other;
    }
    else {
        int k = 3;
        Long m[3];
        Long n[3];
        int B = max(floor(Size / k), floor(other.Size / k)) + 1;
        //cout << B<<"\n";

        m[2] = (*this << (- 2 * B));
        m[1] = (*this - ((*this << (- 2 * B)) << (2 * B))) << -B;
        m[0] = *this - ((*this << (-B)) << (B));

        n[2] = (other << (- 2 * B));
        n[1] = (other - ((other << (- 2 * B)) << (2 * B))) << -B;
        n[0] = other - ((other << (-B)) << (B));

        Long p[5];
        Long q[5];

        int vals[5][3] = {
                {1, 0, 0},
                {1, 1, 1},
                {1, -1, 1},
                {1, -2, 4},
                {0, 0, 1}
        };
        for (int i = 0; i < 5; i++){
            for (int j = 0; j < 3; j++){
                p[i] = p[i] + m[j] * vals[i][j];
                q[i] = q[i] + n[j] * vals[i][j];
            }
        }

//        for (int i = 0; i < 5; i++){
//            cout << q[i] << ' ';
//        }
       // cout << '\n';

        Long rv[5];
        for (int i = 0; i < 5; i++){
            rv[i] = p[i].toom3(q[i]);
        }
        Long r[5];
        r[0] = rv[0];
        r[4] = rv[4];
        r[3] = (rv[3] - rv[1]) / 3;
        r[1] = (rv[1] - rv[2]) / 2;
        r[2] = rv[2] - rv[0];
        r[3] = (r[2] - r[3]) / 2 + rv[4] * 2;
        r[2] = r[2] + r[1] - r[4];
        r[1] = r[1] - r[3];

        for (int i = 0; i < 5; i++){
           // cout << r[i] << ' ';
        }
        //cout << '\n';
        for (int i = 0; i < 5; i++){
            //cout << (r[i] << (5 - i)) << ' ';
        }
        //cout << '\n';

        for (int i = 0; i < 5; i++){
            res = res + (r[i] << (B * i));
        }
        res.trim();
        return res;
    }
}

Long Long::operator<<(int n) const {
    Long res = *this;
    return res << n;
}

Long Long::operator*(int n) const{
    Long res;
    bool neg = (n < 0);
    n = abs(n);
    int m = Size;
    for (int i = 1; i <= m; i++){
       // const Long s = Long(digit(m - i) * n);
        //cout << digit(m - i)  << " " << s << " " << (s << (i- 1)) << "\n";
        res = res +  (Long(digit(m - i) * n) << (i - 1));
    }
    if (!neg xor (sign == '+'))
        res.sign = '-';
    //cout << "------\n";
    return res;
}

Long Long::operator/(int n) const {
    Long res;
    bool neg = (n < 0);
    n = abs(n);
    int m = Size;
    res.digits.resize(m, '0');
    res.Size = m;
    int rem = 0;
    for (int i = 0; i < m; i++){
        rem = rem * 10 + digit(i);
        //rem /= n;
        res.set(i, rem / n);
        rem %= n;
    }
    if (!neg xor (sign == '+'))
        res.sign = '-';
    res.trim();
    return res;
}

Long Long::operator%(const Long &other) const {
    return (*this - (*this / other) * other);
}

Long Long::pow_mod(const Long &other, const Long& modulus) const {
//    Long one = Long(1);
//    Long res = one;
//    Long count = other;
//    while (count > 0){
//        res = res * *this;
//        res = res % modulus;
//        count = count - one;
//    }
//    return res;
    if (modulus == 1)
        return 0;
    Long r = 1;
    Long b = *this % modulus;
    Long e = other;
    while (e > 0){
        if (e % 2 == 1)
            r = (r * b) % modulus;
        e = e / 2;
        b = (b * b) % modulus;
    }
    return r;
}

bool Long::fermat() const {
    if (*this == 2)
        return true;
    if (*this % 2 == 0)
        return false;
    bool prime = true;
    srand(time(0));
    Rand r(rand() % 500, *this);
    int k = Size;
    while (k > 0 && prime){
        k--;
        Long test = r.next();
        test = test.pow_mod(*this - 1, *this);
        //cout << test << ' ';
        if (test != 1)
            prime = false;
    }
    return prime;
}

Long Long::operator/(const Long &other) const {
    if (other.digits == "0")
        return *this;
    Long res;
    if (*this < other)
        return res;
    int n = other.Size;
    int it = Size - other.Size;
    Long tmp = *this << (-Size + n);
    int ind = tmp.Size;
    while (it >= 0){
        //cout << tmp << ' ';
        int k = bin_search(0, base - 1, tmp, other);

        res.digits += to_string(k);
        tmp = tmp - other * k;
        //cout << tmp << ' ';
        tmp.digits += digits[ind];
        ind++;
       // cout << tmp << '\n';
        tmp.Size++;
        tmp.trim();
        it--;
    }

//    while (tmp >= other) {
//        int off = tmp.Size - n - rem;
//        tmp << (-off);
//        int k = bin_search(0, base - 1, tmp, other);
//        cout << k << '\n';
//        res.digits = res.digits + to_string(k);
//        if (k == 0)
//            rem = 1;
//        else
//            rem = 0;
//        res.trim();
//        Long g = other * res;
//        tmp = *this - (g << (Size - g.Size)) ;
//    }
    if (sign != other.sign)
        res.sign = '-';
    res.trim();
    return res;
}

int Long::bin_search(int a, int b, const Long &M, const Long &n) const {
    if (a == b){
        return a;
    }
    int c = (a + b) / 2 + 1;
    if ( n * c > M)
        return bin_search(a, c - 1, M, n);
    else
        return bin_search(c, b, M, n);
}

bool Long::sol_strass() {
    if (*this == 2)
        return true;
    if (*this % 2 == 0) {
        cout << "f";
        return false;
    }
    int k = Size;
    bool prime = true;
    srand(time(0));
    Rand r(rand() % 500, *this);
    while (k > 0 && prime){
        k--;
        Long a = r.next();
        while (a == 1 || a == 0){
            a = r.next();
        }
        Long x = Jacobi(a, *this);
        //cout << a << "a";
        if (x == -1)
            x = x + *this;
        //cout << x << "x";
        Long y = a.pow_mod((*this - 1) / 2, *this);
        //cout << y << "y";
        if (x == 0 || x != y)
            prime = false;
    }
    return prime;
}

Long Long::operator%(int n) const {
    return (*this - (*this / n) * n);
}

bool Long::rab_mil() const {
    if (*this == 2)
        return true;
    if (*this % 2 == 0) {
        cout << "f";
        return false;
    }
    int r = 0;
    Long d = *this - 1;
    while (d % 2 == 0) {
        d = d / 2;
        r++;
    }
    int k = Size;
    bool prime = true;
    srand(time(0));
    Rand rnd(rand() % 500, *this);
    while(k > 0 && prime) {
        k--;
        Long a = rnd.next();
        while (a == 1 || a == 0){
            a = rnd.next();
        }
        Long x =a.pow_mod(d, *this);
        if (x == 1 || x == *this - 1)
            continue;
        for (int i = 0; i < r - 1; i++) {
            x = x.pow_mod(2, *this);
            if (x != *this - 1)
                prime = false;
        }
    }
    return prime;
}

Long GCD(const Long& first, const Long& other){
    if (other == 0)
        return first;
    return GCD(other, first % other);
}


ostream &operator<<(ostream &out, Long i) {
    out << string(i);
    return out;
}

Long Jacobi(Long n, Long k) {
    n = n % k;
    Long t = 1;
    while (n != 0) {
        while (n % 2 == 0) {
            n = n / 2;
            Long r = k % 8;
            if (r == 3 || r == 5)
                t = t * -1;
        }
        Long tmp = n;
        n = k;
        k = tmp;
        if (n % 4 == 3 && k % 4 == 3)
            t = t * -1;
        n = n % k;
    }
    if (k == 1)
        return t;
    else
        return 0;
}
