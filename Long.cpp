//
// Created by Slava Karpii on 9/13/19.
//

#include <iostream>
#include <cmath>
#include<ctime>
#include "Long.h"
#include "Rand.h"

int Long::base = 10;
Mult* Long::t = new Mult();

void Long::bswitch(int b) {
    base = b;
}

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
    binary();
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
        if (s <  0){
            s += base;
            carry = 1;
        }
        res.set(N - i, s);
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
    cout << typeid(*t).name() << endl;
    return t -> mult(*this, other);
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


Long::Long(int n) {
    sign = '+';
    if (n < 0) {
        sign = '-';
        n *= -1;
    }
    digits = to_string(n);
    Size = digits.size();
    binary();
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
        res = res +  (Long(digit(m - i) * n) << (i - 1));
    }
    if (!neg xor (sign == '+'))
        res.sign = '-';
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
        res.set(i, rem / n);
        rem %= n;
    }
    if (!neg xor (sign == '+'))
        res.sign = '-';
    res.trim();
    return res;
}

Long Long::operator%(const Long &other) const {
    Long res = abs(*this) - (abs(*this) / other).naive(other);
    if (sign == '-')
        res.sign = '-';
    return res;
}

Long Long::pow_mod(const Long &other, const Long& modulus) const {
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
        int k = bin_search(0, base - 1, tmp, other);
        res.digits += to_string(k);
        tmp = tmp - other * k;
        tmp.digits += digits[ind];
        ind++;
        tmp.Size++;
        tmp.trim();
        it--;
    }
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
        if (x == -1)
            x = x + *this;
        Long y = a.pow_mod((*this - 1) / 2, *this);
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

void Long::binary()  {
    Long tmp = *this;
    while (tmp.digits > "0") {
        char dig = tmp.digits[tmp.Size - 1];
        if (string("13579").find(dig) != string::npos)
            dig = '1';
        else
            dig = '0';
        bin = dig + bin;
        tmp = tmp / 2;
    }
}

Long Long::pow(const Long &other) const {
        Long r = 1;
        Long b = *this;
        Long e = other;
        while (e > 0){
            if (e % 2 == 1)
                r = (r.naive(b));
            e = e / 2;
            b = (b.naive(b));
        }
        return r;
}


string Long::cook(bool decimal ) const{
    int b = base;
    bswitch(2);
    if (Size < 2) {
        int a = stoi(digits);
        if (sign == '-')
            a *= -1;
        cout << 1 / a << '\n';
        string s = to_string(1 / a);
        s = s.substr(s.find('.') + 1);
        Long res = Long(s);
        if (decimal)
            return s;
        return res.bin;
    }
    Long z;
    int zd; //number of digits after floating point
    double a = floor(32 / (4 * (int(bin[0]) - '0') + 2 * (int(bin[1]) - '0') + (int(bin[2]) - '0'))) * 0.25;
   if (bin[1] == '0') {
       if (bin[2] == '0') {
           z = Long("10");
           zd =0;
       }
       else {
           z = Long("11");
           zd = 1;
       }
   }
   else {
       if (bin[2] == '0') {
           z = Long("101");
           zd = 1;
       }
       else {
           z = Long("1");
           zd = 0;
       }
   }


   int n = bin.size();
   int k = 0;
   while (::pow(2, k) < n) {
       Long zz = z * z;
       int zzd = 2 * zd;

       int Vkd = ::pow(2, k + 1) + 3;
       Long Vk = Long(bin.substr(0,Vkd));

       Long Vzz = Vk * zz;
       int Vzzd = Vkd + zzd;

       Long z2 = z + z;
       int z2d = zd;

       if (z2d > Vzzd) {
           Vzz << (z2d - Vzzd);
           Vzzd = z2d;
       }
       else if (z2d < Vzzd){
           z2 << (-z2d + Vzzd);
           z2d = Vzzd;
       }
       if (Vzz.Size < Vzzd){
           for (int i = 0; i < Vzzd - Vzz.Size; i++)
               Vzz.set(-1, 0);
       }

       if (z2.Size < z2d){
           for (int i = 0; i < z2d - z2.Size; i++)
               z2.set(-1, 0);
       }

       z = z2 - Vzz;
       zd = z2d;
       while (z.digit(z.Size - 1) == 0) {
           z << (-1);
           zd--;
       }

       k++;
   }
//   cout << "0.";
   string s;
   for (int i = 0; i < bin.size() - 1; i++) {
      // cout << '0';
       s = s + '0';
   }
   //cout << z << '\n';
   s += z.digits;
   bswitch(b);
   if (!decimal)
       return s;
   else {
       Long ten = Long(1) << s.size();
       Long decres;
       for (int i = 1; i <= s.size(); i++) {
           decres = decres + (ten / (Long(2).pow(i))) * int(s[i - 1] - '0');
       }
       //cout << "0.";
       string decs;
       for (int i = 0; i < s.size() - decres.Size; i++) {
           decs += '0';
           //cout << '0';
       }
       decs += decres.digits;
       //cout << decres << endl;
       return decs;
   }
}

Long Long::div_cook(const Long &other) {
    string cook = other.cook(false);
    Long inv = Long(cook);
    int n = cook.size();
    int b = base;
    bswitch(2);
    Long res = Long(bin) * inv;
    res << (-n);
    if (sign == '+' xor other.sign == '+')
        res.sign = '-';
    bswitch(b);
    Long decres;
    for (int i = 1; i <= res.Size; i++) {
        decres = decres + Long(2).pow(res.Size - i) * res.digit(i - 1) ;
    }
    return decres;
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

istream &operator>>(istream &in, Long &a) {
    string s;
    in >> s;
    a = s;
    return in;
}

Long abs(const Long &a) {
    Long res = a;
    res.sign = '+';
    return res;
}

Long Long::naive(const Long &other) const {
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
            mult = res.digit(N - i - j + 1) + digit(n - i) * other.digit(m - j) + carry;
            carry = mult / base;
            mult %= base;
            res.set(N - i - j + 1, mult);
        }
        res.set(N - i - m, carry);
    }
    res.trim();
    return res;
}

Long Toom3::toom3(const Long &a, const Long &other) const {
    Long res;
    res.Size = (a.Size + other.Size);
    res.digits.resize(res.Size, '0');

    if (a.Size < 4 || other.Size < 4 || a.Size + other.Size < 8) {
        return naive(a, other);
    }
    else {
        int k = 3;
        Long m[3];
        Long n[3];
        int B = max(floor(a.Size / k), floor(other.Size / k)) + 1;

        m[2] = (a << (- 2 * B));
        m[1] = (a - ((a << (- 2 * B)) << (2 * B))) << -B;
        m[0] = a - ((a << (-B)) << (B));

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

        Long rv[5];
        for (int i = 0; i < 5; i++){
            rv[i] = toom3(p[i], q[i]);
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
            res = res + (r[i] << (B * i));
        }
        res.trim();
        return res;
    }
}

Long Karatsuba::karatsuba(const Long &a, const Long &other) const {
    Long res;
    res.Size = (a.Size + other.Size);
    res.digits.resize(res.Size, '0');

    if (a.Size < 4 || other.Size < 4 || a.Size + other.Size < 8) {
        return naive(a, other);
    }
    else {
        int m = (max(a.Size, other.Size) / 2);
        const Long a1 = (a << -m);
        Long a0 = a - (a1 << m);

        const Long b1 = (other << -m);
        Long b0 = other - (b1 << m);

        Long a0b0 = karatsuba(a0, b0);
        Long a1b1 = karatsuba(a1, b1);
        Long a0a1 = a0 + a1;
        Long b0b1 = b0 + b1;

        res = a0b0 + ((karatsuba(a0a1, b0b1) - a0b0 - a1b1) << m) + (a1b1 << (2 * m));

        res.trim();
        return res;
    }
}

Long Mult::naive(const Long &a, const Long &other) const {
    return a.naive(other);
}

Long Sen_strass::sen_strass(const Long &a, const Long &other) const {
    int n = 1;
    int s = max(a.Size, other.Size);
    while (n < s) {
        n <<= 1;
    }
    n <<= 1;
    string first;
    first.resize( (n - a.Size % n) % n, '0');
    first += a.digits;
    string second;
    second.resize((n - other.Size % n) % n, '0');
    second += other.digits;
    int k1 = first.size() / n;
    int k2 = second.size() / n;
    vector<b> sig1(n);
    vector<b> sig2(n);
    for (int i = 1; i <= n; i++){
        sig1[i - 1] = stoi(first.substr(first.size() - k1 * i, k1));
        sig2[i - 1] = stoi(second.substr(second.size() - k2 * i, k2));
    }
    fft(sig1, false);
    fft(sig2, false);
    for (int i = 0; i < n; i++){
        sig1[i] *= sig2[i];
    }
    fft(sig1, true);
    vector<int> res(n);
    for (int i = 0; i < n; i++) {
        res[i] = int(sig1[i].real() + 0.5);
    }
    int carry = 0;
    for (int i = 0; i < n; i++) {
        res[i] += carry;
        carry = res[i] / Long::base;
        res[i] %= Long::base;
    }
    string re;
    for (int i = 0; i < n; i++) {
        re = to_string(res[i]) + re;
    }
    return re;
}

void Sen_strass::fft(vector<b> &a, bool invert) const {
    int n = (int) a.size();
    if (n == 1)  return;

    vector<b> a0 (n/2),  a1 (n/2);
    for (int i = 0, j = 0; i < n; i += 2, ++j) {
        a0[j] = a[i];
        a1[j] = a[i+1];
    }
    fft (a0, invert);
    fft (a1, invert);

    double ang = 2 * M_PI / n * (invert ? -1 : 1);
    b w (1),  wn (cos(ang), sin(ang));
    for (int i = 0; i < n / 2; ++i) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        if (invert)
            a[i] /= 2,  a[i + n / 2] /= 2;
        w *= wn;
    }
}

Long Modular::modular(const Long &a, const Long &other) const {
    if (a.Size < 4 || other.Size < 4 || a.Size + other.Size < 8) {
        return naive(a, other);
    } else {
        int qk = 1;
        int pk = 18 * qk + 8;
        int k = 0;
        Long x = a;
        Long y = other;
        int mlen = max(x.bin.size(), y.bin.size());
        while (pk < mlen) {
            k++;
            qk = q(k);
            pk = 18 * qk + 8;
        }
        int xoff = pk - x.bin.size();
        int yoff = pk - y.bin.size();
        for (int i = 0; i < xoff; i++)
            x = x * int(2);
        for (int i = 0; i < yoff; i++)
            y = y * int(2);

        vector<Long> mods = m(qk);
        vector<Long> u(6);
        vector<Long> v(6);
        vector<Long> w(6);
        for (int i = 0; i < 6; i++) {
            u[i] = x % mods[i];
            v[i] = y % mods[i];
            w[i] = naive(u[i], v[i]);
        }

        Long r[6][6];
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                r[i][j] = inv_mod(mods[i], mods[j]);
            }
        }

        vector<Long> z(6);
        for (int i = 0; i < 6; ++i) {
            z[i] = w[i];
            for (int j = 0; j < i; ++j) {
                z[i] = naive(r[j][i], (z[i] - z[j]));
                z[i] = z[i] % mods[i];
                if (z[i].sign == '-')
                    z[i] = z[i] + mods[i];
            }
        }

        Long res;
        for (int i = 0; i < 6; i++) {
            Long p = 1;
            for (int j = 0; j < i; j++) {
                p = naive(p, mods[j]);
            }
            res = res + naive(z[i], p);
        }

        for (int i = 0; i < xoff; i++)
            res = res / int(2);
        for (int i = 0; i < yoff; i++)
            res = res / int(2);

        return res;
    }
}

Long Modular::gcdex(Long a, Long b, Long &x, Long &y) const{
    if (a == 0) {
        x = 0; y = 1;
        return b;
    }
    Long x1, y1;
    Long d = gcdex (b % a, a, x1, y1);
    x = y1 - naive((b / a), x1);
    y = x1;
    return d;
}

Long Modular::inv_mod(Long a, Long m) const{
    Long x, y;
    Long g = gcdex (a, m, x, y);
    x = (x % m + m) % m;
    return x;
}

int Modular::q(int k) const {
    if (k == 0)
        return 1;
    return 3 * q(k - 1) - 1;
}

vector<Long> Modular::m(int qk) const {
    vector<Long> mods (6);
    int pows[6];
    pows[0] = 6 * qk - 1;
    pows[1] = 6 * qk + 1;
    pows[2] = 6 * qk + 2;
    pows[3] = 6 * qk + 3;
    pows[4] = 6 * qk + 5;
    pows[5] = 6 * qk + 7;
    for (int i = 0; i < 6; i++) {
        mods[i] = (Long(2).pow(pows[i])) - 1;
    }
    return mods;
}
