//
// Created by Slava Karpii on 9/23/19.
//

#include "Rand.h"

Long Rand::next() {
    Long xnext = xn * a;
    xn = xnext % m;
    return xn;
}

Rand::Rand(Long seed, Long m_) {
    m = m_;
    xn = seed;
}

