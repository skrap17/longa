//
// Created by Slava Karpii on 9/23/19.
//

#ifndef LONGA_RAND_H
#define LONGA_RAND_H

#include "Long.h"


class Rand {
    Long xn;
    int a = 11;
    Long c = Long(1);
    Long m;
public:
    Rand(Long seed = Long("2"), Long m = Long("100000000000000"));
    Long next();
};


#endif //LONGA_RAND_H
