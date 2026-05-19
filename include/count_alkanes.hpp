#pragma once
#ifndef __COUNT_ALKANES_HPP__
#define __COUNT_ALKANES_HPP__

#include <vector>
#include <stdexcept>

// alkyl: generator function A, for alkyls without chiral
// A(x) = 1 + \frac{1}{6} x [A(x) ^ 3 + 3 A(x) A(x^2) + 2 A(x^3)]
// A(x) = \sum_{i=0}^{\infty} a_i x^i
// a_0 = 1

// gen_func_P: generator function P, for all alkanes without chiral with a special carbon atom
// P(x) = \frac{1}{24} x [A(x)^4 + 6 A(x)^2 A(x^2) + 3 A(X^2)^2 + 8 A(x) A(x^3) + 6 A(x^4)]
// p_0 = 0

// gen_func_Q: generator function Q, for all alkanes without chiral with a special carbon-carbon bond
// Q(x) = \frac{1}{2} [A(x)^2 + A(x^2)] - A(x)
// q_0 = q_1 = 0

// gen_func_C: generator function C, for all alkanes without chiral
// C(x) = P(x) - Q(x) + A(x^2) - 1
// c_0 = 0

class Alkanes_counter {
public:
    std::vector<size_t> alkyl, alkane;

    Alkanes_counter();

    void calc_alkyl(int n);

    void calc_alkane(int n);

    size_t get_alkyl(int n);

    size_t get_alkane(int n);
};

#endif // __COUNT_ALKANES_HPP__
