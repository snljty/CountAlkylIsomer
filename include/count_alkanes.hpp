#pragma once
#ifndef __COUNT_ALKANES_HPP__
#define __COUNT_ALKANES_HPP__

#include <vector>
#include <stdexcept>

// A(x): alkyl: generator function A, for alkyls without chiral
// A(x) = 1 + \frac{1}{6} x [A(x)^3 + 3 A(x) A(x^2) + 2 A(x^3)]
// A(x) = \sum_{i=0}^{\infty} a_i x^i
// a_0 = 1

// P(x): generator function P, for all alkanes without chiral with a special carbon atom
// P(x) = \frac{1}{24} x [A(x)^4 + 6 A(x)^2 A(x^2) + 3 A(X^2)^2 + 8 A(x) A(x^3) + 6 A(x^4)]
// P(x) = \sum_{i=0}^{\infty} p_i x^i
// p_0 = 0

// Q(x): generator function Q, for all alkanes without chiral with a special carbon-carbon bond
// Q(x) = \frac{1}{2} [A(x)^2 + A(x^2)] - A(x)
// Q(x) = \sum_{i=0}^{\infty} q_i x^i
// q_0 = q_1 = 0

// C(x): alkane: generator function C, for all alkanes without chiral
// C(x) = P(x) - Q(x) + A(x^2) - 1
// C(x) = \sum_{i=0}^{\infty} c_i x^i
// c_0 = 0

// B(x): stereo alkyl: generator function B, for alkyls with chiral
// B(x) = 1 + \frac{1}{3} x [B(x)^3 + 2 B(x^3)]
// B(x) = \sum_{i=0}^{\infty} b_i x^i

// F(x): generator function F, for all alkanes with chiral with a special carbon atom
// F(x) = \frac{1}{12} x [B(x)^4 + 3 B(X^2)^2 + 8 B(x) B(x^3)]
// F(x) = \sum_{i=0}^{\infty} f_i x^i
// f_0 = 0

// G(x): generator function G, for all alkanes with chiral with a special carbon-carbon bond
// G(x) = \frac{1}{2} [B(x)^2 + B(x^2)] - B(x)
// G(x) = \sum_{i=0}^{\infty} g_i x^i
// g_0 = g_1 = 0

// D(x): alkane: generator function D, for all alkanes with chiral
// D(x) = F(x) - G(x) + B(x^2) - 1
// D(x) = \sum_{i=0}^{\infty} d_i x^i
// d_0 = 0

class Alkanes_counter {
public:
    std::vector<size_t> alkyl, alkane;
    std::vector<size_t> stereo_alkyl, stereo_alkane;
    Alkanes_counter();
    void calc_alkyl(int n);
    void calc_alkane(int n);
    size_t get_alkyl(int n);
    size_t get_alkane(int n);
    void calc_stereo_alkyl(int n);
    void calc_stereo_alkane(int n);
    size_t get_stereo_alkyl(int n);
    size_t get_stereo_alkane(int n);
};

#endif // __COUNT_ALKANES_HPP__
