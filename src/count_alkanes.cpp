#include "count_alkanes.hpp"

Alkanes_counter::Alkanes_counter() : 
    alkyl(1, 1), alkane(1, 0), stereo_alkyl(1, 1), stereo_alkane(1, 0) {};

void Alkanes_counter::calc_alkyl(int n) {
    if (alkyl.size() > n) return;
    if (!n) return;
    calc_alkyl(n - 1);
    size_t a_n = 0;
    // contribution of A(x) ^ 3
    for (int i = 0; i < n; ++ i) {
        for (int j = 0; i + j < n; ++ j) {
            int k = n - 1 - i - j;
            a_n += alkyl[i] * alkyl[j] * alkyl[k];
        }
    }
    // contribution of 3 A(x) A(x^2)
    for (int j = 0; 2 * j < n; ++ j) {
        int i = n - 1 - 2 * j;
        a_n += 3 * alkyl[i] * alkyl[j];
    }
    // contribution of 2 A(x^3)
    if (n % 3 == 1) {
        int i = (n - 1) / 3;
        a_n += 2 * alkyl[i];
    }
    // |S_3| = 3! = 6
    if (a_n % 6 != 0) throw std::runtime_error("Error: A(x) calculation wrong.");
    a_n /= 6;
    alkyl.push_back(a_n);
}

void Alkanes_counter::calc_stereo_alkyl(int n) {
    if (stereo_alkyl.size() > n) return;
    if (!n) return;
    calc_stereo_alkyl(n - 1);
    size_t b_n = 0;
    // contribution of B(x) ^ 3
    for (int i = 0; i < n; ++ i) {
        for (int j = 0; i + j < n; ++ j) {
            int k = n - 1 - i - j;
            b_n += stereo_alkyl[i] * stereo_alkyl[j] * stereo_alkyl[k];
        }
    }
    // contribution of 2 B(x^3)
    if (n % 3 == 1) {
        int i = (n - 1) / 3;
        b_n += 2 * stereo_alkyl[i];
    }
    // |A_3| = 3!/2 = 3
    if (b_n % 3 != 0) throw std::runtime_error("Error: B(x) calculation wrong.");
    b_n /= 3;
    stereo_alkyl.push_back(b_n);
}

void Alkanes_counter::calc_alkane(int n) {
    if (alkane.size() > n) return;
    if (!n) return;
    alkane.resize(n + 1);
    // all unique carbon atoms in all isomers
    size_t p_n = 0;
    // contribution of A(x)^4
    for (int i = 0; i < n; ++ i) {
        for (int j = 0; i + j < n; ++ j) {
            for (int k = 0; i + j + k < n; ++ k) {
                int l = n - 1 - i - j - k;
                p_n += alkyl[i] * alkyl[j] * alkyl[k] * alkyl[l];
            }
        }
    }
    // contribution of 6 A(x)^2 A(x^2)
    for (int k = 0; 2 * k < n; ++ k) {
        for (int i = 0; i + 2 * k < n; ++ i) {
            int j = n - 1 - i - 2 * k;
            p_n += 6 * alkyl[i] * alkyl[j] * alkyl[k];
        }
    }
    // contribution of 3 A(x^2)^2
    if (n % 2 == 1) {
        for (int i = 0; 2 * i < n; ++ i) {
            int j = (n - 1 - 2 * i) / 2;
            p_n += 3 * alkyl[i] * alkyl[j];
        }
    }
    // contribution of 8 A(x) A(x^3)
    for (int j = 0; 3 * j < n; ++ j) {
        int i = n - 1 - 3 * j;
        p_n += 8 * alkyl[i] * alkyl[j];
    }
    // contribution of 6 A(x^4)
    if (n % 4 == 1) {
        int i = (n - 1) / 4;
        p_n += 6 * alkyl[i];
    }
    // |S_4| = 4! = 24
    if (p_n % 24 != 0) throw std::runtime_error("Error: P(x) calculation wrong.");
    p_n /= 24;

    // all unique carbon-carbon bonds in all isomers
    size_t q_n = 0;
    // contribution of A(x)^2
    for (int i = 0; i <= n; ++ i) {
        int j = n - i;
        q_n += alkyl[i] * alkyl[j];
    }
    // contribution of A(x^2)
    if (n % 2 == 0) {
        int i = n / 2;
        q_n += alkyl[i];
    }
    // |S_2| = 2! = 2
    if (q_n % 2 != 0)  throw std::runtime_error("Error: Q(x) calculation wrong.");
    q_n /= 2;
    // contribution of - A(x)
    q_n -= alkyl[n];

    // C(x) = P(x) - Q(x) + A(x^2) - 1
    alkane[n] = p_n - q_n + (n % 2 == 0 ? alkyl[n / 2] : 0);
}

void Alkanes_counter::calc_stereo_alkane(int n) {
    if (stereo_alkane.size() > n) return;
    if (!n) return;
    stereo_alkane.resize(n + 1);
    // all unique carbon atoms in all isomers
    size_t f_n = 0;
    // contribution of B(x)^4
    for (int i = 0; i < n; ++ i) {
        for (int j = 0; i + j < n; ++ j) {
            for (int k = 0; i + j + k < n; ++ k) {
                int l = n - 1 - i - j - k;
                f_n += stereo_alkyl[i] * stereo_alkyl[j] * stereo_alkyl[k] * stereo_alkyl[l];
            }
        }
    }
    // contribution of 3 B(x^2)^2
    if (n % 2 == 1) {
        for (int i = 0; 2 * i < n; ++ i) {
            int j = (n - 1 - 2 * i) / 2;
            f_n += 3 * stereo_alkyl[i] * stereo_alkyl[j];
        }
    }
    // contribution of 8 B(x) B(x^3)
    for (int j = 0; 3 * j < n; ++ j) {
        int i = n - 1 - 3 * j;
        f_n += 8 * stereo_alkyl[i] * stereo_alkyl[j];
    }
    // |A_4| = 4!/2 = 12
    if (f_n % 12 != 0) throw std::runtime_error("Error: F(x) calculation wrong.");
    f_n /= 12;

    // all unique carbon-carbon bonds in all isomers
    size_t g_n = 0;
    // contribution of B(x)^2
    for (int i = 0; i <= n; ++ i) {
        int j = n - i;
        g_n += stereo_alkyl[i] * stereo_alkyl[j];
    }
    // contribution of B(x^2)
    if (n % 2 == 0) {
        int i = n / 2;
        g_n += stereo_alkyl[i];
    }
    // |S_2| = 2! = 2
    if (g_n % 2 != 0)  throw std::runtime_error("Error: G(x) calculation wrong.");
    g_n /= 2;
    // contribution of - B(x)
    g_n -= stereo_alkyl[n];

    // C(x) = P(x) - Q(x) + B(x^2) - 1
    stereo_alkane[n] = f_n - g_n + (n % 2 == 0 ? stereo_alkyl[n / 2] : 0);
}

size_t Alkanes_counter::get_alkyl(int n) {
    calc_alkyl(n);
    return alkyl[n];
}

size_t Alkanes_counter::get_alkane(int n) {
    calc_alkane(n);
    return alkane[n];
}

size_t Alkanes_counter::get_stereo_alkyl(int n) {
    calc_stereo_alkyl(n);
    return stereo_alkyl[n];
}

size_t Alkanes_counter::get_stereo_alkane(int n) {
    calc_stereo_alkane(n);
    return stereo_alkane[n];
}
