#include "count_alkanes.hpp"

#include <iostream>
#include <string>
#include <iomanip>

int main(int argc, const char *argv[]) {
    int max_n_count = 20;
    if (argc - 1 == 1) {
        max_n_count = std::atoi(argv[1]);
    }

    Alkanes_counter counter;

    std::string label;

    std::cout << std::left;

    std::cout << "Alkyl:" << std::endl;
    for (int n = 1; n <= max_n_count; ++ n) {
        int m = 2 * n + 1;
        if (n == 1) {
            label = "CH" + std::to_string(m) + "-:";
        } else {
            label = "C" + std::to_string(n) + "H" + std::to_string(m) + "-:";
        }
        std::cout << std::setw(9) << label << counter.get_alkyl(n) << '\n';
    }
    std::cout << std::endl;

    std::cout << "Alkane:" << std::endl;
    for (int n = 1; n <= max_n_count; ++ n) {
        int m = 2 * n + 2;
        if (n == 1) {
            label = "CH" + std::to_string(m) + ":";
        } else {
            label = "C" + std::to_string(n) + "H" + std::to_string(m) + ":";
        }
        std::cout << std::setw(8) << label << counter.get_alkane(n) << '\n';
    }

    return 0;
}

