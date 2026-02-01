#include <iostream>
#include "riemann.h"
#include "solver.h"
#include <vector>
#include <sstream>
#include <iomanip>
#include <chrono>
#include "tests.cpp"

int main() {
    for (auto& entry : std::filesystem::directory_iterator("../tests")) {
        std::cout << "Processing test file: " << entry.path() << std::endl;
        if (!entry.is_regular_file()) continue;
        auto tc = parse_test_file(entry.path());
        std::cout << "Running test: " << tc.name << std::endl;
        run_test(tc, "../data");
        std::cout << "Test " << tc.name << " completed." << std::endl;
        std::cout << "----------------------------------------" << std::endl;
    }
    return 0;
}