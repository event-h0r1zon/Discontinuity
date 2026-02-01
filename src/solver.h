#pragma once
#include "riemann.h"
#include <vector>

struct PressureFunctionResult {
    double f_k;
    double f_k_prime;
};

Conserved riemann_solver(const State& left, const State& right);
Conserved hllc_solver(const State& left, const State& right);