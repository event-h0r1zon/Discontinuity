#include "solver.h"
#include "riemann.h"
#include <cmath>
#include <algorithm>

PressureFunctionResult compute_rarefaction_p_functions(double p, const State& state) {
    double a = std::sqrt(1.4 * state.p / state.rho);
    double gamma = 1.4;
    double f_k = (2 * a / (gamma - 1)) * (std::pow(p / state.p, (gamma - 1) / (2 * gamma)) - 1);
    double f_k_prime = (1 / (state.rho * a)) * std::pow(p / state.p, -(gamma + 1) / (2 * gamma));

    return {f_k, f_k_prime};
}

PressureFunctionResult compute_shock_p_functions(double p, const State& state) {
    double gamma = 1.4;
    double A = 2 / ((gamma + 1) * state.rho);
    double B = ((gamma - 1) / (gamma + 1)) * state.p;

    double sqrt_term = std::sqrt(A / (p + B));
    double f_k = (p - state.p) * sqrt_term;
    double f_k_prime = sqrt_term * (1 - 0.5 * (p - state.p) / (p + B));

    return {f_k, f_k_prime};
}

State compute_star_state(const State& state, double p_star, double u_star) {
    State star_state;
    star_state.u = u_star;
    star_state.p = p_star;

    const double gamma = 1.4;

    if (p_star > state.p) {
        const double g1 = (gamma - 1.0) / (gamma + 1.0);
        const double pr = p_star / state.p;
        star_state.rho = state.rho * (pr + g1) / (g1 * pr + 1.0);
    } else {
        star_state.rho = state.rho * std::pow(p_star / state.p, 1.0 / gamma);
    }

    return star_state;
}

State rarefaction_fan_state(const State& state, const double sign) {
    State fan_state;
    const double gamma = 1.4;
    const double xi = 0.0;
    double a = std::sqrt(gamma * state.p / state.rho);
    double gamma_term = (2/(gamma + 1));

    double a_local = gamma_term * (a - sign * 0.5 * (gamma - 1) * (state.u - xi));
    double u_local = gamma_term * (- sign * a + 0.5 * (gamma - 1) * state.u + xi);
    double rho_local = state.rho * std::pow(a_local / a, 2 / (gamma - 1));
    double p_local = state.p * std::pow(a_local / a, 2 * gamma / (gamma - 1));

    fan_state.rho = rho_local;
    fan_state.u = u_local;
    fan_state.p = p_local;

    return fan_state;
}

State sample_interface_solution(const State& star_state, const State& state, const double sign) {
    double p_star = star_state.p;
    double u_star = star_state.u;
    double gamma = 1.4;

    double a_K = std::sqrt(1.4 * state.p / state.rho);

    if (p_star > state.p) {

        double sqrt_term = std::sqrt(((gamma + 1) / (2 * gamma)) * (p_star/state.p) + (gamma - 1) / (2 * gamma));
        double S = state.u + sign * a_K * sqrt_term;
        
        if (sign * S <= 0) return state;
        else return star_state;

    } else {

        double S_H = state.u + sign * a_K;
        double a_star = a_K * std::pow(p_star / state.p, (gamma - 1) / (2 * gamma));
        double S_T = u_star + sign * a_star;

        if (sign * S_H <= 0) return state;
        else if (sign * S_T >= 0) return star_state;
        else return rarefaction_fan_state(state, sign);

    }
}

Conserved riemann_solver(const State& left, const State& right) {
    double a_L = std::sqrt(1.4 * left.p / left.rho);
    double a_R = std::sqrt(1.4 * right.p / right.rho);

    double p_star_guess = 0.5 * (left.p + right.p) - 0.125 * (right.u - left.u) * (left.rho + right.rho) * (a_L + a_R);
    p_star_guess = std::max(1e-6, p_star_guess);

    double p_star_next = 0.0;
    double p_star_curr = p_star_guess;
    
    for (int iter = 0; iter < 100; ++iter) {

        if (iter > 0) 
            p_star_curr = p_star_next;

        PressureFunctionResult f_L_result = (p_star_curr > left.p) 
                                                ? compute_shock_p_functions(p_star_curr, left)
                                                : compute_rarefaction_p_functions(p_star_curr, left);

        double f_L = f_L_result.f_k;
        double f_L_prime = f_L_result.f_k_prime;

        PressureFunctionResult f_R_result = (p_star_curr > right.p)
                                                ? compute_shock_p_functions(p_star_curr, right)
                                                : compute_rarefaction_p_functions(p_star_curr, right);

        double f_R = f_R_result.f_k;
        double f_R_prime = f_R_result.f_k_prime;

        double residual = f_L + f_R + (right.u - left.u);
        p_star_next = p_star_curr - residual / (f_L_prime + f_R_prime);
        
        if (p_star_next < 1e-12) p_star_next = 1e-12;

        const double rel = std::abs(p_star_next - p_star_curr) / (p_star_next + p_star_curr);

        if (rel < 1e-8) break;
    }

    double p_star = p_star_next;

    PressureFunctionResult f_L_result = (p_star > left.p) 
                                                ? compute_shock_p_functions(p_star, left)
                                                : compute_rarefaction_p_functions(p_star, left);

    PressureFunctionResult f_R_result = (p_star > right.p) 
                                                ? compute_shock_p_functions(p_star, right)
                                                : compute_rarefaction_p_functions(p_star, right);

    double f_R = f_R_result.f_k;
    double f_L = f_L_result.f_k;

    double u_star = 0.5 * (left.u + right.u) + 0.5 * (f_R - f_L);

    State star_left, star_right;

    star_left = compute_star_state(left, p_star, u_star);
    star_right = compute_star_state(right, p_star, u_star);

    State interface_state;

    interface_state = (u_star >= 0) 
                        ? sample_interface_solution(star_left, left, -1.0)
                        : sample_interface_solution(star_right, right, 1.0);

    return calculate_flux(interface_state);
}

Conserved hlle_flux(const State& left, const State& right) {
    Conserved U_L = state_to_conserved(left);
    Conserved U_R = state_to_conserved(right);
    Conserved flux_L = calculate_flux(left);
    Conserved flux_R = calculate_flux(right);

    double gamma = 1.4;

    double a_L = std::sqrt(gamma * left.p / left.rho);
    double a_R = std::sqrt(gamma * right.p / right.rho);

    double S_L = std::min(left.u - a_L, right.u - a_R);
    double S_R = std::max(left.u + a_L, right.u + a_R);
    
    if (std::abs(S_R - S_L) < 1e-12) return 0.5 * (flux_L + flux_R);

    Conserved flux_star = (S_R * flux_L - S_L * flux_R + S_L * S_R * (U_R - U_L)) / (S_R - S_L);

    if (0 <= S_L) return flux_L;
    else if (S_L <= 0 && 0 <= S_R) return flux_star;
    else return flux_R;
}

Conserved hllc_solver(const State& left, const State& right) {
    Conserved U_L = state_to_conserved(left);
    Conserved U_R = state_to_conserved(right);
    Conserved flux_L = calculate_flux(left);
    Conserved flux_R = calculate_flux(right);

    double gamma = 1.4;
    double eps = 1e-12;

    double a_L = std::sqrt(gamma * left.p / left.rho);
    double a_R = std::sqrt(gamma * right.p / right.rho);

    double S_L = std::min(left.u - a_L, right.u - a_R);
    double S_R = std::max(left.u + a_L, right.u + a_R);

    double A_L = left.rho * (S_L - left.u);
    double A_R = right.rho * (S_R - right.u);

    double denom_S_star = (A_L - A_R);

    if (std::abs(denom_S_star) < eps) return hlle_flux(left, right);

    double S_star = (right.p - left.p + A_L * left.u - A_R * right.u) / denom_S_star;
    
    double denom_L = (S_L - S_star);
    double denom_R = (S_R - S_star);

    if (std::abs(denom_L) < eps || std::abs(denom_R) < eps) return hlle_flux(left, right);

    double p_star_L = left.p + left.rho * (S_L - left.u) * (S_star - left.u);
    double p_star_R = right.p + right.rho * (S_R - right.u) * (S_star - right.u);

    double p_star = 0.5 * (p_star_L + p_star_R);

    double rho_star_L = left.rho * (S_L - left.u) / (S_L - S_star);
    double rho_star_R = right.rho * (S_R - right.u) / (S_R - S_star);

    double momentum_star_L = rho_star_L * S_star;
    double momentum_star_R = rho_star_R * S_star;

    double energy_star_L = ((S_L - left.u) * U_L.energy - left.p * left.u + p_star * S_star) / denom_L;
    double energy_star_R = ((S_R - right.u) * U_R.energy - right.p * right.u + p_star * S_star) / denom_R;

    Conserved U_star_L = {rho_star_L, momentum_star_L, energy_star_L};
    Conserved U_star_R = {rho_star_R, momentum_star_R, energy_star_R};

    Conserved flux_star_L = flux_L + S_L * (U_star_L - U_L);
    Conserved flux_star_R = flux_R + S_R * (U_star_R - U_R);

    if (0 <= S_L) return flux_L;
    else if (S_L <= 0 && 0 <= S_star) return flux_star_L;
    else if (S_star <= 0 && 0 <= S_R) return flux_star_R;
    else return flux_R;
}