#pragma once
#include <vector>
#include <string>

struct State {
    double rho;     // Density
    double u;       // Velocity
    double p;       // Pressure

    double getEnergy() const {
        return p / ((1.4 - 1.0)) + 0.5 * rho * u * u;
    }

    double getInternalEnergy() const {
        return p / ((1.4 - 1.0) * rho) ;
    }
};

struct Conserved {
    double rho;        // Mass density
    double momentum;    // Momentum density
    double energy;      // Energy density
};

// Operator overloads for Conserved struct.

inline Conserved& operator+=(Conserved& a, const Conserved& b) {
    a.rho += b.rho;
    a.momentum += b.momentum;
    a.energy += b.energy;
    return a;
}

inline Conserved& operator-=(Conserved& a, const Conserved& b) {
    a.rho -= b.rho;
    a.momentum -= b.momentum;
    a.energy -= b.energy;
    return a;
}

inline Conserved operator+(Conserved a, const Conserved& b) {
    a += b;
    return a;
}

inline Conserved operator-(Conserved a, const Conserved& b) {
    a -= b;
    return a;
}

inline Conserved operator*(double s, Conserved a) {
    a.rho *= s;
    a.momentum *= s;
    a.energy *= s;
    return a;
}

inline Conserved operator*(Conserved a, double s) {
    return s * a;
}

inline Conserved& operator/=(Conserved& a, double s) {
    a.rho /= s;
    a.momentum /= s;
    a.energy /= s;
    return a;
}

inline Conserved operator/(Conserved a, double s) {
    a /= s;
    return a;
}

std::vector<State> initialize_mesh(int num_cells, const State& left_state, const State& right_state);
std::vector<State> update_states(const std::vector<State>& mesh, const std::vector<Conserved>& fluxes, double dt, double dx);
std::vector<Conserved> calculate_fluxes(const std::vector<State>& mesh, const std::string& solver_type);
double compute_time_step(const std::vector<State>& mesh, double cfl, double dx);
void write_snapshot_csv(const std::vector<State>& mesh, double dx, double time, const std::string& filename);
Conserved calculate_flux(const State& state);
State conserved_to_state(const Conserved& conserved);
Conserved state_to_conserved(const State& state);