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
};

struct Conserved {
    double rho;        // Mass density
    double momentum;    // Momentum density
    double energy;      // Energy density
};

std::vector<State> initialize_mesh(int num_cells, const State& left_state, const State& right_state);
std::vector<State> update_states(const std::vector<State>& mesh, const std::vector<Conserved>& fluxes, double dt, double dx);
std::vector<Conserved> calculate_fluxes(const std::vector<State>& mesh);
double compute_time_step(const std::vector<State>& mesh, double cfl, double dx);
void write_snapshot_csv(const std::vector<State>& mesh, double dx, double time, const std::string& filename);
Conserved calculate_flux(const State& state);