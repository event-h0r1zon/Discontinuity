#include "riemann.h"
#include "solver.h"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <filesystem>

std::vector<State> initialize_mesh(int num_cells, const State& left_state, const State& right_state) {
    
    if (num_cells <= 0) {
        throw std::invalid_argument("Number of cells must be positive.");
    }

    if (num_cells % 2 != 0) {
        throw std::invalid_argument("Number of cells must be even for symmetric initialization.");
    }

    std::vector<State> mesh(num_cells);
    int mid = num_cells / 2;

    for (int i = 0; i < num_cells; ++i) {
        if (i < mid) {
            mesh[i] = left_state;
        } else {
            mesh[i] = right_state;
        }
    }

    return mesh;
}

Conserved state_to_conserved(const State& state) {
    Conserved conserved;

    conserved.rho = state.rho;
    conserved.momentum = state.rho * state.u;
    conserved.energy = state.getEnergy();

    return conserved;
}

State conserved_to_state(const Conserved& conserved) {
    State state;

    state.rho = conserved.rho;
    state.u = conserved.momentum / conserved.rho;
    state.p = (1.4 - 1.0) * (conserved.energy - 0.5 * conserved.momentum * state.u);

    return state;
}

std::vector<Conserved> calculate_fluxes(const std::vector<State>& mesh) {
    int N = mesh.size();
    
    std::vector<Conserved> fluxes(N - 1);

    for (int i = 0; i < N - 1; i++) {
        State left = mesh[i];
        State right = mesh[i + 1];

        fluxes[i] = riemann_solver(left, right);
    }

    return fluxes;
}

Conserved calculate_flux(const State& state) {
    Conserved flux;

    flux.rho = state.rho * state.u;
    flux.momentum = state.rho * state.u * state.u + state.p;
    flux.energy = (state.getEnergy() + state.p) * state.u;

    return flux;
}

double compute_time_step(const std::vector<State>& mesh, double cfl, double dx) {
    int N = mesh.size();
    double max_speed = 0.0;

    for (int i = 0; i < N; i++) {
        double a = std::sqrt(1.4 * mesh[i].p / mesh[i].rho); // Speed of sound
        double speed = std::abs(mesh[i].u) + a;

        if (speed > max_speed) {
            max_speed = speed;
        }
    }

    return cfl * dx / max_speed;
}

std::vector<State> update_states(const std::vector<State>& mesh, const std::vector<Conserved>& fluxes, double dt, double dx) {
    int N = mesh.size();
    std::vector<State> new_mesh(N);

    for (int i = 1; i < N - 1; i++) {
        Conserved conserved = state_to_conserved(mesh[i]);
        
        conserved.rho -= (dt / dx) * (fluxes[i].rho - fluxes[i - 1].rho);
        conserved.momentum -= (dt / dx) * (fluxes[i].momentum - fluxes[i - 1].momentum);
        conserved.energy -= (dt / dx) * (fluxes[i].energy - fluxes[i - 1].energy);

        new_mesh[i] = conserved_to_state(conserved);
    }

    new_mesh[0] = mesh[0];
    new_mesh[N - 1] = mesh[N - 1];

    return new_mesh;
}

void write_snapshot_csv(const std::vector<State>& mesh, double dx, double t, const std::string& filename)
{
    namespace fs = std::filesystem;

    fs::path path(filename);
    if (path.has_parent_path()) {
        fs::create_directories(path.parent_path());
    }

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("write_snapshot_csv: failed to open " + filename);
    }

    out << std::setprecision(17);
    out << "t," << t << "\n";
    out << "i,x,rho,u,p\n";

    for (std::size_t i = 0; i < mesh.size(); ++i) {
        const double x = (static_cast<double>(i) + 0.5) * dx; // cell center
        out << i << "," << x << "," << mesh[i].rho << "," << mesh[i].u << "," << mesh[i].p << "\n";
    }
}