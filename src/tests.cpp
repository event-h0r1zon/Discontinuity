#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include "riemann.h"
#include "solver.h"
#include <vector>

struct TestCase {
    std::string name;
    std::string solver;
    int N = 100;
    double t_end = 0.25;
    double cfl = 0.5;
    double dt_out = 0.005;
    State left{}, right{};
};

static std::string trim(const std::string& str) {
    const char* whitespace = " \t\n\r";
    size_t start = str.find_first_not_of(whitespace);
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(whitespace);
    return str.substr(start, end - start + 1);
}

static TestCase parse_test_file(const std::filesystem::path& path) {
    TestCase tc;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open test file: " + path.string());

    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto eq = line.find('=');
        if (eq == std::string::npos) continue;

        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));

        std::istringstream iss(val);
        if (key == "name") tc.name = val;
        else if (key == "solver") iss >> tc.solver;
        else if (key == "N") iss >> tc.N;
        else if (key == "t_end") iss >> tc.t_end;
        else if (key == "cfl") iss >> tc.cfl;
        else if (key == "dt_out") iss >> tc.dt_out;
        else if (key == "left") iss >> tc.left.rho >> tc.left.u >> tc.left.p;
        else if (key == "right") iss >> tc.right.rho >> tc.right.u >> tc.right.p;
    }

    if (tc.name.empty()) tc.name = path.stem().string();
    return tc;
}

static void run_test(const TestCase& tc, const std::filesystem::path& out_base_dir) {
    namespace fs = std::filesystem;
    fs::path out_dir = out_base_dir / tc.name;
    fs::create_directories(out_dir);

    const double dx = 1.0 / tc.N;
    double t = 0.0;
    int step = 0;

    std::vector<State> mesh = initialize_mesh(tc.N, tc.left, tc.right);

    write_snapshot_csv(mesh, dx, t, (out_dir / "step_000000.csv").string());

    double next_out_t = tc.dt_out;

    while (t < tc.t_end) {
        auto fluxes = calculate_fluxes(mesh, tc.solver);
        double dt = compute_time_step(mesh, tc.cfl, dx);
        if (t + dt > tc.t_end) dt = tc.t_end - t;

        mesh = update_states(mesh, fluxes, dt, dx);
        t += dt;
        ++step;

        if (t + 1e-12 >= next_out_t || t + 1e-12 >= tc.t_end) {
            std::ostringstream name;
            name << "step_" << std::setw(6) << std::setfill('0') << step << ".csv";
            write_snapshot_csv(mesh, dx, t, (out_dir / name.str()).string());
            next_out_t += tc.dt_out;
        }
    }
}