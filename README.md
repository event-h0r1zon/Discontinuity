# Discontinuity
Inviscid &amp; Compressible 1D FVM solver containing implementations of exact and approximate Riemann solvers.

## Requirements
- CMake >= 3.20
- A C++20-capable compiler (e.g., GCC 10+, Clang 12+)
- Python 3.9+ (for post-processing)

## Install Python requirements
From the project root:

```bash
python3 -m pip install -r requirements.txt
```

## Build
This project is built with CMake. Configure and build from a dedicated build folder:

```bash
mkdir -p build
cd build
cmake ..
cmake --build .
```

## Run (must be from /build)
The executable expects to find inputs in ../tests and writes outputs to ../data, so it **must be run from the /build directory**:

```bash
cd build
./discontinuity
```

This will iterate through all test files in tests/ and write step_*.csv snapshots into data/.

## Post-processing utilities (Python)
The script [post_process.py](post_process.py) turns the CSV snapshots into plots and animations.

### Common usages
- Generate a static plot at a target time (closest snapshot chosen automatically) and an animation:

```bash
python3 post_process.py --data-dir data/test_1 --out-dir data/test_1 --t-target 0.25 --anim animation.mp4
```

- Use every Nth frame for faster animations:

```bash
python3 post_process.py --data-dir data/test_1 --out-dir data/test_1 --stride 4 --fps 30 --anim animation.gif
```

- Please make sure the post processing utility script is ran from the projects root.

### Options
- `--data-dir`: Directory containing step_*.csv files (default: data)
- `--pattern`: Glob pattern for snapshots (default: step_*.csv)
- `--out-dir`: Output directory for plots/animations (default: data)
- `--t-target`: Target time for static plot (default: 0.25)
- `--stride`: Use every Nth frame (default: 1)
- `--fps`: Frames per second for animation (default: 30)
- `--anim`: Animation file name (.mp4 or .gif, default: animation.mp4)

## Tests
There are three predefined test cases in tests/ that correspond to the Riemann problems described in *Riemann Solvers and Numerical Methods for Fluid Dynamics* by Eleuterio F. Toro, Chapter 4.

Qualitative comparison of the analytical solution against the numerical one for Test 1 against Fig. 4.7 “Test 1: Exact solution for density, velocity, pressure, and specific internal energy at time t = 0.25” from Toro’s *Riemann Solvers for Fluid Dynamics*:

![Test 1 Comparison - Error Loading](https://i.imgur.com/yS4XTAc.png)

Qualitative comparison of the analytical solution against the numerical one for Test 3 against Fig. 4.9 “Test 3: Exact solution for density, velocity, pressure, and specific internal energy at time t = 0.012” from Toro’s *Riemann Solvers for Fluid Dynamics*:

![Test 3 Comparison - Error Loading](https://i.imgur.com/jh07ayi.png)

### Test file format
Each test is a plain-text file in tests/ (e.g., [tests/test_2.txt](tests/test_2.txt)) with key-value pairs, one per line:

- `name`: Output subfolder name under data/
- `solver`: Riemann solver type to use (see below)
- `N`: Number of cells
- `t_end`: Final simulation time
- `cfl`: CFL number
- `dt_out`: Output interval
- `left`: Left state as `rho u p`
- `right`: Right state as `rho u p`

### Solver types
Two solver types are available:
- `exact`
- `hllc`
