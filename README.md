## Direct Gravitational N-Body Integration Benchmark

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23566896.svg)](https://doi.org/10.5281/zenodo.23566896)

[![Release](https://img.shields.io/badge/release-v1.2.0-blue)](https://github.com/DurodolaDaniel/nbody-gravitational-integration-benchmark/releases/tag/v1.2.0)

A reproducible computational-physics benchmark comparing numerical methods for direct Newtonian gravitational N-body integration. The project studies the distinction between conservation behaviour and trajectory accuracy using analytical Kepler solutions, timestep refinement, long-duration eccentric-orbit integrations, and multi-body examples.

Latest release: v1.2.0
Latest release DOI: https://doi.org/10.5281/zenodo.23566896

The earlier releases remain available:

v1.1.0: https://doi.org/10.5281/zenodo.23480336

v1.0.0: https://doi.org/10.5281/zenodo.23286569

Each version-specific DOI identifies a particular archived software release. Changes made to the GitHub repository after a release are not automatically included in that archived version.

GitHub repository: https://github.com/DurodolaDaniel/nbody-gravitational-integration-benchmark

> The DOI identifies the archived software release. It is not a journal-publication DOI for the manuscript.

## Numerical methods

The benchmark code includes:

• Forward Euler
• Velocity-Verlet
• Classical fourth-order Runge-Kutta (RK4)
• Fourth-order Yoshida composition (Yoshida4)
• SciPy’s adaptive DOP853 integrator as a numerical baseline in the eccentric Kepler work–precision experiment

The direct-force model uses Newtonian pairwise gravitational interactions. Diagnostics include relative energy error, angular-momentum error, position error against an analytical reference, centre-of-mass drift, linear-momentum error, convergence behaviour, and computational work.

## Repository contents

• code/nbody_research.py — core direct-force model, initial conditions, integrators, diagnostics, and the main set of Sun–Earth, Sun–Earth–Mars, binary-star, scaling, and softening experiments.

• code/five_year_integrator_comparison.py — regenerates the four-method, five-year Sun–Earth comparison table.

• code/make_main_yoshida_figures.py — regenerates the four-method comparison and timestep-convergence figures, including Yoshida4.

• code/final_kepler_benchmark.py — compares fixed-step methods with an analytical eccentric Kepler solution; also generates work–precision, convergence, and velocity-Verlet energy-error analyses.

• code/euler_fine_convergence.py — performs a fine-grid Euler convergence check.

• code/long_duration_and_figure8.py — runs the 1,000-orbit eccentric Kepler resolution sweep and the figure-eight return-error experiment. This is a longer run and uses Numba-compiled kernels.

• results/tables/ — CSV benchmark tables and run metadata.

• results/figures/ — plots used for analysis and the manuscript.

• manuscript/ — PDF and editable Word manuscript files.

• CITATION.cff — software citation metadata.

• requirements.txt — pinned Python dependencies.

## Requirements

The archived experiments were prepared with Python 3.13 and these pinned dependencies:

• NumPy 2.3.5
• SciPy 1.17.0
• Matplotlib 3.10.8
• pandas 2.2.3
• Numba 0.65.1

Other compatible environments may work, but the exact versions above are the reference environment recorded for this release. Runtime measurements depend on hardware and software environment.

## Installation

Run these commands from the repository root.

Windows (PowerShell or VS Code terminal)

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If PowerShell blocks activation, you can open a Command Prompt terminal and use:

```bat
.venv\Scripts\activate.bat
```

macOS or Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Reproducing the outputs

Open a terminal in the top-level repository directory—the directory containing README.md, requirements.txt, code/, and results/. Run scripts from that directory so relative output paths resolve as expected.

1. Quick first check: five-year integrator comparison

```bash
python code/five_year_integrator_comparison.py
```

This prints a table comparing Euler, Velocity-Verlet, RK4, and Yoshida4 over the same nominal five-year Sun–Earth experiment and writes:

• results/tables/five_year_integrator_comparison.csv

Check that the CSV contains four rows, one for each method, and columns for energy error, angular-momentum error, and Earth position error. The simulated duration is 4.999316 years because 913 complete steps of two days do not equal exactly five Julian years.

2. Regenerate the four-method comparison and convergence figures

```bash
python code/make_main_yoshida_figures.py
```

Writes or updates:

• results/figures/fig02_energy_error_integrators.png
• results/figures/fig03_angular_momentum_error_integrators.png
• results/figures/fig04_timestep_convergence_position.png
• results/figures/fig05_timestep_convergence_energy.png
• results/tables/convergence_results.csv

Open the PNGs and check the axis labels, legends, log scales, and whether the methods separate as expected. This script may display non-fatal Matplotlib layout warnings on some versions.

3. Analytical eccentric-Kepler benchmark and work–precision results

```bash
python code/final_kepler_benchmark.py
```

This compares fixed-step methods against an analytical Kepler solution for eccentricity 0.9 and runs a separate velocity-Verlet energy-error sweep. It writes or updates:

• results/tables/kepler_exact_work_precision.csv
• results/tables/kepler_exact_convergence_slopes.csv
• results/tables/verlet_energy_slope_exact_grid.csv
• results/tables/verlet_energy_slope_exact_grid_slopes.csv
• results/tables/final_benchmark_metadata.json
• results/figures/kepler_exact_work_precision.png
• results/figures/kepler_exact_convergence.png
• results/figures/verlet_energy_slope_exact_grid.png

Inspect the CSVs as well as the plots. Position error is measured against the analytical Kepler solution; DOP853 is an adaptive numerical baseline, not the exact reference. Runtime values can vary between machines.

4. Fine-grid Euler convergence

```bash
python code/euler_fine_convergence.py
```

Writes or updates:

• results/tables/euler_fine_convergence.csv
• results/figures/euler_fine_convergence.png

This run includes up to 1,000,000 steps for a one-year circular two-body orbit. It prints the tested step counts, timestep sizes, maximum position errors, and fitted slope. Check that the error decreases as the number of steps increases.

5. Long-duration eccentric orbit and figure-eight return experiment

```bash
python code/long_duration_and_figure8.py
```

This is the longer run. It tests 1,000 eccentric orbits at 5,000, 10,000, and 20,000 steps per orbit for Velocity-Verlet, RK4, and Yoshida4. It also performs a figure-eight one-period return test.

It writes or updates:

• results/tables/eccentric_1000_orbits_resolution_comparison.csv
• Dense sampled diagnostic CSV files for each method and resolution
• results/tables/figure8_return_error.csv
• results/tables/longrun_environment.txt
• results/figures/eccentric_long_duration.png

The first run may take longer because Numba compiles the numerical kernels. Check the terminal output for all nine method/resolution combinations and inspect the resolution-comparison CSV. The figure-eight experiment uses rounded initial conditions and is not treated as an external validation or as a basis for ranking the integrators.

6. Run the original/main experiment suite

```bash
python code/nbody_research.py
```

This runs the core experiments: Sun–Earth orbit comparison, long-term energy and angular-momentum diagnostics, timestep convergence, long-term orbital error, Mars perturbation, binary stars, momentum and centre-of-mass diagnostics, computational scaling, gravitational softening sensitivity, and a summary calculation.

When launched from the repository root, it creates a paper_outputs/ directory in the root and writes the figures and CSV tables there, including summary_results.csv, convergence_results.csv, scaling_results.csv, and softening_results.csv. This output location is separate from the curated results/ directory.

## How to check that a run worked

1. Run one script at a time and wait for the terminal prompt to return.
2. Read the printed summary and check for a traceback or error message.
3. Confirm the expected CSV/PNG files exist in the paths listed above.
4. Open CSV files and inspect their column names, row counts, and numerical values.
5. Open the PNG files and check that the plots have axes, labels, legends, and finite plotted values.
6. Compare reproduced numerical results with the bundled reference CSVs. Small differences in runtime are normal; numerical differences should be investigated rather than silently accepted.
7. Before committing regenerated outputs, run git status and review the changes. Timing metadata and generated files can change when scripts are rerun.

These scripts are reproducibility/benchmark programs; the repository does not currently include a dedicated automated unit-test suite. For a stronger publication package, add tests for the force calculation, conservation/centre-of-mass properties, and convergence against analytical solutions, and run them in a clean environment.

## Interpretation and limitations

• Good energy conservation does not guarantee small trajectory error over very long integrations; phase error must also be examined.

• In the 1,000-orbit eccentric test, Velocity-Verlet’s trajectory discrepancy remains large at the tested resolutions, while RK4 and Yoshida4 show smaller trajectory errors. Do not rank methods from energy error alone.

• The figure-eight experiment uses rounded initial conditions, so its return-error floor is influenced by the precision of those initial conditions.

• The Pythagorean close-encounter experiment is under-resolved and is not used to rank methods.

• Comparisons against REBOUND and JPL Horizons have not been completed. The project therefore does not claim external validation against those tools.

• Benchmark timings are hardware- and environment-dependent.

## Citation

If you use this software or its benchmark outputs, cite the specific archived release:

Durodola, D. (2026). Conservation versus Trajectory Accuracy in Direct Gravitational N-Body Integration: A Benchmark of Euler, Velocity-Verlet, RK4, and Yoshida4 (Version v1.2.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23566896

The DOI above identifies the archived v1.2.0 software release. It does not imply that the associated manuscript has been peer-reviewed or published. Cite the manuscript separately according to its publication status.

## License

The Zenodo record currently identifies the release as MIT License. Before reusing or redistributing the code, confirm that a matching LICENSE file is present in the GitHub repository and that the license metadata is consistent.
