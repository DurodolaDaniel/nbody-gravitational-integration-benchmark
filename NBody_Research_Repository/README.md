# Direct Gravitational N-Body Integration Benchmark

This repository accompanies **“Conservation versus Trajectory Accuracy in Direct Gravitational N-Body Integration: A Benchmark of Euler, Velocity-Verlet, RK4, and Yoshida4.”** It is a numerical-methods benchmark, not a new integration algorithm.

## Contents

- `manuscript/Durodola_NBody_Benchmark.pdf` — corrected manuscript PDF.
- `manuscript/Durodola_NBody_Benchmark.docx` — editable manuscript source.
- `code/nbody_research.py` — direct-force N-body model, integrators, diagnostics, and original experiments.
- `code/final_kepler_benchmark.py` — fixed-step work–precision and convergence tests against the analytical two-body Kepler solution; DOP853 is an adaptive numerical baseline.
- `code/euler_fine_convergence.py` — fine-grid Euler convergence check.
- `code/make_main_yoshida_figures.py` — main four-integrator comparison and convergence figures.
- `code/five_year_integrator_comparison.py` — regenerates the common five-year, four-integrator comparison in Table 2.
- `code/long_duration_and_figure8.py` — reruns the 1,000-orbit eccentric sweep at 5,000, 10,000, and 20,000 steps/orbit, with 200 diagnostic samples per orbit, and reruns the figure-eight return test.
- `results/tables/` — CSV tables, including the five-year comparison and long-duration resolution sweep.
- `results/figures/` — figures used in the manuscript and numerical diagnostics. The long-duration eccentric Figure 8 has been regenerated and visually checked.
- `CITATION.cff` — citation metadata with the author's ORCID.

## Environment

The experiments were run with Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, Numba 0.65.1, and pandas 2.2.3 on Linux x86_64. The runtime reported an AMD EPYC 9V74 80-Core Processor, with three logical CPUs available to the process. Timing results are environment-dependent and should not be interpreted as hardware-independent performance claims.

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

Reproduce the main experiments from the repository root:

```bash
python code/final_kepler_benchmark.py
python code/euler_fine_convergence.py
python code/five_year_integrator_comparison.py
python code/long_duration_and_figure8.py
python code/make_main_yoshida_figures.py
```

The long-duration script may take time to compile and run its numerical kernels. It writes per-resolution sampled diagnostic CSVs, the resolution-comparison table, the figure-eight return table, and the long-duration figure. Full sampled traces are generated on demand rather than bundled, to keep the source package compact; summary tables and the principal plot are included.

## Numerical conventions

The eccentric fixed-step benchmark uses the analytical Kepler equation as the reference. Long-duration diagnostics are sampled 200 times per orbit to resolve the high angular speed near periapsis and avoid phase-unwrapping ambiguity. Figure 8 displays one sample per orbit for readability; its resolution-sweep table reports maxima from the denser sampling. DOP853 is an adaptive numerical baseline, not exact truth.

Force-evaluation accounting uses cached acceleration: one new acceleration evaluation per velocity-Verlet step, four per RK4 step, and three per Yoshida4 composition step. The adaptive DOP853 cost is reported using SciPy's solver-reported right-hand-side evaluation count.

The Pythagorean close-encounter experiment remains under-resolved and is not used to rank methods. Comparisons with REBOUND and JPL Horizons have not been completed. No claim of external validation is made.

## Reference audit

The Harfst entry was corrected to the article *Performance analysis of direct N-body algorithms on special-purpose supercomputers*, DOI `10.1016/j.newast.2006.11.003`. The previous DOI `10.1016/j.parco.2007.01.001` belongs to a different article by Gualandris and colleagues, not the Harfst paper. The other DOI-bearing references were cross-checked against publisher or bibliographic records; no additional DOI mismatch was identified in this review. The IAU Resolution B2 entry is an institutional URL, not a DOI.

## Public archive status — action required before posting

This prepared package has **not** been uploaded to a public GitHub repository and has **no Zenodo DOI**. A DOI is intentionally not supplied or fabricated. The code and manuscript must be uploaded to an authenticated public GitHub repository, then a tagged release must be archived through Zenodo. After Zenodo issues the DOI, update the DOI and public repository URL in this README, `CITATION.cff`, and the manuscript before restoring any claim that a permanent public archive is available. Until then, the manuscript does not claim a public archive.
