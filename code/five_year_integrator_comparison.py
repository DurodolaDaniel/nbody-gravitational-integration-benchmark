"""Reproduce the common five-year, four-integrator comparison used in Table 2."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import nbody_research as nbody

r0, v0, masses = nbody.sun_earth()
rows = []
for method in ["Euler", "Velocity-Verlet", "RK4", "Yoshida4"]:
    times, positions, velocities = nbody.integrate(
        r0, v0, masses, 5 * nbody.YEAR, 2 * nbody.DAY, method, store_every=1
    )
    analytic = nbody.analytic_two_body(times, masses[0], masses[1])
    position_error_au = np.linalg.norm(
        positions[:, 1] - analytic[:, 1], axis=1
    ) / nbody.AU
    diagnostics = nbody.diagnostics(times, positions, velocities, masses)
    rows.append({
        "method": method,
        "nominal_duration_years": 5.0,
        "actual_duration_years": float(times[-1] / nbody.YEAR),
        "steps": len(times) - 1,
        "dt_days": 2.0,
        "max_relative_energy_error": float(np.max(diagnostics["energy_error"])),
        "max_relative_angular_momentum_error": float(np.max(diagnostics["angular_momentum_error"])),
        "max_earth_position_error_au": float(np.max(position_error_au)),
        "final_earth_position_error_au": float(position_error_au[-1]),
    })

out = ROOT / "results" / "tables" / "five_year_integrator_comparison.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(pd.DataFrame(rows).to_string(index=False))
print(f"Saved {out}")
