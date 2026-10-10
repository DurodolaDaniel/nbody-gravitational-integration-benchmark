import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

G = 6.67430e-11
AU = 1.495978707e11
DAY = 86400.0
YEAR = 365.25 * DAY
MSUN = 1.98847e30
MEARTH = 3.0034896e-6 * MSUN
MMARS = 3.227151e-7 * MSUN
OUTPUT = "paper_outputs"
os.makedirs(OUTPUT, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 350,
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.22,
    "figure.constrained_layout.use": True
})

COLORS = {
    "Euler": "#D55E00",
    "Verlet": "#0072B2",
    "RK4": "#009E73",
    "Analytic": "#111111",
    "Earth": "#4C78A8",
    "Mars": "#E45756",
    "Sun": "#F2C14E",
    "BinaryA": "#6F4E7C",
    "BinaryB": "#F28E2B",
    "Perturbation": "#59A14F"
}

def acceleration(r, m, eps=0.0):
    d = r[None, :, :] - r[:, None, :]
    d2 = np.einsum("ijk,ijk->ij", d, d)
    np.fill_diagonal(d2, np.inf)
    if eps:
        d2 = d2 + eps * eps
        np.fill_diagonal(d2, np.inf)
    inv_r3 = d2 ** -1.5
    return G * np.einsum("j,ijk,ij->ik", m, d, inv_r3)

def euler_step(r, v, m, dt, eps=0.0):
    a = acceleration(r, m, eps)
    return r + v * dt, v + a * dt

def velocity_verlet_step(r, v, m, dt, eps=0.0):
    """Standalone velocity-Verlet step (two force calls when no acceleration is supplied)."""
    a = acceleration(r, m, eps)
    rn = r + v * dt + 0.5 * a * dt * dt
    an = acceleration(rn, m, eps)
    vn = v + 0.5 * (a + an) * dt
    return rn, vn

def velocity_verlet_cached_step(r, v, m, dt, a, eps=0.0):
    """Velocity-Verlet step that reuses the acceleration from the previous step."""
    rn = r + v * dt + 0.5 * a * dt * dt
    an = acceleration(rn, m, eps)
    vn = v + 0.5 * (a + an) * dt
    return rn, vn, an

def yoshida4_cached_step(r, v, m, dt, a, eps=0.0):
    """Fourth-order Yoshida composition of three cached Verlet substeps."""
    c = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
    d = -(2.0 ** (1.0 / 3.0)) / (2.0 - 2.0 ** (1.0 / 3.0))
    for coeff in (c, d, c):
        r, v, a = velocity_verlet_cached_step(r, v, m, coeff * dt, a, eps)
    return r, v, a

def rk4_step(r, v, m, dt, eps=0.0):
    def f(y):
        rr = y[:3 * len(m)].reshape(len(m), 3)
        vv = y[3 * len(m):].reshape(len(m), 3)
        aa = acceleration(rr, m, eps)
        return np.concatenate((vv.ravel(), aa.ravel()))
    y = np.concatenate((r.ravel(), v.ravel()))
    k1 = f(y)
    k2 = f(y + 0.5 * dt * k1)
    k3 = f(y + 0.5 * dt * k2)
    k4 = f(y + dt * k3)
    yn = y + dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0
    n = len(m)
    return yn[:3*n].reshape(n, 3), yn[3*n:].reshape(n, 3)

INTEGRATORS = {
    "Euler": euler_step,
    "Verlet": velocity_verlet_step,
    "Velocity-Verlet": velocity_verlet_step,
    "RK4": rk4_step,
    "Yoshida4": yoshida4_cached_step,
}

def integrate(r0, v0, m, duration, dt, method="Verlet", eps=0.0, store_every=1):
    stepper = INTEGRATORS[method]
    nsteps = int(np.floor(duration / dt))
    times = np.arange(nsteps + 1, dtype=float) * dt
    stride = max(1, int(store_every))
    indices = np.arange(0, nsteps + 1, stride)
    if indices[-1] != nsteps:
        indices = np.append(indices, nsteps)
    rr = np.empty((len(indices), len(m), 3))
    vv = np.empty_like(rr)
    rr[0] = r0
    vv[0] = v0
    r = r0.copy()
    v = v0.copy()
    # Cache acceleration for the symplectic methods. This makes the force-work
    # accounting consistent with the benchmark: one new evaluation per Verlet
    # macro-step and three per Yoshida4 composition macro-step (after startup).
    a = acceleration(r, m, eps) if method in ("Verlet", "Velocity-Verlet", "Yoshida4") else None
    j = 1
    for k in range(1, nsteps + 1):
        if method in ("Verlet", "Velocity-Verlet"):
            r, v, a = velocity_verlet_cached_step(r, v, m, dt, a, eps)
        elif method == "Yoshida4":
            r, v, a = yoshida4_cached_step(r, v, m, dt, a, eps)
        else:
            r, v = stepper(r, v, m, dt, eps)
        if j < len(indices) and k == indices[j]:
            rr[j] = r
            vv[j] = v
            j += 1
    return times[indices], rr, vv

def total_energy(r, v, m, eps=0.0):
    kinetic = 0.5 * np.sum(m[None, :] * np.sum(v * v, axis=2), axis=1)
    d = r[:, :, None, :] - r[:, None, :, :]
    d2 = np.sum(d * d, axis=3)
    iu = np.triu_indices(len(m), 1)
    dist = np.sqrt(d2[:, iu[0], iu[1]] + eps * eps)
    potential = -G * np.sum((m[iu[0]] * m[iu[1]])[None, :] / dist, axis=1)
    return kinetic + potential

def angular_momentum(r, v, m):
    return np.sum(np.cross(r, m[None, :, None] * v), axis=1)

def linear_momentum(v, m):
    return np.sum(m[None, :, None] * v, axis=1)

def center_of_mass(r, m):
    return np.sum(r * m[None, :, None], axis=1) / np.sum(m)

def diagnostics(times, r, v, m, eps=0.0):
    E = total_energy(r, v, m, eps)
    L = angular_momentum(r, v, m)
    P = linear_momentum(v, m)
    Rcm = center_of_mass(r, m)
    Eerr = np.abs((E - E[0]) / E[0])
    L0 = np.linalg.norm(L[0])
    Lerr = np.linalg.norm(L - L[0], axis=1) / L0 if L0 else np.zeros(len(times))
    Pscale = np.sum(m[None, :] * np.linalg.norm(v, axis=2), axis=1)
    Perr = np.linalg.norm(P - P[0], axis=1) / max(Pscale[0], np.finfo(float).tiny)
    COMerr = np.linalg.norm(Rcm - Rcm[0], axis=1) / AU
    return {"energy": E, "angular_momentum": L, "momentum": P, "com": Rcm,
            "energy_error": Eerr, "angular_momentum_error": Lerr,
            "momentum_error": Perr, "com_drift_au": COMerr}

def circular_two_body(m1, m2, separation_au=1.0):
    a = separation_au * AU
    M = m1 + m2
    r1 = a * m2 / M
    r2 = a * m1 / M
    omega = np.sqrt(G * M / a**3)
    r = np.array([[-r1, 0, 0], [r2, 0, 0]], float)
    v = np.array([[0, -omega*r1, 0], [0, omega*r2, 0]], float)
    return r, v, np.array([m1, m2], float)

def sun_earth():
    return circular_two_body(MSUN, MEARTH, 1.0)

def sun_earth_mars():
    m = np.array([MSUN, MEARTH, MMARS])
    r = np.array([[0,0,0], [AU,0,0], [1.524*AU,0,0]], float)
    vE = np.sqrt(G * np.sum(m[:2]) / AU)
    vM = np.sqrt(G * np.sum(m[[0,2]]) / (1.524*AU))
    v = np.array([[0,0,0], [0,vE,0], [0,vM,0]], float)
    v[0] = -(m[1]*v[1] + m[2]*v[2]) / np.sum(m)
    rcm = np.sum(r*m[:,None], axis=0) / np.sum(m)
    vcm = np.sum(v*m[:,None], axis=0) / np.sum(m)
    return r-rcm, v-vcm, m

def binary_stars():
    return circular_two_body(MSUN, 0.8*MSUN, 1.0)

def analytic_two_body(times, m1, m2, separation_au=1.0):
    a = separation_au * AU
    M = m1 + m2
    r1 = a * m2 / M
    r2 = a * m1 / M
    omega = np.sqrt(G*M/a**3)
    theta = omega * times
    r = np.zeros((len(times), 2, 3))
    r[:,0,0] = -r1*np.cos(theta)
    r[:,0,1] = -r1*np.sin(theta)
    r[:,1,0] = r2*np.cos(theta)
    r[:,1,1] = r2*np.sin(theta)
    return r

def plot_save(name, title):
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    return fig, ax

def save(fig, name):
    fig.savefig(os.path.join(OUTPUT, name), bbox_inches="tight")
    plt.close(fig)

def experiment_two_body_orbit():
    r0,v0,m = sun_earth()
    t,r,v = integrate(r0,v0,m,YEAR,DAY/2,"Verlet")
    ra = analytic_two_body(t,m[0],m[1])
    fig,ax = plot_save("", "")
    ax.plot(ra[:,1,0]/AU, ra[:,1,1]/AU, color=COLORS["Analytic"], lw=2.2, ls="--", label="Analytical Earth orbit")
    ax.plot(r[:,1,0]/AU, r[:,1,1]/AU, color=COLORS["Earth"], lw=1.8, label="Velocity-Verlet")
    ax.scatter(r[0,1,0]/AU,r[0,1,1]/AU,s=80,color=COLORS["Earth"],edgecolor="white",zorder=5,label="Initial state")
    ax.set_xlabel("x position (AU)")
    ax.set_ylabel("y position (AU)")
    ax.set_title("Sun–Earth Orbit: Numerical and Analytical Solutions")
    ax.set_aspect("equal")
    ax.legend(frameon=False)
    save(fig,"fig01_sun_earth_validation.png")
    return t,r,v,m

def experiment_integrator_errors():
    r0,v0,m = sun_earth()
    methods = ["Euler","Verlet","RK4"]
    out = {}
    for method in methods:
        t,r,v = integrate(r0,v0,m,5*YEAR,2*DAY,method,store_every=2)
        out[method] = (t, diagnostics(t,r,v,m))
    fig,ax = plot_save("","")
    for method in methods:
        t,d = out[method]
        ax.semilogy(t/YEAR,d["energy_error"],lw=2,label=method,color=COLORS[method])
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Relative energy error")
    ax.set_title("Long-Term Energy Conservation Across Integrators")
    ax.legend(frameon=False)
    save(fig,"fig02_energy_error_integrators.png")
    fig,ax = plot_save("","")
    for method in methods:
        t,d = out[method]
        ax.semilogy(t/YEAR,np.maximum(d["angular_momentum_error"],1e-18),lw=2,label=method,color=COLORS[method])
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Relative angular-momentum error")
    ax.set_title("Long-Term Angular-Momentum Conservation Across Integrators")
    ax.legend(frameon=False)
    save(fig,"fig03_angular_momentum_error_integrators.png")
    return out

def experiment_convergence():
    r0,v0,m = sun_earth()
    dts = np.array([4,2,1,0.5,0.25,0.125])*DAY
    methods = ["Euler","Verlet","RK4"]
    rows=[]
    for method in methods:
        for dtday in dts/DAY:
            t,r,v = integrate(r0,v0,m,YEAR,dtday*DAY,method)
            ra = analytic_two_body(t,m[0],m[1])
            poserr = np.linalg.norm(r[:,1]-ra[:,1],axis=1)/AU
            d = diagnostics(t,r,v,m)
            rows.append([method,dtday,np.max(poserr),np.max(d["energy_error"])])
    df=pd.DataFrame(rows,columns=["method","dt_days","max_position_error_au","max_energy_error"])
    df.to_csv(os.path.join(OUTPUT,"convergence_results.csv"),index=False)
    fig,ax=plot_save("","")
    for method in methods:
        q=df[df.method==method]
        ax.loglog(q.dt_days,q.max_position_error_au,"o-",lw=2,label=method,color=COLORS[method])
    ax.set_xlabel("Timestep (days)")
    ax.set_ylabel("Maximum Earth position error (AU)")
    ax.set_title("Timestep Convergence of the Numerical Solutions")
    ax.invert_xaxis()
    ax.legend(frameon=False)
    save(fig,"fig04_timestep_convergence_position.png")
    fig,ax=plot_save("","")
    for method in methods:
        q=df[df.method==method]
        ax.loglog(q.dt_days,q.max_energy_error,"o-",lw=2,label=method,color=COLORS[method])
    ax.set_xlabel("Timestep (days)")
    ax.set_ylabel("Maximum relative energy error")
    ax.set_title("Timestep Dependence of Energy Conservation")
    ax.invert_xaxis()
    ax.legend(frameon=False)
    save(fig,"fig05_timestep_convergence_energy.png")
    return df

def experiment_orbital_error():
    r0,v0,m = sun_earth()
    t,r,v = integrate(r0,v0,m,10*YEAR,0.5*DAY,"Verlet")
    ra=analytic_two_body(t,m[0],m[1])
    err=np.linalg.norm(r[:,1]-ra[:,1],axis=1)/AU
    radial=np.linalg.norm(r[:,1]-r[:,0],axis=1)/AU
    fig,ax=plot_save("","")
    ax.plot(t/YEAR,err,lw=2,color=COLORS["Perturbation"])
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Earth position error (AU)")
    ax.set_title("Long-Term Phase and Position Error of Velocity-Verlet")
    save(fig,"fig06_long_term_orbital_error.png")
    fig,ax=plot_save("","")
    ax.plot(t/YEAR,radial,lw=1.8,color=COLORS["Earth"])
    ax.axhline(1.0,color=COLORS["Analytic"],ls="--",lw=1.4,label="Reference separation")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Earth–Sun separation (AU)")
    ax.set_title("Long-Term Orbital Radius Stability")
    ax.legend(frameon=False)
    save(fig,"fig07_orbital_radius_stability.png")

def experiment_mars_perturbation():
    r2,v2,m2=sun_earth()
    t2,r2h,v2h=integrate(r2,v2,m2,5*YEAR,1*DAY,"Verlet")
    r3,v3,m3=sun_earth_mars()
    t3,r3h,v3h=integrate(r3,v3,m3,5*YEAR,1*DAY,"Verlet")
    e2=r2h[:,1]
    e3=r3h[:,1]
    delta=np.linalg.norm(e3-e2,axis=1)/AU
    fig,ax=plot_save("","")
    ax.plot(t3/YEAR,delta,lw=2,color=COLORS["Perturbation"])
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Earth trajectory separation (AU)")
    ax.set_title("Earth Orbital Perturbation Introduced by Mars")
    save(fig,"fig08_mars_earth_perturbation.png")
    fig,ax=plot_save("","")
    ax.plot(e2[:,0]/AU,e2[:,1]/AU,lw=1.7,color=COLORS["Earth"],label="Sun–Earth")
    ax.plot(e3[:,0]/AU,e3[:,1]/AU,lw=1.7,color=COLORS["Perturbation"],label="Sun–Earth–Mars")
    ax.set_xlabel("x position (AU)")
    ax.set_ylabel("y position (AU)")
    ax.set_title("Earth Orbit With and Without the Martian Perturber")
    ax.set_aspect("equal")
    ax.legend(frameon=False)
    save(fig,"fig09_earth_orbit_mars_comparison.png")

def experiment_binary():
    r0,v0,m=binary_stars()
    t,r,v=integrate(r0,v0,m,5*YEAR,0.5*DAY,"Verlet")
    ra=analytic_two_body(t,m[0],m[1],1.0)
    fig,ax=plot_save("","")
    ax.plot(ra[:,0,0]/AU,ra[:,0,1]/AU,ls="--",lw=2,color=COLORS["Analytic"],label="Analytical Star A")
    ax.plot(ra[:,1,0]/AU,ra[:,1,1]/AU,ls="--",lw=2,color=COLORS["Analytic"],label="Analytical Star B")
    ax.plot(r[:,0,0]/AU,r[:,0,1]/AU,lw=1.8,color=COLORS["BinaryA"],label="Numerical Star A")
    ax.plot(r[:,1,0]/AU,r[:,1,1]/AU,lw=1.8,color=COLORS["BinaryB"],label="Numerical Star B")
    ax.set_xlabel("x position (AU)")
    ax.set_ylabel("y position (AU)")
    ax.set_title("Equal-Plane Binary-Star Validation")
    ax.set_aspect("equal")
    ax.legend(frameon=False,ncol=2)
    save(fig,"fig10_binary_validation.png")
    return t,r,v,m

def experiment_momentum_com():
    r0,v0,m=sun_earth_mars()
    t,r,v=integrate(r0,v0,m,5*YEAR,1*DAY,"Verlet")
    d=diagnostics(t,r,v,m)
    fig,ax=plot_save("","")
    ax.semilogy(t/YEAR,np.maximum(d["com_drift_au"],1e-18),lw=2,color=COLORS["Earth"])
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Centre-of-mass displacement (AU)")
    ax.set_title("Centre-of-Mass Stability in the Three-Body System")
    save(fig,"fig11_center_of_mass_drift.png")
    fig,ax=plot_save("","")
    ax.semilogy(t/YEAR,np.maximum(d["momentum_error"],1e-18),lw=2,color=COLORS["Perturbation"])
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Relative linear-momentum error")
    ax.set_title("Linear-Momentum Conservation")
    save(fig,"fig12_linear_momentum_error.png")

def make_cluster(n, seed=42):
    rng=np.random.default_rng(seed)
    m=rng.uniform(0.1,2.0,n)*MSUN
    r=rng.normal(0,0.6*AU,(n,3))
    v=rng.normal(0,5000,(n,3))
    r-=np.sum(r*m[:,None],axis=0)/np.sum(m)
    v-=np.sum(v*m[:,None],axis=0)/np.sum(m)
    return r,v,m

def experiment_scaling():
    ns=np.array([10,25,50,100,200,400])
    rows=[]
    for n in ns:
        r,v,m=make_cluster(n)
        dt=0.25*DAY
        t0=time.perf_counter()
        repeats=3
        elapsed_values=[]
        for _ in range(repeats):
            t0=time.perf_counter()
            integrate(r,v,m,10*dt,dt,"Verlet")
            elapsed_values.append(time.perf_counter()-t0)
        elapsed=float(np.median(elapsed_values))
        rows.append([n,elapsed,elapsed/10])
    df=pd.DataFrame(rows,columns=["N","runtime_seconds","seconds_per_step"])
    df.to_csv(os.path.join(OUTPUT,"scaling_results.csv"),index=False)
    fig,ax=plot_save("","")
    ax.loglog(df.N,df.runtime_seconds,"o-",lw=2,color=COLORS["BinaryA"])
    slope=np.polyfit(np.log(df.N),np.log(df.runtime_seconds),1)[0]
    ax.set_xlabel("Number of bodies, N")
    ax.set_ylabel("Runtime for two timesteps (s)")
    ax.set_title(f"Direct N-Body Computational Scaling (measured slope = {slope:.2f})")
    save(fig,"fig13_computational_scaling.png")
    return df

def experiment_softening():
    r0,v0,m=binary_stars()
    eps=np.array([0,1e-4,5e-4,1e-3,5e-3])*AU
    rows=[]
    for e in eps:
        t,r,v=integrate(r0,v0,m,2*YEAR,0.25*DAY,"Verlet",e)
        d=diagnostics(t,r,v,m,e)
        rows.append([e/AU,np.max(d["energy_error"])])
    df=pd.DataFrame(rows,columns=["softening_au","max_energy_error"])
    df.to_csv(os.path.join(OUTPUT,"softening_results.csv"),index=False)
    fig,ax=plot_save("","")
    ax.loglog(np.maximum(df.softening_au,1e-8),df.max_energy_error,"o-",lw=2,color=COLORS["Mars"])
    ax.set_xlabel("Softening length (AU)")
    ax.set_ylabel("Maximum relative energy error")
    ax.set_title("Sensitivity of Numerical Conservation to Gravitational Softening")
    save(fig,"fig14_softening_sensitivity.png")
    return df

def summary():
    r0,v0,m=sun_earth()
    t,r,v=integrate(r0,v0,m,YEAR,0.5*DAY,"Verlet")
    d=diagnostics(t,r,v,m)
    ra=analytic_two_body(t,m[0],m[1])
    poserr=np.linalg.norm(r[:,1]-ra[:,1],axis=1)/AU
    result=pd.DataFrame([{
        "system":"Sun-Earth",
        "integrator":"Velocity-Verlet",
        "dt_days":0.5,
        "duration_years":1.0,
        "max_relative_energy_error":np.max(d["energy_error"]),
        "max_relative_angular_momentum_error":np.max(d["angular_momentum_error"]),
        "max_position_error_AU":np.max(poserr),
        "max_COM_drift_AU":np.max(d["com_drift_au"])
    }])
    result.to_csv(os.path.join(OUTPUT,"summary_results.csv"),index=False)
    return result

def run_all():
    experiment_two_body_orbit()
    experiment_integrator_errors()
    experiment_convergence()
    experiment_orbital_error()
    experiment_mars_perturbation()
    experiment_binary()
    experiment_momentum_com()
    experiment_scaling()
    experiment_softening()
    result=summary()
    print(result.to_string(index=False))
    print(f"\nFigures and tables saved to: {OUTPUT}")

if __name__=="__main__":
    run_all()
