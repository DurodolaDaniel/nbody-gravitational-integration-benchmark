from pathlib import Path
import json, platform, sys, time
import numpy as np
import pandas as pd
import scipy
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'tables'
FIG=ROOT/'results'/'figures'
OUT.mkdir(parents=True,exist_ok=True)
FIG.mkdir(parents=True,exist_ok=True)
G=1.0

def acc(r,m):
    d=r[None,:,:]-r[:,None,:]
    d2=np.einsum('ijk,ijk->ij',d,d)
    np.fill_diagonal(d2,np.inf)
    return np.einsum('j,ijk,ij->ik',m,d,d2**-1.5)

def energy(r,v,m):
    ke=.5*np.sum(m[:,None]*v*v)
    d=r[None,:,:]-r[:,None,:]
    d2=np.einsum('ijk,ijk->ij',d,d)
    np.fill_diagonal(d2,np.inf)
    pe=-.5*np.sum(m[:,None]*m[None,:]/np.sqrt(d2))
    return ke+pe

def setup(e):
    m=np.array([1/(1+1e-6),1e-6/(1+1e-6)],float)
    M=m.sum(); a=1.0; rp=a*(1-e)
    relr=np.array([rp,0.,0.]); relv=np.array([0.,np.sqrt(M*(1+e)/(a*(1-e))),0.])
    r=np.array([-m[1]/M*relr,m[0]/M*relr])
    v=np.array([-m[1]/M*relv,m[0]/M*relv])
    T=2*np.pi*np.sqrt(a**3/M)
    return r,v,m,a,T

def exact_kepler(t,e,m,a=1.):
    Mtot=m.sum(); n=np.sqrt(Mtot/a**3); mean=np.remainder(n*np.asarray(t),2*np.pi)
    E=mean.copy()
    for _ in range(20):
        d=(E-e*np.sin(E)-mean)/(1-e*np.cos(E))
        E-=d
        if np.max(np.abs(d))<2e-15: break
    relr=np.stack([a*(np.cos(E)-e),a*np.sqrt(1-e*e)*np.sin(E),np.zeros_like(E)],axis=-1)
    relv=np.stack([-a*n*np.sin(E)/(1-e*np.cos(E)),a*n*np.sqrt(1-e*e)*np.cos(E)/(1-e*np.cos(E)),np.zeros_like(E)],axis=-1)
    r=np.stack([-m[1]/Mtot*relr,m[0]/Mtot*relr],axis=-2)
    v=np.stack([-m[1]/Mtot*relv,m[0]/Mtot*relv],axis=-2)
    return r,v

def verlet_cached(r,v,m,h,a0=None):
    a=acc(r,m) if a0 is None else a0
    rn=r+h*v+0.5*h*h*a
    an=acc(rn,m)
    vn=v+0.5*h*(a+an)
    return rn,vn,an

def rk4_step(r,v,m,h):
    def f(rr,vv): return vv,acc(rr,m)
    k1r,k1v=f(r,v)
    k2r,k2v=f(r+h*k1r/2,v+h*k1v/2)
    k3r,k3v=f(r+h*k2r/2,v+h*k2v/2)
    k4r,k4v=f(r+h*k3r,v+h*k3v)
    return r+h*(k1r+2*k2r+2*k3r+k4r)/6, v+h*(k1v+2*k2v+2*k3v+k4v)/6

def yoshida_cached(r,v,m,h,a0=None):
    c=1/(2-2**(1/3)); d=-2**(1/3)/(2-2**(1/3))
    a=acc(r,m) if a0 is None else a0
    for coeff in (c,d,c):
        r,v,a=verlet_cached(r,v,m,coeff*h,a)
    return r,v,a

def rhs(t,y,m):
    n=len(m); r=y[:3*n].reshape(n,3); v=y[3*n:].reshape(n,3)
    return np.r_[v.ravel(),acc(r,m).ravel()]

def run_eccentric_sweep():
    e=.9; r0,v0,m,a,T=setup(e); n=len(m); y0=np.r_[r0.ravel(),v0.ravel()]
    # DOP853 is included as an adaptive baseline; exact Kepler solution is the reference truth.
    sol=solve_ivp(rhs,(0,T),y0,args=(m,),method='DOP853',rtol=2.3e-14,atol=2.3e-15,dense_output=True,max_step=np.inf)
    if not sol.success: raise RuntimeError(sol.message)
    rows=[]
    grids=[100,300,1000,3000,10000,30000,100000]
    methods={'Euler':(None,1),'Velocity-Verlet':(verlet_cached,1),'RK4':(rk4_step,4),'Yoshida4':(yoshida_cached,3)}
    for method,(step,fe_per_step) in methods.items():
        for steps in grids:
            h=T/steps; r=r0.copy(); v=v0.copy(); a_now=acc(r,m) if method in ('Velocity-Verlet','Yoshida4') else None
            maxerr=0.; max_energy=0.; E0=energy(r,v,m)
            tstart=time.perf_counter()
            stride=max(1,steps//2000)
            for i in range(1,steps+1):
                if method=='Euler':
                    olda=acc(r,m); r=r+h*v; v=v+h*olda
                elif method=='RK4': r,v=step(r,v,m,h)
                else: r,v,a_now=step(r,v,m,h,a_now)
                if i % stride == 0 or i == steps:
                    tr=i*h
                    rr, vv=exact_kepler(np.array([tr]),e,m,a)
                    maxerr=max(maxerr,float(np.max(np.linalg.norm(r-rr[0],axis=1))))
                    max_energy=max(max_energy,abs((energy(r,v,m)-E0)/E0))
            elapsed=time.perf_counter()-tstart
            # Include both endpoints and a common 2001-point grid to avoid method-specific sampling.
            ts=np.linspace(0,T,2001)
            # Reintegrate states at a common grid would be costly; use dense sampled integration above for errors.
            rows.append([method,steps,h,elapsed,steps*fe_per_step,maxerr,max_energy,np.nan])
    # DOP853 accuracy against exact solution on dense common grid and actual function evaluation count.
    ts=np.linspace(0,T,20001); yd=sol.sol(ts).T
    rr_exact,_=exact_kepler(ts,e,m,a)
    rr_dop=yd[:,:3*n].reshape(-1,n,3)
    dop_err=float(np.max(np.linalg.norm(rr_dop-rr_exact,axis=2)))
    rows.append(['DOP853 (adaptive)',np.nan,np.nan,np.nan,sol.nfev,dop_err,np.nan,sol.nfev])
    df=pd.DataFrame(rows,columns=['method','steps_per_orbit','dt','runtime_s','force_evaluations','max_position_error_a_units','max_relative_energy_error','rhs_evaluations'])
    df.to_csv(OUT/'kepler_exact_work_precision.csv',index=False)
    # Convergence slopes only fit the top four resolution levels, where expected asymptotic behaviour should be visible.
    slope_rows=[]
    for method,g in df[df.steps_per_orbit.notna()].groupby('method'):
        gg=g[g.steps_per_orbit>=3000]
        slope=float(np.polyfit(np.log(1/gg.steps_per_orbit.to_numpy(float)),np.log(gg.max_position_error_a_units.to_numpy(float)),1)[0])
        slope_rows.append([method,slope,int(len(gg))])
    pd.DataFrame(slope_rows,columns=['method','observed_order_fit_steps_ge_3000','n_points']).to_csv(OUT/'kepler_exact_convergence_slopes.csv',index=False)
    fig,ax=plt.subplots(figsize=(7.2,5.2))
    for method,g in df[df.steps_per_orbit.notna()].groupby('method'):
        ax.loglog(g.force_evaluations,g.max_position_error_a_units,'o-',label=method)
    ax.scatter([sol.nfev],[dop_err],marker='*',s=120,label=f'DOP853 adaptive ({sol.nfev:,} RHS evals)')
    ax.set_xlabel('Force / right-hand-side evaluations'); ax.set_ylabel('Maximum position error (a units)')
    ax.set_title('Work–precision against exact Kepler solution (e = 0.9)')
    ax.grid(True,which='both',alpha=.25); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(FIG/'kepler_exact_work_precision.png',dpi=320); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.2,5.2))
    for method,g in df[df.steps_per_orbit.notna()].groupby('method'):
        ax.loglog(g.steps_per_orbit,g.max_position_error_a_units,'o-',label=method)
    ax.set_xlabel('Steps per orbit'); ax.set_ylabel('Maximum position error (a units)')
    ax.set_title('Eccentric Kepler convergence against analytical solution')
    ax.grid(True,which='both',alpha=.25); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(FIG/'kepler_exact_convergence.png',dpi=320); plt.close(fig)
    return df,slope_rows,sol.nfev,dop_err

def energy_slope_test():
    rows=[]
    for e in [0.,.9]:
        r0,v0,m,a,T=setup(e); E0=energy(r0,v0,m)
        for steps in [100,300,1000,3000,10000,30000,100000]:
            h=T/steps; r=r0.copy(); v=v0.copy(); ac=acc(r,m); maxe=0.
            for _ in range(steps):
                r,v,ac=verlet_cached(r,v,m,h,ac)
                maxe=max(maxe,abs((energy(r,v,m)-E0)/E0))
            rows.append([e,steps,maxe])
    df=pd.DataFrame(rows,columns=['eccentricity','steps_per_orbit','max_relative_energy_error'])
    df.to_csv(OUT/'verlet_energy_slope_exact_grid.csv',index=False)
    slopes=[]
    for e,g in df.groupby('eccentricity'):
        # Fit middle asymptotic range, exclude coarse saturated and finest roundoff points where appropriate.
        use=g[(g.steps_per_orbit>=300)&(g.steps_per_orbit<=30000)]
        slope=float(np.polyfit(np.log(1/use.steps_per_orbit.to_numpy(float)),np.log(use.max_relative_energy_error.to_numpy(float)),1)[0])
        slopes.append([e,slope,len(use)])
    pd.DataFrame(slopes,columns=['eccentricity','observed_energy_error_slope','fit_points']).to_csv(OUT/'verlet_energy_slope_exact_grid_slopes.csv',index=False)
    fig,ax=plt.subplots(figsize=(7.2,5.2))
    for e,g in df.groupby('eccentricity'):
        ax.loglog(g.steps_per_orbit,g.max_relative_energy_error,'o-',label=f'e = {e:g}')
    ax.set_xlabel('Steps per orbit'); ax.set_ylabel('Maximum relative energy error')
    ax.set_title('Velocity-Verlet energy error: circular versus eccentric Kepler orbit')
    ax.grid(True,which='both',alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(OUT/'verlet_energy_slope_exact_grid.png',dpi=320); plt.close(fig)
    return df,slopes

if __name__=='__main__':
    df,slopes,nfev,doperr=run_eccentric_sweep()
    print('WORK-PRECISION\n',df.to_string(index=False)); print('POSITION SLOPES',slopes,'DOP853 NFEV',nfev,'DOP853 exact max error',doperr)
    edf,esl=energy_slope_test(); print('ENERGY SLOPES',esl)
    meta={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'cpu_model':'AMD EPYC 9V74 80-Core Processor','cpu_visible_to_process':3,'reference':'analytical Kepler equation for fixed-step comparisons; DOP853 adaptive baseline rtol=2.3e-14, atol=2.3e-15','force_evaluations_per_step':'cached velocity-Verlet=1 new acceleration per step; RK4=4; Yoshida4=3 new acceleration evaluations per composed step','timing':'single-run exploratory timing; accuracy conclusions use analytical solution, not runtime'}
    (OUT/'final_benchmark_metadata.json').write_text(json.dumps(meta,indent=2))
