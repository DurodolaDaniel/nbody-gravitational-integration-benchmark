import numpy as np, pandas as pd, matplotlib.pyplot as plt, time, platform, os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'tables'
FIG = ROOT / 'results' / 'figures'
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
from numba import njit

G=1.0; PERIOD=2*np.pi; E=0.9; A=1.0; MU=1.0
@njit
def acc(x,y):
    r2=x*x+y*y
    inv=r2**(-1.5)
    return -x*inv,-y*inv
@njit
def verlet_step(x,y,vx,vy,h,ax,ay):
    xn=x+h*vx+0.5*h*h*ax
    yn=y+h*vy+0.5*h*h*ay
    axn,ayn=acc(xn,yn)
    vxn=vx+0.5*h*(ax+axn)
    vyn=vy+0.5*h*(ay+ayn)
    return xn,yn,vxn,vyn,axn,ayn
@njit
def run_verlet(nsteps, sample_stride, h):
    x,y,vx,vy=0.1,0.0,0.0,np.sqrt(19.0)
    ax,ay=acc(x,y)
    nout=nsteps//sample_stride+1
    out=np.zeros((nout,6)); j=0
    out[j]=[0,x,y,vx,vy,0.5*(vx*vx+vy*vy)-1/np.sqrt(x*x+y*y)]; j+=1
    for i in range(1,nsteps+1):
        x,y,vx,vy,ax,ay=verlet_step(x,y,vx,vy,h,ax,ay)
        if i%sample_stride==0:
            out[j]=[i,x,y,vx,vy,0.5*(vx*vx+vy*vy)-1/np.sqrt(x*x+y*y)]; j+=1
    return out[:j]
@njit
def run_rk4(nsteps, sample_stride, h):
    x,y,vx,vy=0.1,0.0,0.0,np.sqrt(19.0)
    nout=nsteps//sample_stride+1
    out=np.zeros((nout,6)); j=0
    out[j]=[0,x,y,vx,vy,0.5*(vx*vx+vy*vy)-1/np.sqrt(x*x+y*y)]; j+=1
    for i in range(1,nsteps+1):
        ax1,ay1=acc(x,y); k1x,k1y=vx,vy; k1vx,k1vy=ax1,ay1
        ax2,ay2=acc(x+0.5*h*k1x,y+0.5*h*k1y); k2x,k2y=vx+0.5*h*k1vx,vy+0.5*h*k1vy; k2vx,k2vy=ax2,ay2
        ax3,ay3=acc(x+0.5*h*k2x,y+0.5*h*k2y); k3x,k3y=vx+0.5*h*k2vx,vy+0.5*h*k2vy; k3vx,k3vy=ax3,ay3
        ax4,ay4=acc(x+h*k3x,y+h*k3y); k4x,k4y=vx+h*k3vx,vy+h*k3vy; k4vx,k4vy=ax4,ay4
        x += h*(k1x+2*k2x+2*k3x+k4x)/6
        y += h*(k1y+2*k2y+2*k3y+k4y)/6
        vx += h*(k1vx+2*k2vx+2*k3vx+k4vx)/6
        vy += h*(k1vy+2*k2vy+2*k3vy+k4vy)/6
        if i%sample_stride==0:
            out[j]=[i,x,y,vx,vy,0.5*(vx*vx+vy*vy)-1/np.sqrt(x*x+y*y)]; j+=1
    return out[:j]
@njit
def run_yoshida(nsteps, sample_stride, h):
    w1=1.0/(2.0-2.0**(1.0/3.0)); w0=-2.0**(1.0/3.0)/(2.0-2.0**(1.0/3.0))
    x,y,vx,vy=0.1,0.0,0.0,np.sqrt(19.0); ax,ay=acc(x,y)
    nout=nsteps//sample_stride+1; out=np.zeros((nout,6)); j=0
    out[j]=[0,x,y,vx,vy,0.5*(vx*vx+vy*vy)-1/np.sqrt(x*x+y*y)]; j+=1
    for i in range(1,nsteps+1):
        x,y,vx,vy,ax,ay=verlet_step(x,y,vx,vy,w1*h,ax,ay)
        x,y,vx,vy,ax,ay=verlet_step(x,y,vx,vy,w0*h,ax,ay)
        x,y,vx,vy,ax,ay=verlet_step(x,y,vx,vy,w1*h,ax,ay)
        if i%sample_stride==0:
            out[j]=[i,x,y,vx,vy,0.5*(vx*vx+vy*vy)-1/np.sqrt(x*x+y*y)]; j+=1
    return out[:j]
def exact_kepler(t):
    M=np.mod(t,PERIOD)/PERIOD*2*np.pi
    Ean=M.copy()
    for _ in range(30):
        Ean -= (Ean-E*np.sin(Ean)-M)/(1-E*np.cos(Ean))
    x=A*(np.cos(Ean)-E); y=A*np.sqrt(1-E*E)*np.sin(Ean)
    # velocities for mu=1, mean motion n=1
    dEdt=1/(1-E*np.cos(Ean))
    vx=-A*np.sin(Ean)*dEdt; vy=A*np.sqrt(1-E*E)*np.cos(Ean)*dEdt
    return x,y,vx,vy

def main():
    # Resolve the long-run diagnostics at 200 samples per orbit. The 10,000-step/orbit
    # run remains the main figure; 5,000 and 20,000 provide resolution checks.
    orbits=1000
    all_rows=[]
    plot_frames=[]
    for steps_per_orbit in [5000,10000,20000]:
        nsteps=orbits*steps_per_orbit
        stride=steps_per_orbit//200
        h=PERIOD/steps_per_orbit
        for name,fun in [('Velocity-Verlet',run_verlet),('RK4',run_rk4),('Yoshida4',run_yoshida)]:
            t0=time.perf_counter(); arr=fun(nsteps,stride,h); elapsed=time.perf_counter()-t0
            t=arr[:,0]*h; exact=np.array(exact_kepler(t)).T
            poserr=np.sqrt((arr[:,1]-exact[:,0])**2+(arr[:,2]-exact[:,1])**2)
            # Dense output prevents angular branch ambiguity in the unwrapping step.
            phn=np.unwrap(np.arctan2(arr[:,2],arr[:,1]))
            phe=np.unwrap(np.arctan2(exact[:,1],exact[:,0]))
            phase=np.abs(phn-phe)
            relE=np.abs((arr[:,5]-(-0.5))/(-0.5))
            row={'steps_per_orbit':steps_per_orbit,'orbits':orbits,'method':name,
                 'samples_per_orbit':200,'max_relative_energy_error':float(relE.max()),
                 'max_position_error_a':float(poserr.max()),'final_position_error_a':float(poserr[-1]),
                 'max_unwrapped_true_longitude_error_rad':float(phase.max()),
                 'final_unwrapped_true_longitude_error_rad':float(phase[-1]),'runtime_s':elapsed}
            all_rows.append(row)
            frame=pd.DataFrame({'time_orbits':t/PERIOD,'relative_energy_error':relE,
                                'position_error_a':poserr,'unwrapped_true_longitude_error_rad':phase})
            suffix=name.lower().replace('-','_')
            frame.to_csv(OUT / f'eccentric_1000_orbits_{steps_per_orbit}_steps_per_orbit_{suffix}_dense200.csv',index=False)
            if steps_per_orbit==10000:
                plot_frames.append(frame.assign(method=name))
            print(row,flush=True)
    comp=pd.DataFrame(all_rows)
    comp.to_csv(OUT / 'eccentric_1000_orbits_resolution_comparison.csv',index=False)
    df=pd.concat(plot_frames,ignore_index=True)
    fig,axs=plt.subplots(2,1,figsize=(8.2,7.2),sharex=True)
    for name,g in df.groupby('method'):
        # Plot one sample per orbit for legibility; all 200 samples/orbit are retained in CSVs.
        g=g.iloc[::200]
        axs[0].plot(g.time_orbits,g.relative_energy_error,label=name,lw=1.0)
        axs[1].plot(g.time_orbits,g.unwrapped_true_longitude_error_rad,label=name,lw=1.0)
    axs[0].set_yscale('log'); axs[0].set_ylabel('Relative energy error')
    axs[1].set_yscale('log'); axs[1].set_ylabel('Absolute unwrapped true-longitude error (rad)'); axs[1].set_xlabel('Elapsed orbital periods')
    axs[0].legend(ncol=3); axs[0].grid(True,which='both',alpha=.25); axs[1].grid(True,which='both',alpha=.25)
    fig.suptitle('Long-duration eccentric Kepler integration (e = 0.9; 10,000 steps/orbit)')
    fig.tight_layout(); fig.savefig(FIG / 'eccentric_long_duration.png',dpi=220,bbox_inches='tight'); plt.close(fig)
    # Figure-eight RK4 one period return test with truncated published ICs
    @njit
    def acc3(r):
        a=np.zeros_like(r)
        for i in range(3):
            for j in range(3):
                if i!=j:
                    dx=r[j,0]-r[i,0]; dy=r[j,1]-r[i,1]; q=(dx*dx+dy*dy)**(-1.5)
                    a[i,0]+=dx*q; a[i,1]+=dy*q
        return a
    @njit
    def deriv(y):
        r=y[:6].reshape((3,2)); v=y[6:].reshape((3,2)); a=acc3(r)
        return np.concatenate((v.ravel(),a.ravel()))
    @njit
    def rkstep(y,h):
        k1=deriv(y); k2=deriv(y+0.5*h*k1); k3=deriv(y+0.5*h*k2); k4=deriv(y+h*k3)
        return y+h*(k1+2*k2+2*k3+k4)/6
    r0=np.array([[.97000436,-.24308753],[-.97000436,.24308753],[0.,0.]])
    v0=np.array([[.466203685,.432365730],[.466203685,.432365730],[-.932407370,-.864731460]])
    y0=np.concatenate((r0.ravel(),v0.ravel()))
    T=6.32591398
    rows=[]
    for steps in [1000,2000,4000,8000,16000]:
        y=y0.copy(); hh=T/steps
        for i in range(steps): y=rkstep(y,hh)
        dr=np.linalg.norm(y[:6]-y0[:6]); dv=np.linalg.norm(y[6:]-y0[6:])
        rows.append({'steps_per_period':steps,'position_return_error':dr,'velocity_return_error':dv})
    pd.DataFrame(rows).to_csv(OUT / 'figure8_return_error.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    with open(OUT / 'longrun_environment.txt','w') as f:
        f.write('CPU model recorded for benchmark: AMD EPYC 9V74 80-Core Processor (3 logical CPUs available to the process)\n')
        f.write('Platform: '+platform.platform()+'\n')
        import numpy, scipy, numba
        f.write(f'Python: {platform.python_version()}\nNumPy: {numpy.__version__}\nSciPy: {scipy.__version__}\nNumba: {numba.__version__}\n')
if __name__=='__main__': main()
