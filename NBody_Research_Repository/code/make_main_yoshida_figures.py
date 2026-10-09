import sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import nbody_research as nb
OUT=str(ROOT/'results'/'figures')+'/'
G,AU,DAY,YEAR,MSUN,MEARTH=nb.G,nb.AU,nb.DAY,nb.YEAR,nb.MSUN,nb.MEARTH

def acc(r,m):
 d=r[None,:,:]-r[:,None,:]; d2=np.einsum('ijk,ijk->ij',d,d); np.fill_diagonal(d2,np.inf); return G*np.einsum('j,ijk,ij->ik',m,d,d2**-1.5)
def verlet_cached(r,v,m,h,a=None):
 if a is None:a=acc(r,m)
 rn=r+h*v+.5*h*h*a; an=acc(rn,m); vn=v+.5*h*(a+an); return rn,vn,an
def yoshida(r,v,m,h,a=None):
 c=1/(2-2**(1/3)); d=-2**(1/3)/(2-2**(1/3))
 if a is None:a=acc(r,m)
 for q in [c,d,c]:r,v,a=verlet_cached(r,v,m,q*h,a)
 return r,v,a
def rk4(r,v,m,h):
 def f(rr,vv):return vv,acc(rr,m)
 k1r,k1v=f(r,v);k2r,k2v=f(r+h*k1r/2,v+h*k1v/2);k3r,k3v=f(r+h*k2r/2,v+h*k2v/2);k4r,k4v=f(r+h*k3r,v+h*k3v)
 return r+h*(k1r+2*k2r+2*k3r+k4r)/6,v+h*(k1v+2*k2v+2*k3v+k4v)/6
def energy(r,v,m):return nb.total_energy(r[None,:,:],v[None,:,:],m)[0]
def lvec(r,v,m):return np.sum(np.cross(r,m[:,None]*v),axis=0)

# Five-year 2-day comparison, all four methods.
r0,v0,m=nb.sun_earth(); dt=2*DAY; nsteps=round(5*YEAR/dt); methods=['Euler','Velocity-Verlet','RK4','Yoshida4']; rows={}
for method in methods:
 r=r0.copy();v=v0.copy();a=acc(r,m); E0=energy(r,v,m);L0=lvec(r,v,m); ts=[]; ee=[]; ll=[]
 for i in range(nsteps+1):
  if i%5==0 or i==nsteps:
   ts.append(i*dt/YEAR);ee.append(abs((energy(r,v,m)-E0)/E0));ll.append(np.linalg.norm(lvec(r,v,m)-L0)/np.linalg.norm(L0))
  if i==nsteps:break
  if method=='Euler':r,v=r+dt*v,v+dt*acc(r,m)
  elif method=='Velocity-Verlet':r,v,a=verlet_cached(r,v,m,dt,a)
  elif method=='RK4':r,v=rk4(r,v,m,dt)
  else:r,v,a=yoshida(r,v,m,dt,a)
 rows[method]=(np.array(ts),np.array(ee),np.array(ll))
fig,ax=plt.subplots(figsize=(7.2,5.2))
for method,(t,e,l) in rows.items():ax.semilogy(t,np.maximum(e,1e-18),label=method)
ax.set_xlabel('Time (years)');ax.set_ylabel('Relative energy error');ax.set_title('Five-year integrator comparison (Δt = 2 days)');ax.grid(True,which='both',alpha=.25);ax.legend();fig.tight_layout();fig.savefig(OUT+'fig02_energy_error_integrators.png',dpi=320);plt.close(fig)
fig,ax=plt.subplots(figsize=(7.2,5.2))
for method,(t,e,l) in rows.items():ax.semilogy(t,np.maximum(l,1e-18),label=method)
ax.set_xlabel('Time (years)');ax.set_ylabel('Relative angular-momentum error');ax.set_title('Five-year angular-momentum comparison (Δt = 2 days)');ax.grid(True,which='both',alpha=.25);ax.legend();fig.tight_layout();fig.savefig(OUT+'fig03_angular_momentum_error_integrators.png',dpi=320);plt.close(fig)
# Circular Sun-Earth convergence; exact circular analytical relative orbit at each tested dt.
conv=[]
for dtday in [4,2,1,.5,.25,.125]:
 h=dtday*DAY; steps=round(YEAR/h); times=np.arange(steps+1)*h
 for method in methods:
  r=r0.copy();v=v0.copy();a=acc(r,m); maxpos=0.;E0=energy(r,v,m);maxE=0.
  for i in range(1,steps+1):
   if method=='Euler':r,v=r+h*v,v+h*acc(r,m)
   elif method=='Velocity-Verlet':r,v,a=verlet_cached(r,v,m,h,a)
   elif method=='RK4':r,v=rk4(r,v,m,h)
   else:r,v,a=yoshida(r,v,m,h,a)
   if i%max(1,steps//1000)==0 or i==steps:
    t=i*h; ra=nb.analytic_two_body(np.array([t]),MSUN,MEARTH)
    maxpos=max(maxpos,np.linalg.norm(r[1]-ra[0,1])/AU)
    maxE=max(maxE,abs((energy(r,v,m)-E0)/E0))
  conv.append([method,dtday,maxpos,maxE])
df=pd.DataFrame(conv,columns=['method','dt_days','max_position_error_au','max_energy_error']);df.to_csv(ROOT/'results'/'tables'/'convergence_results.csv',index=False)
fig,ax=plt.subplots(figsize=(7.2,5.2))
for method,g in df.groupby('method'):ax.loglog(g.dt_days,g.max_position_error_au,'o-',label=method)
ax.set_xlabel('Timestep (days)');ax.set_ylabel('Maximum Earth position error (au)');ax.set_title('Timestep convergence including Yoshida4');ax.invert_xaxis();ax.grid(True,which='both',alpha=.25);ax.legend();fig.tight_layout();fig.savefig(OUT+'fig04_timestep_convergence_position.png',dpi=320);plt.close(fig)
fig,ax=plt.subplots(figsize=(7.2,5.2))
for method,g in df.groupby('method'):ax.loglog(g.dt_days,g.max_energy_error,'o-',label=method)
ax.set_xlabel('Timestep (days)');ax.set_ylabel('Maximum relative energy error');ax.set_title('Energy-error convergence including Yoshida4');ax.invert_xaxis();ax.grid(True,which='both',alpha=.25);ax.legend();fig.tight_layout();fig.savefig(OUT+'fig05_timestep_convergence_energy.png',dpi=320);plt.close(fig)
print('wrote figures and convergence CSV')
