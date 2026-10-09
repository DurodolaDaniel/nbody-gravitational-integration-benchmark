import sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code')); import nbody_research as nb
G,AU,DAY,YEAR,MSUN,MEARTH=nb.G,nb.AU,nb.DAY,nb.YEAR,nb.MSUN,nb.MEARTH
r0,v0,m=nb.sun_earth(); rows=[]
for steps in [1000,3000,10000,30000,100000,300000,1000000]:
 h=YEAR/steps;r=r0.copy();v=v0.copy();maxerr=0.
 stride=max(1,steps//3000)
 for i in range(1,steps+1):
  a=nb.acceleration(r,m); rn=r+h*v; vn=v+h*a; r,v=rn,vn
  if i%stride==0 or i==steps:
   exact=nb.analytic_two_body(np.array([i*h]),MSUN,MEARTH)[0]
   maxerr=max(maxerr,np.linalg.norm(r[1]-exact[1])/AU)
 rows.append([steps,h/DAY,maxerr])
df=pd.DataFrame(rows,columns=['steps_per_year','timestep_days','max_position_error_au']);df.to_csv(ROOT/'results'/'tables'/'euler_fine_convergence.csv',index=False)
# fit fine steps >= 10000
q=df[df.steps_per_year>=10000];slope=np.polyfit(np.log(1/q.steps_per_year),np.log(q.max_position_error_au),1)[0]
fig,ax=plt.subplots(figsize=(7.2,5.2));ax.loglog(df.steps_per_year,df.max_position_error_au,'o-');ax.set_xlabel('Steps per year');ax.set_ylabel('Maximum Earth position error (au)');ax.set_title(f'Forward Euler convergence on circular two-body orbit (slope={slope:.3f})');ax.grid(True,which='both',alpha=.25);fig.tight_layout();fig.savefig(ROOT/'results'/'figures'/'euler_fine_convergence.png',dpi=320);plt.close(fig)
print(df.to_string(index=False));print('slope',slope)
