"""Optional: pip install matplotlib. All figures use Python reference data."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from charging_reference import Thermo, write_csv


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--results',type=Path,default=Path('reference_results'))
    root=ap.parse_args().results
    rows=list(csv.DictReader((root/'summary.csv').open(encoding='utf-8')))
    fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for row,(ids,xkey,ykey,xlabel,ylabel) in enumerate([
        (range(1,5),'P_bar','eta_ch_pct','Pressure [bar]','Charge efficiency [%]'),
        (range(1,5),'P_bar','n4_kmol_s','Pressure [bar]','Cell fuel inlet [kmol/s]'),
        (range(5,9),'Tcell_C','eta_ch_pct','Cell temperature [degC]','Charge efficiency [%]'),
        (range(5,9),'Tcell_C','Pin_MW','Cell temperature [degC]','Charging input [MW]')]):
        a=ax.flat[row]; a.plot([float(rows[i][xkey]) for i in ids],[float(rows[i][ykey]) for i in ids],'o-',color='#126b82')
        a.set(xlabel=xlabel,ylabel=ylabel); a.grid(alpha=.25)
    fig.suptitle('Charging sweeps | independent Python / Cantera 3.0.0 reference')
    fig.savefig(root/'charging_sweeps.png',dpi=180); plt.close(fig)
    p=json.loads((root/'base.json').read_text()); th=Thermo(); n=np.array(p['n']); t=p['T']
    q=np.linspace(0,p['metrics']['Q_hex'],201)
    hot=[th.invH(n[6],th.H(n[6],t[7])+v,t[1],t[6]) for v in q]
    cold=[th.invH(n[1],th.H(n[1],t[1])+v,t[1],t[6]) for v in q]
    write_csv(root/'base_hex_profile.csv',[dict(q_MW=float(v/1e6),T_hot_C=h-273.15,T_cold_C=c-273.15,dT_K=h-c) for v,h,c in zip(q,hot,cold)])
    fig,ax=plt.subplots(figsize=(8,4.8),constrained_layout=True)
    ax.plot(q/1e6,np.array(hot)-273.15,label='Hot: stream 8 end to stream 7 end',color='#be4b38')
    ax.plot(q/1e6,np.array(cold)-273.15,label='Cold: stream 2 end to stream 3 end',color='#126b82')
    ax.fill_between(q/1e6,np.array(cold)-273.15,np.array(hot)-273.15,alpha=.08,color='gray')
    ax.set(xlabel='Cumulative heat from cold-inlet end [MW]',ylabel='Temperature [degC]',
           title='HEX 1 counterflow profile | reference, 201 points')
    ax.legend(fontsize=8); ax.grid(alpha=.25); fig.savefig(root/'base_hex_profile.png',dpi=180); plt.close(fig)


if __name__=='__main__': main()
