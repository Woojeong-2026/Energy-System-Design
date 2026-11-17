"""Independent numerical reference for the assignment, not a MATLAB runtime.

Uses Cantera 3.0.0 Gibbs properties and a bracketed equilibrium-composition root.
The MATLAB model uses an analytic recycle root. Property data are shared, so this
checks implementation consistency, not the accuracy of the thermochemical data.
"""
import argparse
import csv
import json
import math
from functools import lru_cache
from pathlib import Path

import cantera as ct
import numpy as np

DEFAULT = dict(P_bar=5., Tcell_C=850., TBch_C=650., T0_C=25.,
               nCO2_feed=1., X_soec=.70, dV_loss=.10, eff_hex=.85,
               f_link=.85, dT_link=20., Ru=8314., F=96485., p0=1e5, LHV_C=393.5e6)
METRICS = ('K_B','yCO','E0','E_rev','V','V_tn','I','W_cell','Q_hex',
           'Qmax_hex','Q_link','Qel_tot','Qrej_tot','P_in','nC','eta_ch')


def bisect(f, lo, hi):
    flo, fhi = f(lo), f(hi)
    if flo == 0: return lo
    if fhi == 0: return hi
    if not (math.isfinite(flo) and math.isfinite(fhi)) or flo*fhi > 0:
        raise ValueError('Root is not bracketed')
    for _ in range(100):
        mid = (lo+hi)/2
        fm = f(mid)
        if not math.isfinite(fm): raise ValueError('Nonfinite root residual')
        if fm == 0 or hi-lo < 1e-12*max(1., abs(mid)): return mid
        if flo*fm < 0: hi = mid
        else: lo, flo = mid, fm
    return (lo+hi)/2


def config(overrides=None):
    repo_only={'P0_bar','TBdis_C','nC_dis','U_sofc','lambda_O2','xO2_air','xN2_air'}
    unknown=set(overrides or {})-set(DEFAULT)-repo_only
    if unknown: raise ValueError('Unknown config fields: '+', '.join(sorted(unknown)))
    c = DEFAULT | {k: v for k, v in (overrides or {}).items() if k in DEFAULT}
    for k,v in c.items():
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):
            raise ValueError(f'{k} must be a finite real scalar')
    if not (c['P_bar']>0 and c['nCO2_feed']>0 and 0<c['X_soec']<1):
        raise ValueError('Require P>0, feed>0, 0<X_soec<1')
    if not (0<=c['eff_hex']<=1 and 0<=c['f_link']<=1 and c['dV_loss']>=0 and c['dT_link']>=0):
        raise ValueError('Invalid heat recovery or voltage loss')
    if any(c[k]<=0 for k in ('Ru','F','LHV_C')) or abs(c['p0']-1e5)>1e-6:
        raise ValueError('Invalid constants or standard pressure')
    ts=[c[k]+273.15 for k in ('T0_C','Tcell_C','TBch_C')]
    if min(ts)<200 or max(ts)>3500 or min(ts[1:])<=ts[0]:
        raise ValueError('Temperature outside fixed-process model domain')
    return c


class Thermo:
    def __init__(self):
        self.gas=ct.Solution('gri30.yaml'); self.gr=ct.Solution('graphite.yaml')

    @lru_cache(maxsize=30000)
    def props(self, t):
        h=[]; g=[]; mw=[]
        for name in ('CO','CO2','O2','N2'):
            self.gas.TPX=t,1e5,{name:1}
            h.append(self.gas.enthalpy_mole); g.append(self.gas.gibbs_mole)
            mw.append(self.gas.mean_molecular_weight)
        self.gr.TP=t,1e5
        h.append(self.gr.enthalpy_mole); g.append(self.gr.gibbs_mole)
        mw.append(self.gr.mean_molecular_weight)
        return np.array(h),np.array(g),np.array(mw)

    def H(self,n,t): return float(np.dot(n,self.props(float(t))[0]))

    def invH(self,n,target,lo,hi):
        return bisect(lambda t:(self.H(n,t)-target)/1e6,lo-1,hi+1)


def solve(overrides=None, th=None):
    c=config(overrides); th=th or Thermo()
    t0=c['T0_C']+273.15; tc=c['Tcell_C']+273.15; tb=c['TBch_C']+273.15
    p=c['P_bar']*1e5; pr=p/c['p0']; f=c['nCO2_feed']; x=c['X_soec']
    gb=th.props(tb)[1]; lnK=-(gb[4]+gb[1]-2*gb[0])/(c['Ru']*tb)
    y=bisect(lambda y:math.log1p(-y)-2*math.log(y)-math.log(pr)-lnK,1e-14,1-1e-14)
    # Overall carbon and CO balances determine extent and CO2 recycle.
    co2in=2*f/x; co2rec=co2in-f; corec=co2rec*y/(1-y)
    n=np.zeros((13,5)); n[0,1]=f; n[8,:2]=corec,co2rec
    n[1]=n[0]+n[8]; n[2]=n[3]=n[1]
    n[4,:2]=corec+2*f,co2in*(1-x); n[5]=n[4]
    n[6]=n[7]=n[8]; n[9,4]=n[10,4]=f; n[11,2]=n[12,2]=f
    t=np.array([t0,t0,0,tc,tc,tb,tb,0,t0,tb,t0,tc,t0])
    qmax=min(th.H(n[6],tb)-th.H(n[6],t0),th.H(n[1],tb)-th.H(n[1],t0))
    qhex=c['eff_hex']*qmax
    t[2]=th.invH(n[1],th.H(n[1],t0)+qhex,t0,tb)
    t[7]=th.invH(n[6],th.H(n[6],tb)-qhex,t0,tb)
    gc=th.props(tc)[1]; e0=-(gc[1]-gc[0]-.5*gc[2])/(2000*c['F'])
    er=e0+c['Ru']*tc/(2000*c['F'])*math.log(n[4,0]/n[4,1]*math.sqrt(pr))
    v=er+c['dV_loss']; current=4000*c['F']*f; w=current*v
    h=np.array([th.H(nn,tt) for nn,tt in zip(n,t)])
    q=np.array([h[3]-h[2],h[4]+h[11]-h[3]-w,h[5]-h[4],
                h[6]+h[9]-h[5],h[8]-h[7],h[10]-h[9],h[12]-h[11]])
    ql=np.zeros(7); link=0.
    if q[3]<0 and q[1]>0 and tb-tc>=c['dT_link']-1e-9:
        link=min(-c['f_link']*q[3],q[1]); ql[3]=-link; ql[1]=link
    elif q[1]<0 and q[3]>0 and tc-tb>=c['dT_link']-1e-9:
        link=min(-c['f_link']*q[1],q[3]); ql[1]=-link; ql[3]=link
    qe=np.maximum(q,0)-np.maximum(ql,0); qr=np.maximum(-q,0)-np.maximum(-ql,0)
    pin=w+sum(qe)
    m=dict(zip(METRICS,[math.exp(lnK),y,e0,er,v,(h[4]+h[11]-h[3])/current,
          current,w,qhex,qmax,link,sum(qe),sum(qr),pin,f,f*c['LHV_C']/pin]))
    return dict(schema_version=1,cfg=c,n=n.tolist(),T=t.tolist(),P=[p]*13,
                q=q.tolist(),ql=ql.tolist(),qe=qe.tolist(),qr=qr.tolist(),metrics=m)


def audit(packet,th=None):
    """Compare supplied raw outputs to a newly solved, same-input reference.

    Rejects malformed schema/nonfinite values rather than reporting a partial pass.
    Adds conservation and internal HEX checks independent of the output fields.
    """
    if packet.get('schema_version')!=1: raise ValueError('Unsupported schema')
    if not isinstance(packet.get('cfg'),dict) or not set(DEFAULT)<=set(packet['cfg']):
        raise ValueError('All required charging configuration fields must be present')
    th=th or Thermo(); c=config(packet['cfg']); ref=solve(c,th); rows=[]
    def check(name,res,limit):
        rows.append(dict(Check=name,Residual=float(res),Limit=float(limit),
                         Pass=bool(math.isfinite(res) and res<=limit)))
    # Validation tolerances: 0.1 W + relative 1e-7, 1e-5 K, 1e-4 Pa.
    abs_tol=dict(n=1e-9,T=1e-5,P=1e-4,q=.1,ql=.1,qe=.1,qr=.1)
    for key,atol in abs_tol.items():
        a=np.asarray(packet[key],dtype=float); b=np.asarray(ref[key],dtype=float)
        if a.shape!=b.shape or not np.isfinite(a).all(): raise ValueError(f'Invalid {key}')
        check(key+' reference scaled error',float(np.max(np.abs(a-b)/(atol+1e-7*np.abs(b)))),1.)
    for key in METRICS:
        a=packet['metrics'][key]; b=ref['metrics'][key]
        if not isinstance(a,(int,float)) or not math.isfinite(a): raise ValueError(f'Invalid {key}')
        atol=.1 if key in ('I','W_cell','Q_hex','Qmax_hex','Q_link','Qel_tot','Qrej_tot','P_in') else 1e-9
        check(key+' reference',abs(a-b),atol+1e-7*abs(b))
    n=np.asarray(packet['n']); t=np.asarray(packet['T']); m=packet['metrics']
    h=np.array([th.H(nn,tt) for nn,tt in zip(n,t)]); qr=np.asarray(packet['qr'])
    A=np.array([[1,1,0,0,1],[1,2,2,0,0],[0,0,0,2,0]])
    check('Overall C/O/N balance',max(abs(A@(n[0]-n[10]-n[12]))),1e-8*c['nCO2_feed']+1e-9)
    check('Overall energy balance W',abs(h[0]+m['P_in']-sum(qr)-h[10]-h[12]),.1+1e-8*abs(m['P_in']))
    check('HEX energy balance W',abs(h[6]+h[1]-h[7]-h[2]),.1)
    q=np.asarray(packet['q'])
    units=[('Mixer',[0,8],[1],0.,0.),('HEX1',[1,6],[2,7],0.,0.),
        ('H1',[2],[3],q[0],0.),('SOEC',[3],[4,11],q[1],m['W_cell']),
        ('Trim',[4],[5],q[2],0.),('Reactor',[5],[6,9],q[3],0.),
        ('Recycle cooler',[7],[8],q[4],0.),('C cooler',[9],[10],q[5],0.),
        ('O2 cooler',[11],[12],q[6],0.)]
    for label,ins,outs,qq,ww in units:
        atoms=A@(n[ins].sum(axis=0)-n[outs].sum(axis=0))
        check(label+' C/O/N',max(abs(atoms)),1e-8*c['nCO2_feed']+1e-9)
        check(label+' energy W',abs(sum(h[ins])+qq+ww-sum(h[outs])),.1+1e-8*abs(m['P_in']))
    qhex=h[6]-h[7]; gaps=[]
    for q in np.linspace(0,qhex,201):
        hot=th.invH(n[6],h[7]+q,t[1],t[6]); cold=th.invH(n[1],h[1]+q,t[1],t[6])
        gaps.append(hot-cold)
    check('HEX positive temperature difference',float(min(gaps)<=0),0.)
    def entropy(nn,tt,pp):
        hh,gg,_=th.props(float(tt)); s=(hh-gg)/tt; ng=sum(nn[:4]); value=nn[4]*s[4]
        for i in range(4):
            if nn[i]>0: value+=nn[i]*(s[i]-ct.gas_constant*math.log(nn[i]/ng*pp/1e5))
        return value
    ss=[entropy(nn,tt,pp) for nn,tt,pp in zip(n,t,packet['P'])]
    sgen=ss[2]+ss[7]-ss[1]-ss[6]
    check('HEX nonnegative entropy generation',max(0,-sgen),1e-3)
    return rows,dict(dTmin_K=min(gaps),H_storage_W=h[10]+h[12]-h[0],
                    energy_residual_W=h[0]+m['P_in']-sum(qr)-h[10]-h[12],Sgen_hex_W_K=sgen)


def write_csv(path, rows):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def case_list():
    return [('base',{})]+[(f'P{p}',dict(P_bar=p)) for p in (1,2,5,10)]+[
        (f'T{t}',dict(Tcell_C=t)) for t in (750,800,850,900)]


def generate(out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True); th=Thermo(); summary=[]
    for label,overrides in case_list():
        p=solve(overrides,th); rows,extra=audit(p,th); m=p['metrics']
        if not all(r['Pass'] for r in rows): raise AssertionError(label)
        (out/(label+'.json')).write_text(json.dumps(p,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        write_csv(out/(label+'_checks.csv'),rows)
        streamrows=[]
        mw=th.props(298.15)[2]
        for k,(nn,tt,pp) in enumerate(zip(p['n'],p['T'],p['P']),1):
            row=dict(No=k,T_C=tt-273.15,P_bar=pp/1e5,n_kmol_s=sum(nn),m_kg_s=float(np.dot(nn,mw)))
            for name,val in zip(('CO','CO2','O2','N2','C'),nn):
                row['n_'+name]=val; row['x_'+name]=val/sum(nn)
            streamrows.append(row)
        write_csv(out/(label+'_streams.csv'),streamrows)
        names=('Heater H1','SOEC','Trim 5->6','Carbon reactor','Recycle cooler','C cooler','O2 cooler')
        heatrows=[dict(Unit=name,Q_MW=p['q'][j]/1e6,Qlink_MW=p['ql'][j]/1e6,
                      Qel_MW=p['qe'][j]/1e6,Qrej_MW=p['qr'][j]/1e6) for j,name in enumerate(names)]
        write_csv(out/(label+'_heat.csv'),heatrows)
        summary.append(dict(Case=label,P_bar=p['cfg']['P_bar'],Tcell_C=p['cfg']['Tcell_C'],
            K_B=m['K_B'],yCO7=m['yCO'],n4_kmol_s=sum(p['n'][3]),V_ch=m['V'],
            Wcell_MW=m['W_cell']/1e6,Qel_MW=m['Qel_tot']/1e6,Pin_MW=m['P_in']/1e6,
            eta_ch_pct=100*m['eta_ch'],Qrej_MW=m['Qrej_tot']/1e6,Qlink_MW=m['Q_link']/1e6,
            **extra,AllChecksPass=True))
    write_csv(out/'summary.csv',summary)
    (out/'provenance.json').write_text(json.dumps(dict(runtime='Python independent reference',
        cantera_version=ct.__version__,model_baseline_commit='fb28621885af0df18ac6d20412ad4db788992efe',
        reference_pressure_Pa=1e5,matlab_executed=False),indent=2)+'\n',encoding='utf-8')
    return summary


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=Path('reference_results'))
    ap.add_argument('--compare',type=Path,help='raw.json exported by MATLAB')
    args=ap.parse_args()
    if args.compare:
        p=json.loads(args.compare.read_text(encoding='utf-8')); rows,extra=audit(p)
        args.out.mkdir(parents=True,exist_ok=True)
        write_csv(args.out/'comparison.csv',rows)
        print(json.dumps(extra,indent=2)); failures=[r for r in rows if not r['Pass']]
        print(f'{len(rows)-len(failures)}/{len(rows)} checks passed')
        return int(bool(failures))
    for row in generate(args.out):
        print(f"{row['Case']:>5}: Pin={row['Pin_MW']:.4f} MW, eta_ch={row['eta_ch_pct']:.4f}%, dTmin={row['dTmin_K']:.3f} K")
    return 0


if __name__=='__main__':
    try: raise SystemExit(main())
    except (ValueError,KeyError,TypeError) as exc:
        raise SystemExit('Invalid input: '+str(exc)) from exc
