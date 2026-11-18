"""Reproduce teaching examples using the existing independent Cantera reference."""
from pathlib import Path
import sys,json,copy,math,csv
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from charging_reference import solve,audit,Thermo
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out',type=Path,default=Path('results/charging_lessons'))
ROOT=parser.parse_args().out
ROOT.mkdir(parents=True,exist_ok=True)
th=Thermo(); base=solve(th=th)
experiments=[
 ('pressure','압력','P_bar',[1,2,5,10],5,'bar','과제 필수'),
 ('cell','셀 온도','Tcell_C',[750,800,850,900],850,'°C','과제 필수'),
 ('loss','전압 손실','dV_loss',[.05,.10,.15],.10,'V','민감도'),
 ('reactor','충전 반응기 온도','TBch_C',[650,700,869,870,900],650,'°C','개선 탐색'),
 ('hex','HEX 유용도','eff_hex',[.75,.85,.95],.85,'','민감도'),
 ('conversion','SOEC 전환율','X_soec',[.5,.7,.9],.7,'','민감도')]
lessons=[]; flat=[]
for key,label,param,values,default,unit,scope in experiments:
    cases=[]
    for v in values:
        p=solve({param:v},th); checks,extra=audit(p,th)
        assert all(c['Pass'] for c in checks),(key,v)
        m=p['metrics']
        row=dict(value=v,Pin=m['P_in']/1e6,eta=100*m['eta_ch'],n4=sum(p['n'][3]),
            W=m['W_cell']/1e6,Qcell=p['q'][1]/1e6,Qrx=p['q'][3]/1e6,
            L=m['Q_link']/1e6,Qhex=m['Q_hex']/1e6,Qel=m['Qel_tot']/1e6,
            Qrej=m['Qrej_tot']/1e6,Tcell=p['cfg']['Tcell_C'],TB=p['cfg']['TBch_C'],
            K=m['K_B'],y=m['yCO'],V=m['V'],I=m['I']/1e6,dTmin=extra['dTmin_K'])
        row={k:round(float(vv),8) for k,vv in row.items()}
        cases.append(row); flat.append(dict(experiment=key,**row))
    lessons.append(dict(id=key,label=label,param=param,default=default,unit=unit,scope=scope,cases=cases))

def settle(p):
    q=np.array(p['q']); ql=np.array(p['ql'])
    p['qe']=(np.maximum(q,0)-np.maximum(ql,0)).tolist()
    p['qr']=(np.maximum(-q,0)-np.maximum(-ql,0)).tolist()
    m=p['metrics'];m['Qel_tot']=sum(p['qe']);m['Qrej_tot']=sum(p['qr'])
    m['P_in']=m['W_cell']+m['Qel_tot'];m['eta_ch']=m['nC']*p['cfg']['LHV_C']/m['P_in']

faults=[]
for key,label in [('normal','정상 계산'),('false_link','온도 조건을 무시한 열회수'),
                  ('current','kmol → mol 변환 누락'),('nernst','Nernst를 스트림 4로 평가'),
                  ('double_hex','HEX 열을 입력 전력에서 또 차감'),('carbon','탄소 출구 유량 10% 누락')]:
    p=copy.deepcopy(base);m=p['metrics']
    if key=='false_link':
        L=min(-p['q'][3]*.85,p['q'][1]);p['ql'][1]=L;p['ql'][3]=-L;m['Q_link']=L;settle(p)
    if key=='current':
        old=m['W_cell'];m['I']/=1000;m['W_cell']=m['I']*m['V'];p['q'][1]+=old-m['W_cell'];settle(p)
    if key=='nernst':
        c=p['cfg'];n=p['n'][3];old=m['W_cell']
        m['E_rev']=m['E0']+c['Ru']*(c['Tcell_C']+273.15)/(2000*c['F'])*math.log(n[0]/n[1]*math.sqrt(c['P_bar']))
        m['V']=m['E_rev']+c['dV_loss'];m['W_cell']=m['I']*m['V'];p['q'][1]+=old-m['W_cell'];settle(p)
    if key=='double_hex':
        m['P_in']-=m['Q_hex'];m['eta_ch']=m['nC']*p['cfg']['LHV_C']/m['P_in']
    if key=='carbon':
        p['n'][10][4]=.9;m['nC']=.9;m['eta_ch']=m['nC']*p['cfg']['LHV_C']/m['P_in']
    results,extra=audit(p,th);byname={r['Check']:r for r in results}
    rows=[]
    for name,title,unit in [
        ('Overall C/O/N balance','전체 원소 수지','kmol/s'),
        ('Overall energy balance W','전체 에너지 수지','W'),
        ('I reference','Faraday 전류','A'),
        ('E_rev reference','지정 조성의 Nernst 전압','V'),
        ('Q_link reference','열회수 적격 조건','W')]:
        r=byname[name];rows.append(dict(label=title,residual=r['Residual'],limit=r['Limit'],passed=r['Pass'],unit=unit))
    passed=all(r['Pass'] for r in results)
    assert passed==(key=='normal')
    if key in ('false_link','current','nernst'):
        assert byname['Overall energy balance W']['Pass'],'Instructional balanced-fault construction failed'
    faults.append(dict(id=key,label=label,Pin=m['P_in']/1e6,eta=m['eta_ch']*100,
        I=m['I']/1e6,L=m['Q_link']/1e6,passed=passed,rows=rows,
        pass_count=sum(r['Pass'] for r in results),check_count=len(results)))

(ROOT/'physical_lessons.json').write_text(json.dumps(lessons,ensure_ascii=False,separators=(',',':')),encoding='utf-8')
(ROOT/'fault_lessons.json').write_text(json.dumps(faults,ensure_ascii=False,separators=(',',':')),encoding='utf-8')
with (ROOT/'physical_lessons.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
summary={'physical_cases':len(flat),'checks_per_case':46,'fault_scenarios':len(faults),
         'all_intended_faults_detected':True,'balanced_faults_detected':['false_link','current','nernst'],
         'runtime':'Python 3.11 / Cantera 3.0.0','matlab_executed':False}
(ROOT/'lesson_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
for f in faults: print(f['id'],f['pass_count'],f['check_count'],round(f['Pin'],3),round(f['eta'],3))
