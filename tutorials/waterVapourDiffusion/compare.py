#!/usr/bin/env python3
"""Compare native H2O/N2 simulations with the dilute-background diffusion law."""
from pathlib import Path
import csv,json,math,re
import numpy as np

ROOT=Path(__file__).resolve().parent

def values(path,n=None):
    text=path.read_text()
    uniform=re.search(r'internalField\s+uniform\s+([^;]+);',text)
    if uniform:
        a=np.fromstring(uniform[1].replace('(',' ').replace(')',' '),sep=' ')
        if n is None:raise ValueError('Need mesh size for uniform field: '+str(path))
        return np.full(n,a[0]) if len(a)==1 else np.tile(a,(n,1))
    match=re.search(r'internalField\s+nonuniform\s+List<(\w+)>\s+\d+\s*\((.*?)\)\s*;',text,re.S)
    if not match:raise ValueError('Cannot read native field '+str(path))
    a=np.fromstring(match[2].replace('(',' ').replace(')',' '),sep=' ')
    return a.reshape(-1,3) if match[1]=='vector' else a

def block(text,name):
    text=re.sub(r'//[^\n]*','',text)
    start=re.search(r'\b'+name+r'\s*\{',text).end();i=start;depth=1
    while depth:
        if text[i]=='{':depth+=1
        if text[i]=='}':depth-=1
        i+=1
    return text[start:i-1]

def scalar(text,key):return float(re.search(r'\b'+key+r'\s+(?:\[[^\]]+\]\s*)?([\deE.+-]+)\s*;',text)[1])

def main():
    single=(ROOT/'legacy/constant/foam.dat').exists()
    region=Path('.') if single else Path('fluid')
    nitrogen=block((ROOT/'legacy/constant'/region/'foam.dat').read_text(),'N2')
    T=float(values(ROOT/'legacy/0'/region/'T',1)[0])
    p=float(values(ROOT/'legacy/0'/region/'p',1)[0])
    length=.020;R=8314.46261815324
    mu=scalar(nitrogen,'As')*T**1.5/(T+scalar(nitrogen,'Ts'))
    rho=p*scalar(nitrogen,'molWeight')/(R*T)
    dLegacy=mu/rho
    water=block((ROOT/'speciesDependent/constant'/region/'speciesDiffusionProperties').read_text(),'H2O')
    dSpecies=scalar(water,'D')*(T/scalar(water,'Tref'))**scalar(water,'temperatureExponent')*(scalar(water,'pref')/p)**scalar(water,'pressureExponent')
    expected={'legacy':dLegacy,'speciesDependent':dSpecies}
    out=ROOT/'postProcessing/comparison';out.mkdir(parents=True,exist_ok=True)
    summary=dict(temperature_K=T,pressure_Pa=p,background='N2',mu_N2_Pa_s=mu,rho_N2_kg_m3=rho,D_legacy_background_limit=dLegacy,D_H2O_N2_reference=dSpecies,reference_to_legacy_ratio=dSpecies/dLegacy,models={})
    histories={};profiles={};allPassed=True
    for name,d in expected.items():
        case=ROOT/name;positions=values(case/'0'/region/'C');n=len(positions)
        y=positions[:,1];cosine=np.cos(math.pi*y/length)
        initial=values(case/'0'/region/'H2O',n)
        initialAmplitude=np.dot(initial-initial.mean(),cosine)/np.dot(cosine,cosine)
        rows=[];closure=0.;meanChange=0.
        times=sorted([q for q in case.iterdir() if q.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',q.name)],key=lambda q:float(q.name))
        for time in times:
            waterValues=values(time/region/'H2O',n);nitrogenValues=values(time/region/'N2',n)
            amplitude=np.dot(waterValues-waterValues.mean(),cosine)/np.dot(cosine,cosine)
            t=float(time.name);theory=initialAmplitude*math.exp(-d*(math.pi/length)**2*t)
            rows.append((t,amplitude,theory))
            closure=max(closure,float(np.max(abs(waterValues+nitrogenValues-1))))
            meanChange=max(meanChange,float(abs(waterValues.mean()-initial.mean())))
        finalTime=rows[-1][0];finalAmplitude=rows[-1][1]
        assert finalTime==1, 'Run both cases to physical time 1 s'
        uniqueY=np.unique(y);dy=uniqueY[1]-uniqueY[0]
        eigen=4*math.sin(math.pi/(2*len(uniqueY)))**2/dy**2
        inferred=((initialAmplitude/finalAmplitude)**(.001/finalTime)-1)/(.001*eigen)
        relative=abs(inferred-d)/d
        passed=relative<.002 and closure<1e-9 and meanChange<3e-8
        allPassed &= passed
        summary['models'][name]=dict(passed=passed,inferredD_m2_s=inferred,referenceD_m2_s=d,inferredD_relative_error=relative,initialAmplitude=initialAmplitude,finalAmplitude=finalAmplitude,massFractionSumMaxError=closure,meanH2OMaxChange=meanChange)
        histories[name]=np.asarray(rows)
        final=values(case/f'{finalTime:g}'/region/'H2O',n)
        profiles[name]=(uniqueY,np.asarray([final[y==yp].mean() for yp in uniqueY]))
        with (out/(name+'-amplitude.csv')).open('w',newline='') as f:
            w=csv.writer(f, lineterminator="\n");w.writerow(['time_s','computedAmplitude','diluteAnalyticalAmplitude']);w.writerows(rows)
    actualRatio=summary['models']['speciesDependent']['inferredD_m2_s']/summary['models']['legacy']['inferredD_m2_s']
    contrastDifference=abs(summary['models']['legacy']['finalAmplitude']-summary['models']['speciesDependent']['finalAmplitude'])/summary['models']['legacy']['initialAmplitude']
    summary.update(actual_diffusivity_ratio=actualRatio,final_profile_difference_over_initial_amplitude=contrastDifference,passed=bool(allPassed and actualRatio>1.5 and contrastDifference>.10))
    with (out/'profiles.csv').open('w',newline='') as f:
        w=csv.writer(f, lineterminator="\n");w.writerow(['y_m','legacy_H2O_massFraction','speciesDependent_H2O_massFraction'])
        w.writerows(zip(profiles['legacy'][0],profiles['legacy'][1],profiles['speciesDependent'][1]))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
        for name,label in [('legacy','Legacy: Sc = 1'),('speciesDependent','H2O/N2 reference D')]:
            yy,field=profiles[name];hist=histories[name]
            axes[0].plot(field*1e6,yy*1000,label=label,lw=2)
            axes[1].plot(hist[:,0],hist[:,1]/hist[0,1],label=label,lw=2)
            axes[1].plot(hist[:,0],hist[:,2]/hist[0,1],'k:',alpha=.5)
        axes[0].set(xlabel='H2O mass fraction (ppm)',ylabel='Height (mm)',title='Computed profiles at 1 s')
        axes[1].set(xlabel='Time (s)',ylabel='Amplitude / initial amplitude',title='Decay rate (dotted: dilute theory)')
        for ax in axes:ax.grid(alpha=.25);ax.legend(fontsize=8)
        fig.savefig(out/'comparison.png',dpi=180);plt.close(fig)
    except ImportError:summary['plot']='Install matplotlib to create comparison.png; CSV/results remain available.'
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    raise SystemExit(0 if summary['passed'] else 1)

if __name__=='__main__':main()
