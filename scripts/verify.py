"""Independently recompute the retained participant-table associations and FDR.

No EEG is accessed. Run from any directory. Figures use residual ranks, not IQ
predictions. This script does not re-select the historical feature family.
"""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import rankdata, t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results/retained'
OUT = ROOT / 'results/verification'

def read(name):
    with (DATA / name).open(newline='') as f:
        return list(csv.DictReader(f))

def array(rows, name):
    return np.array([float(r[name]) for r in rows])

def residual(x, controls):
    design = np.column_stack([np.ones(len(x))] + [rankdata(c) for c in controls])
    ranked = rankdata(x)
    return ranked - design @ np.linalg.lstsq(design, ranked, rcond=None)[0]

def correlate(x, y, controls):
    mask = np.isfinite(x) & np.isfinite(y)
    for c in controls:
        mask &= np.isfinite(c)
    x, y, controls = x[mask], y[mask], [c[mask] for c in controls]
    rx, ry = residual(x, controls), residual(y, controls)
    r = float(np.corrcoef(rx, ry)[0, 1])
    df = len(x) - len(controls) - 2
    p = float(2 * t.sf(abs(r) * np.sqrt(df / max(1-r*r, 1e-9)), df))
    return r, p, len(x)

def fdr(p):
    p = np.array(p)
    order = np.argsort(p)
    qsort = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    q = np.empty(len(p)); q[order] = np.minimum(qsort, 1)
    return q

def verify_family(features, evaluations, benchmark=False):
    computed = []
    errors = []
    for row in evaluations:
        group = features
        if benchmark:
            group = [r for r in group if all(r[k] == row[k] for k in ('condition','band','channel_set'))]
            if row['analysis'].startswith('young'):
                group = [r for r in group if r['cohort'] == 'young']
        x, y, w = (array(group,k) for k in (row['feature'],'LPS','WST'))
        controls = [array(group,k) for k in row['controls'].split(';')]
        r,p,n = correlate(x,y,controls)
        rw,pw,_ = correlate(x,w,controls)
        rlw,plw,_ = correlate(x,y,controls+[w])
        values = dict(n=n,r_lps=r,p_lps=p,r_wst=rw,p_wst=pw,r_lps_plus_wst=rlw,p_lps_plus_wst=plw)
        computed.append(values)
    for pk,qk in [('p_lps','q_lps'),('p_lps_plus_wst','q_lps_plus_wst')]:
        if qk in evaluations[0]:
            for row,q in zip(computed,fdr([r[pk] for r in computed])):
                row[qk]=float(q)
    for expected,actual in zip(evaluations,computed):
        for key,val in actual.items():
            ref=float(expected[key])
            if not np.isclose(val,ref,rtol=1e-9,atol=1e-10,equal_nan=True):
                raise AssertionError((expected['feature'],expected['analysis'],key,val,ref))
            errors.append(abs(val-ref))
    return dict(rows=len(evaluations),values_compared=len(errors),max_absolute_difference=max(errors))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    spatial=read('spatial_subjects_full_full.csv')
    benchmark=read('bfocus_features_full_grid24.csv')
    seval=read('spatial_eval_full_full.csv')
    beval=read('bfocus_eval_full_grid24.csv')
    assert len(spatial)==111 and len({r['sid'] for r in spatial})==111
    report={'scope':'Statistical recomputation from retained participant features; no raw-EEG rerun',
            'spatial':verify_family(spatial,seval),
            'benchmark':verify_family(benchmark,beval,True)}
    rel=read('spatial_reliability_full_full.csv')
    for row in rel:
        a,b=('odd','even') if row['split']=='odd_even' else ('first','second')
        prefix=row['transform']+'_'+row['aperture']+'__'
        r,_,n=correlate(array(spatial,prefix+a+'_score'),array(spatial,prefix+b+'_score'),[])
        assert n==int(row['n'])
        np.testing.assert_allclose([r,2*r/(1+r)],[float(row['split_r']),float(row['spearman_brown'])],rtol=1e-9,atol=1e-10)
    report['reliability_rows_verified']=len(rel)
    feature='raw_posterior_broad__full_score'
    row=next(r for r in seval if r['feature']==feature and r['analysis']=='strict')
    controls=[array(spatial,k) for k in row['controls'].split(';')]
    rx=residual(array(spatial,feature),controls)
    ry=residual(array(spatial,'LPS'),controls)
    rw=residual(array(spatial,'WST'),controls)
    report['headline']=row
    # Participant bootstrap, with rank residualization recomputed in every sample.
    # Descriptive uncertainty conditional on the historically selected score.
    rng=np.random.default_rng(20260905)
    boots=[]
    x,y=array(spatial,feature),array(spatial,'LPS')
    for _ in range(3000):
        ix=rng.integers(0,len(x),len(x))
        boots.append(correlate(x[ix],y[ix],[c[ix] for c in controls])[0])
    report['headline_bootstrap']={'replicates':3000,'seed':20260905,'percentile_95_interval':np.quantile(boots,[.025,.975]).tolist(),'interpretation':'Conditional on historical score choice; not adjusted for discovery selection'}
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(13,4.1))
    colors=['#277a83' if r['cohort']=='young' else '#b86143' for r in spatial]
    for ax,target,label,r in [(axes[0],ry,'Fluid reasoning (LPS)',row['r_lps']),(axes[1],rw,'Vocabulary (WST)',row['r_wst'])]:
        ax.scatter(rx,target,c=colors,s=22,alpha=.75,edgecolors='none')
        slope,intercept=np.polyfit(rx,target,1)
        domain=np.array([rx.min(),rx.max()]); ax.plot(domain,slope*domain+intercept,color='#273544',lw=1.5)
        ax.set(xlabel='Theta concentration: residual rank',ylabel=label+': residual rank',title=f'{label} | partial r = {float(r):.3f}')
    labels=['Raw','Smooth .35','Smooth .50','High-boost .10','High-boost .25','High-boost .50','High-boost .75']
    trans=['raw','smooth35','smooth50','deblur10','deblur25','deblur50','deblur75']
    values=[float(next(r for r in seval if r['feature']==f'{tr}_posterior_broad__full_score' and r['analysis']=='strict')['r_lps']) for tr in trans]
    axes[2].barh(labels[::-1],values[::-1],color=['#aaa']*6+['#277a83'])
    axes[2].set(xlim=(0,.5),xlabel='Partial rank correlation with LPS',title='Retained spatial ablations')
    fig.suptitle('Eyes-open theta lagged-field concentration and fluid reasoning | 111 LEMON participants',fontsize=13)
    fig.tight_layout(rect=(0,0,1,.94))
    for ext in ('png','svg'):
        fig.savefig(ROOT/f'figures/correlate.{ext}',dpi=180,bbox_inches='tight')
    report['passed']=True
    (OUT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
