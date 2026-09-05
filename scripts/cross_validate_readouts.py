"""New participant-level cross-validation of fixed, retained EEG readouts.

This calibrates scalar scores; it does not re-fit or replace their EEG operators.
All calibration and predictor scaling use training subjects only. Historical
feature selection is not nested. No target values are imputed.
"""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results/retained'
OUT=ROOT/'results/cross_validation'
REPEATS=20
SEED=20260905

def read(name):
    with (DATA/name).open(newline='') as f:return list(csv.DictReader(f))

def write(path, rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def get(rows,key):return np.array([float(r[key]) for r in rows])

def fit_predict(x,y,train,test):
    if x.shape[1]==0:return np.full(len(test),y[train].mean())
    mean=x[train].mean(0);sd=x[train].std(0)
    sd=np.where(sd>1e-12,sd,1)
    a=np.column_stack([np.ones(len(train)),(x[train]-mean)/sd])
    b=np.column_stack([np.ones(len(test)),(x[test]-mean)/sd])
    return b @ np.linalg.lstsq(a,y[train],rcond=None)[0]

def metrics(y,p,b):
    denom=float(np.sum((y-y.mean())**2))
    r2=1-float(np.sum((y-p)**2))/denom
    base=1-float(np.sum((y-b)**2))/denom
    return dict(pearson_r=float(np.corrcoef(y,p)[0,1]),spearman_r=float(spearmanr(y,p).statistic),r2=r2,baseline_r2=base,delta_r2=r2-base,rmse=float(np.sqrt(np.mean((y-p)**2))),baseline_rmse=float(np.sqrt(np.mean((y-b)**2))))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    spatial=read('spatial_subjects_full_full.csv')
    bench=[r for r in read('bfocus_features_full_grid24.csv') if (r['condition'],r['band'],r['channel_set'])==('EO','theta','posterior')]
    specs=[('theta_concentration_broad',spatial,'raw_posterior_broad__full_score',['cohort_code','age','raw_posterior_broad__theta_power','raw_posterior_broad__theta_static_plv']),
           ('theta_concentration_lateral',spatial,'raw_lateral_posterior__full_score',['cohort_code','age','raw_lateral_posterior__theta_power','raw_lateral_posterior__theta_static_plv']),
           ('theta_focus_ratio',bench,'bfocus_ratio',['cohort_code','age','band_power','static_plv'])]
    sids=sorted(r['sid'] for r in spatial)
    assert len(sids)==111 and len(set(sids))==111
    n=len(sids);allidx=np.arange(n)
    plans=[('leave_one_out',0,[[i] for i in range(n)])]
    rng=np.random.default_rng(SEED)
    for rep in range(REPEATS):
        plans.append(('repeated_5fold',rep,[list(map(int,x)) for x in np.array_split(rng.permutation(n),5)]))
    split_records=[]
    for protocol,rep,folds in plans:
        for fold,idx in enumerate(folds):
            split_records.extend(dict(protocol=protocol,repeat=rep,fold=fold,sid=sids[i]) for i in idx)
    write(OUT/'fold_assignments.csv',split_records)
    predictions=[];summaries=[]
    leakage_checks=0
    for name,rows,feature,ckeys in specs:
        byid={r['sid']:r for r in rows};assert set(byid)==set(sids)
        rows=[byid[sid] for sid in sids]
        score=get(rows,feature)[:,None];cov=np.column_stack([get(rows,k) for k in ckeys])
        for target in ('LPS','WST'):
            y=get(rows,target)
            for model in ('score_only','score_plus_covariates'):
                b=cov if model=='score_plus_covariates' else np.empty((n,0))
                x=np.column_stack([b,score])
                assert np.isfinite(x).all() and np.isfinite(y).all()
                tr=allidx[1:];te=allidx[:1]
                perturbed=y.copy();perturbed[0]+=10000
                np.testing.assert_array_equal(fit_predict(x,y,tr,te),fit_predict(x,perturbed,tr,te))
                leakage_checks+=1
                for protocol,rep,folds in plans:
                    pred=np.full(n,np.nan);base=np.full(n,np.nan)
                    for fold,test0 in enumerate(folds):
                        test=np.array(test0);train=np.setdiff1d(allidx,test)
                        pred[test]=fit_predict(x,y,train,test);base[test]=fit_predict(b,y,train,test)
                        predictions.extend(dict(score=name,target=target,model=model,protocol=protocol,repeat=rep,fold=fold,sid=sids[i],observed=y[i],prediction=pred[i],baseline_prediction=base[i]) for i in test)
                    assert np.isfinite(pred).all() and np.isfinite(base).all()
                    summaries.append(dict(score=name,target=target,model=model,protocol=protocol,repeat=rep,n=n,**metrics(y,pred,base)))
    write(OUT/'predictions.csv',predictions);write(OUT/'run_metrics.csv',summaries)
    aggregate=[]
    for name,_,_,_ in specs:
        for target in ('LPS','WST'):
            for model in ('score_only','score_plus_covariates'):
                for protocol in ('leave_one_out','repeated_5fold'):
                    rows=[r for r in summaries if (r['score'],r['target'],r['model'],r['protocol'])==(name,target,model,protocol)]
                    item=dict(score=name,target=target,model=model,protocol=protocol,n=n,repeats=len(rows))
                    for key in ('pearson_r','spearman_r','r2','baseline_r2','delta_r2','rmse','baseline_rmse'):
                        v=[r[key] for r in rows]
                        item[key+'_median']=float(np.median(v));item[key+'_q05']=float(np.quantile(v,.05));item[key+'_q95']=float(np.quantile(v,.95))
                    aggregate.append(item)
    write(OUT/'summary.csv',aggregate)
    manifest={'scope':'New cross-validation from retained fixed participant EEG scores; no EEG extraction rerun or external cohort',
        'subjects':n,'targets':['LPS','WST'],'scores':[dict(name=name,feature=f,controls=c) for name,_,f,c in specs],
        'protocols':{'leave_one_out':111,'repeated_5fold':{'folds':5,'repeats':REPEATS,'seed':SEED}},
        'models':['score_only versus training-mean baseline','score_plus_covariates versus covariate-only baseline'],
        'calibration':'ordinary least squares on training subjects; linear scalar readout with intercept; scaling fitted in training fold',
        'selection':'scores fixed from published historical lineage before this run; no choice of best model by outcome; historical score discovery not nested',
        'missing_data':'all inputs and targets finite; no target imputation',
        'target_perturbation_invariance_checks':leakage_checks,'prediction_rows':len(predictions),'run_metric_rows':len(summaries),'summary_rows':len(aggregate),
        'r2_definition':'1 - pooled out-of-fold SSE / full evaluation sample centered SST; baseline_r2 uses same denominator',
        'repeat_intervals':'5th to 95th percentiles across random partitions; not confidence intervals or new participants'}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
    for r in aggregate:
        if r['protocol']=='leave_one_out':print(r['score'],r['target'],r['model'],'r',round(r['pearson_r_median'],3),'R2',round(r['r2_median'],3),'deltaR2',round(r['delta_r2_median'],3))

if __name__=='__main__':main()
