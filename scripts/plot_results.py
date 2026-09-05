"""Publication figures from complete retained results and new CV outputs."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'figures'
DOMAIN={'blood_chemistry':'Blood chemistry','blood_pressure':'Blood pressure','body_size':'Body measurements','cognition':'Cognition','personality_affect':'Personality / affect'}
COLORS=['#087e8b','#c75b39','#5167ad','#8b5aa6','#ad8325']
SCORES={'theta_concentration_broad':'Broad theta concentration','theta_concentration_lateral':'Lateral theta concentration','theta_focus_ratio':'Theta focus ratio'}
KS=[3,6,10,20]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','figure.facecolor':'white','savefig.facecolor':'white','svg.hashsalt':'flint-v030'})

def read(path):
    with (ROOT/path).open(newline='') as f:return list(csv.DictReader(f))

def save(fig,name):
    for ext in ('png','svg'):
        fig.savefig(FIG/f'{name}.{ext}',dpi=170,bbox_inches='tight',metadata={'Date':None} if ext=='svg' else None)
    plt.close(fig)

def finish(fig,title,subtitle):
    fig.suptitle(title,x=.07,ha='left',fontsize=17,y=.995)
    fig.text(.07,.918,subtitle,ha='left',fontsize=10,color='#555')
    fig.tight_layout(rect=(0,0,1,.90))

def main():
    FIG.mkdir(exist_ok=True)
    rows=read('results/context/phenotypes/domain_label_reconstruction.csv')
    domains=read('results/context/phenotypes/domain_recovery_summary.csv')
    cv=read('results/cross_validation/summary.csv')
    runs=read('results/cross_validation/run_metrics.csv')
    preds=read('results/cross_validation/predictions.csv')
    assert len(rows)==1164
    # Every phenotype, same K. Signed values retain inverted predictions.
    fig,ax=plt.subplots(figsize=(10,5.8));rng=np.random.default_rng(20260905)
    for i,(d,title) in enumerate(DOMAIN.items()):
        group=[r for r in rows if r['domain']==d and int(r['k_identity_dims'])==6]
        vals=np.array([float(r['cv_corr']) for r in group])
        ax.scatter(vals,i+rng.uniform(-.24,.24,len(vals)),s=20,alpha=.67,c=COLORS[i],edgecolors='none')
        ax.plot([np.median(vals)]*2,[i-.34,i+.34],c='#172c38',lw=2)
    ax.axvline(0,color='#bbb',lw=1)
    ax.set(yticks=range(5),yticklabels=[f'{v} (n={sum(r["domain"]==d and int(r["k_identity_dims"])==6 for r in rows)})' for d,v in DOMAIN.items()],xlabel='Reported leave-one-out Pearson r',xlim=(-.65,.65))
    ax.invert_yaxis()
    finish(fig,'The full phenotype map','All 291 recorded fields at K = 6 · ticks mark medians · original mean-imputed target scoring')
    save(fig,'phenotype_distribution')
    # Domain R2 paths, full set.
    fig,axes=plt.subplots(1,2,figsize=(12,5.5))
    for (d,title),color in zip(DOMAIN.items(),COLORS):
        dr={int(r['k_identity_dims']):r for r in domains if r['domain']==d}
        for ax,key in zip(axes,['loocv_r2','in_sample_r2']):
            ax.plot(KS,[float(dr[k][key]) for k in KS],marker='o',label=title,color=color)
            ax.axhline(0,color='#bbb',lw=.6);ax.set(xticks=KS,xlabel='Retained identity coordinates K',ylabel='Pooled standardized R²')
    axes[0].set_title('Leave-one-out');axes[1].set_title('In sample');axes[1].legend(frameon=False,fontsize=9)
    finish(fig,'Generalization across whole domains','Five domains · every tested K · domain R² includes every field, not only favorable rows')
    save(fig,'domain_generalization')
    # Full matrices, one figure per domain, no label removal.
    for d,title in DOMAIN.items():
        labels=sorted({r['label'] for r in rows if r['domain']==d})
        lookup={(r['label'],int(r['k_identity_dims'])):float(r['cv_corr']) for r in rows if r['domain']==d}
        matrix=np.array([[lookup[(lab,k)] for k in KS] for lab in labels])
        fig,ax=plt.subplots(figsize=(12,max(3.8,.235*len(labels)+1.6)))
        im=ax.imshow(matrix,aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-.65,vcenter=0,vmax=.65))
        ax.set(xticks=range(4),xticklabels=KS,yticks=range(len(labels)),yticklabels=[' / '.join(l.split('__')[-2:]) for l in labels],xlabel='Retained identity coordinates K')
        ax.tick_params(axis='y',labelsize=8)
        ax.set_title(f'{title}: all {len(labels)} fields\nReported leave-one-out Pearson r · source-name order',loc='left',pad=13)
        if len(labels)<=40:
            for i in range(len(labels)):
                for j in range(4):ax.text(j,i,f'{matrix[i,j]:+.2f}',ha='center',va='center',fontsize=8,color='white' if abs(matrix[i,j])>.42 else '#111')
        fig.colorbar(im,ax=ax,pad=.025,shrink=min(1,5/fig.get_figheight()),label='Signed prediction correlation')
        fig.tight_layout();save(fig,f'atlas_{d}')
    # Same published illustrative outcomes across K; no newly selected best K.
    chosen=[('Waist','__Waist_cm'),('Systolic pressure (BP2, left)','__BP2_left_systole'),('LPS reasoning','__LPS_1'),('CKD-EPI','__CKDEPI_in_ml_min_1.73m'),('HbA1c (%)','__HBA1C_in_%'),('Optimism','__LOT_Optimism')]
    fig,axes=plt.subplots(2,3,figsize=(12,8))
    for ax,(title,suffix) in zip(axes.flat,chosen):
        matches=[r for r in rows if r['label'].endswith(suffix)]
        assert len(matches)==4,(title,len(matches))
        dr={int(r['k_identity_dims']):r for r in matches}
        ax.plot(KS,[float(dr[k]['cv_corr']) for k in KS],'-o',c='#087e8b',label='Leave-one-out')
        ax.plot(KS,[float(dr[k]['in_sample_corr']) for k in KS],'--o',c='#bd7048',label='In sample')
        ax.set(title=title,xticks=KS,xlabel='K',ylabel='Pearson r',ylim=(0,.85));ax.axhline(0,c='#aaa',lw=.5)
    axes[0,0].legend(frameon=False,fontsize=9)
    finish(fig,'How the readouts change with dimensionality','Six previously reported examples · all four K settings · 111 cohort rows with original target imputation')
    save(fig,'phenotype_dimension_curves')
    # New fresh CV metrics vs baseline.
    fig,axes=plt.subplots(1,2,figsize=(12,5.8))
    for ax,target in zip(axes,['LPS','WST']):
        for offset,model,color,label in [(-.13,'score_only','#087e8b','EEG score vs mean'),(.13,'score_plus_covariates','#bd7048','EEG + covariates vs covariates')]:
            for i,name in enumerate(SCORES):
                loo=next(r for r in cv if (r['score'],r['target'],r['model'],r['protocol'])==(name,target,model,'leave_one_out'))
                rep=next(r for r in cv if (r['score'],r['target'],r['model'],r['protocol'])==(name,target,model,'repeated_5fold'))
                y=i+offset
                ax.plot([float(rep['delta_r2_q05']),float(rep['delta_r2_q95'])],[y,y],c=color,lw=4,alpha=.35)
                ax.scatter(float(loo['delta_r2_median']),y,c=color,s=55,label=label if i==0 else None,zorder=4)
        ax.axvline(0,c='#555',lw=1);ax.set(yticks=range(3),yticklabels=list(SCORES.values()),xlabel='Improvement in out-of-fold R² over matched baseline',title=target);ax.invert_yaxis();ax.legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.18))
    finish(fig,'What does the EEG readout add?','New analyses · 111 participants · points: leave-one-out · bars: 5–95% across 20 five-fold partitions (not CIs)')
    save(fig,'cv_incremental_performance')
    # out of fold scatter primary fixed broad score.
    fig,axes=plt.subplots(2,2,figsize=(10,9))
    for ax,(target,model) in zip(axes.flat,[('LPS','score_only'),('LPS','score_plus_covariates'),('WST','score_only'),('WST','score_plus_covariates')]):
        group=[r for r in preds if (r['score'],r['target'],r['model'],r['protocol'])==('theta_concentration_broad',target,model,'leave_one_out')]
        x=np.array([float(r['observed']) for r in group]);y=np.array([float(r['prediction']) for r in group])
        row=next(r for r in cv if (r['score'],r['target'],r['model'],r['protocol'])==('theta_concentration_broad',target,model,'leave_one_out'))
        ax.scatter(x,y,c='#087e8b',s=26,alpha=.7,edgecolors='none')
        bounds=[min(x.min(),y.min()),max(x.max(),y.max())];ax.plot(bounds,bounds,'--',c='#bbb',lw=1)
        ax.set(xlabel=f'Observed {target}',ylabel=f'Out-of-fold predicted {target}',title=f'{target}: '+('EEG only' if model=='score_only' else 'EEG + covariates'))
        ax.text(.03,.97,f"r = {float(row['pearson_r_median']):+.3f}\nR² = {float(row['r2_median']):+.3f}\nBaseline R² = {float(row['baseline_r2_median']):+.3f}",transform=ax.transAxes,va='top',fontsize=9)
    finish(fig,'Participant-level predictions','Broad posterior theta concentration · leave-one-out calibration · identity line shown · scales differ by target')
    save(fig,'cv_predictions')
    # repeated CV distributions for LPS all partitions, each independent whole run summarized once.
    fig,axes=plt.subplots(1,2,figsize=(12,5.6))
    for ax,model in zip(axes,['score_only','score_plus_covariates']):
        for i,(name,label) in enumerate(SCORES.items()):
            group=[r for r in runs if (r['score'],r['target'],r['model'],r['protocol'])==(name,'LPS',model,'repeated_5fold')]
            values=[float(r['r2']) for r in group]
            ax.scatter(np.repeat(i,20)+np.linspace(-.14,.14,20),values,s=30,c=COLORS[i],alpha=.8)
            baseline=np.median([float(r['baseline_r2']) for r in group]);ax.plot([i-.22,i+.22],[baseline,baseline],c='#333',lw=2,label='Median matched baseline' if i==0 else None)
        ax.set(xticks=range(3),xticklabels=['Broad','Lateral','Focus ratio'],ylabel='Out-of-fold R²',title='EEG only' if model=='score_only' else 'EEG + covariates');ax.legend(frameon=False,fontsize=8)
    finish(fig,'Reasoning prediction across 20 data partitions','Five-fold CV · all 111 participants per partition · every run shown · black marks: matched baselines')
    save(fig,'cv_partition_stability')
    outputs=['phenotype_distribution','domain_generalization',*[f'atlas_{d}' for d in DOMAIN],'phenotype_dimension_curves','cv_incremental_performance','cv_predictions','cv_partition_stability']
    (FIG/'figure_manifest.json').write_text(json.dumps({'new_figures':outputs,'count':len(outputs),'formats':['png','svg'],'phenotype_source':'retained aggregate outputs','cv_source':'new participant-level cross-validation','all_phenotype_rows_plotted':len(rows)},indent=2)+'\n')
    print(json.dumps(outputs))

if __name__=='__main__':main()
