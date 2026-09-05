"""Render all retained phenotype summary rows; does not refit predictions."""
import csv
import json
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results/context/phenotypes'
NAMES={'blood_chemistry':'Blood chemistry','blood_pressure':'Blood pressure','body_size':'Body measurements','cognition':'Cognition','personality_affect':'Personality and affect'}

def main():
    with (DATA/'domain_label_reconstruction.csv').open(newline='') as f:
        rows=list(csv.DictReader(f))
    with (DATA/'domain_recovery_summary.csv').open(newline='') as f:
        domains=list(csv.DictReader(f))
    grouped=defaultdict(dict)
    for r in rows:
        key=(r['domain'],r['label']); k=int(r['k_identity_dims'])
        assert k not in grouped[key],(key,k)
        grouped[key][k]=r
    assert len(rows)==1164 and len(grouped)==291
    assert all(set(rs)=={3,6,10,20} for rs in grouped.values())
    lines=['# Flint phenotype atlas','',
        '**291 recorded phenotype fields · 1,164 label-by-dimension evaluations · five domains**','',
        'This is the complete retained domain-recovery result table, rendered for browsing. Every label and all four tested coordinate counts are included. Fields can be individual questionnaire items, repeated measurements, or alternate units of the same assay; 291 fields are not 291 independent biological discoveries.','',
        'These values are Pearson correlations between pooled leave-one-out standardized predictions and standardized outcomes. The broader identity-coordinate pipeline refits the EEG reliability basis and phenotype map within each fold. It does not apply the nuisance adjustment used in the primary fluid-reasoning result.','',
        '**Reading the columns:** `Observed n` counts original nonmissing labels in the 111-person cohort. Missing standardized targets are mean-imputed in the original scoring; the correlation therefore includes imputed outcomes. K is the number of retained identity coordinates. A negative prediction correlation is an inverted prediction relationship, not evidence that a particular EEG feature is a protective or adverse biological factor. The values have no per-label significance or multiplicity correction attached.','',
        'This catalog is generated from retained aggregate outputs. Their underlying participant predictions have not been freshly recomputed. The statistical verifier for the primary fluid-reasoning tables covers a separate experiment.','',
        '[Research note](paper.md) · [Full original label table](results/context/phenotypes/domain_label_reconstruction.csv) · [Full original domain table](results/context/phenotypes/domain_recovery_summary.csv)','',
        '## Domain-level performance','',
        'Pooled standardized leave-one-out R² across every label in each domain. These values assess the full domain reconstruction, rather than its most favorable individual row. Negative R² means the prediction errors exceed the globally centered mean baseline used by the original calculation.','',
        '| Domain | Fields | K = 3 | K = 6 | K = 10 | K = 20 |','|---|---:|---:|---:|---:|---:|']
    for d in NAMES:
        rs={int(r['k_identity_dims']):r for r in domains if r['domain']==d}
        assert set(rs)=={3,6,10,20}
        lines.append('| '+NAMES[d]+' | '+rs[3]['n_labels']+' | '+' | '.join(f"{float(rs[k]['loocv_r2']):+.4f}" for k in (3,6,10,20))+' |')
    lines+=['','## All label-level results','',
        'Labels appear in source-name order within each domain. Display names retain the final source group and exact field key; the linked CSV preserves complete original paths. No ordering or filtering uses the correlation magnitude.']
    for d,title in NAMES.items():
        keys=sorted(key for key in grouped if key[0]==d)
        lines+=['',f'### {title} · {len(keys)} fields','',
            '| Source group / field | Observed n | K = 3 | K = 6 | K = 10 | K = 20 |','|---|---:|---:|---:|---:|---:|']
        for key in keys:
            rs=grouped[key]
            assert len({r['n'] for r in rs.values()})==1
            label=' / '.join(key[1].split('__')[-2:]).replace('|','\\|')
            lines.append('| '+label+' | '+rs[3]['n']+' | '+' | '.join(f"{float(rs[k]['cv_corr']):+.4f}" for k in (3,6,10,20))+' |')
    lines+=['','## Rebuild this catalog','','```sh','python scripts/build_phenotype_catalog.py','```','',
        'The builder checks unique label/dimension keys and complete coverage of all 291 fields at K = 3, 6, 10, and 20. It formats the original values without changing their signs or selecting a best K.']
    (ROOT/'PHENOTYPE_ATLAS.md').write_text('\n'.join(lines)+'\n')
    receipt={'operation':'Retained-table catalog generation, not statistical recomputation','source_label_rows':len(rows),'unique_fields':len(grouped),'dimension_settings':[3,6,10,20],'domain_summary_rows':len(domains),'complete_label_dimension_coverage':True,'counts_by_domain':{d:sum(key[0]==d for key in grouped) for d in NAMES}}
    (DATA/'catalog_coverage.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
