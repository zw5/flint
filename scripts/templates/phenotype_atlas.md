# Flint phenotype atlas

**291 recorded phenotype fields · 1,164 label-by-dimension evaluations · five domains**

This is the complete retained domain-recovery result table, rendered for browsing. Every label and all four tested coordinate counts are included. Fields can be individual questionnaire items, repeated measurements, or alternate units of the same assay; 291 fields are not 291 independent biological discoveries.

These values are Pearson correlations between pooled leave-one-out standardized predictions and standardized outcomes. The broader identity-coordinate pipeline refits the EEG reliability basis and phenotype map within each fold. It does not apply the nuisance adjustment used in the primary fluid-reasoning result.

**Reading the columns:** `Observed n` counts original nonmissing labels in the 111-person cohort. Missing standardized targets are mean-imputed in the original scoring; the correlation therefore includes imputed outcomes. K is the number of retained identity coordinates. A negative prediction correlation is an inverted prediction relationship, not evidence that a particular EEG feature is a protective or adverse biological factor. The values have no per-label significance or multiplicity correction attached.

This catalog is generated from retained aggregate outputs. Their underlying participant predictions have not been freshly recomputed. The statistical verifier for the primary fluid-reasoning tables covers a separate experiment.

[Research note](paper.md) · [Full original label table](results/context/phenotypes/domain_label_reconstruction.csv) · [Full original domain table](results/context/phenotypes/domain_recovery_summary.csv)

## Domain-level performance

Pooled standardized leave-one-out R² across every label in each domain. These values assess the full domain reconstruction, rather than its most favorable individual row. Negative R² means the prediction errors exceed the globally centered mean baseline used by the original calculation.

${domain_table}

## All label-level results

Labels appear in source-name order within each domain. Display names retain the final source group and exact field key; the linked CSV preserves complete original paths. No ordering or filtering uses the correlation magnitude.

${label_tables}

## Rebuild this catalog

```sh
python scripts/build_phenotype_catalog.py
```

The builder checks unique label/dimension keys and complete coverage of all 291 fields at K = 3, 6, 10, and 20. It formats the original values without changing their signs or selecting a best K.
