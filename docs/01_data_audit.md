# GSE44770 audit

Every number below is copied from console output saved in docs/ and results/. Anything not under "Corrected results" is pending.

Environment: WSL Ubuntu; R 4.6.0, GEOquery 2.80.1, limma 3.68.5; Python 3.13.13, scikit-learn 1.9.1, numpy 2.5.4, pandas 3.1.1.

## Data facts
- GEO series GSE44770, platform GPL4372: 39,280 probes x 230 samples (129 AD, 101 control), dorsolateral prefrontal cortex, Harvard Brain Tissue Resource Center. GEO lists PubMed ID 23622250 (verify it is the source paper).
- Two-colour custom Agilent 44K array processed by Rosetta (sample header text). Channel 1 (Cy3) is a common pool (source_name_ch1 = HBTRC_PF_Pool for all 230 samples); channel 2 (Cy5) is the sample.
- GEO column description for VALUE: normalized log10 ratio (Cy5/Cy3), test/reference. So effect sizes are log10 ratios; a cutoff of 0.3 is a 2-fold change. Processing text: Rosetta Resolver error-weighted summarisation. Any further normalization is not stated in the header.
- Value range -2.000 to 2.000 (2,262 values at |value| >= 1.999), median 0.005, 342 NA.
- Covariates: age (age:ch2; age:ch1 is NA for all samples because channel 1 is the pool), sex, PMI, pH, RIN, batch (1/2/3), preservation (dry ice / LNV).

| | AD (n=129) | Control (n=101) |
|---|---|---|
| Mean age | 80.2 | 62.1 |
| Mean PMI (h) | 13.7 | 22.4 |
| Mean pH | 6.37 | 6.55 |
| Mean RIN | 7.12 | 7.31 |
| Female | 67 | 19 |
| Batch 1/2/3 | 91/28/10 | 81/17/3 |
| Dry ice / LNV | 28/101 | 7/94 |

PC1 (36.5% of variance) is associated with diagnosis (p=1.5e-29), age (2e-25), RIN (1.3e-16), PMI, pH and sex (p 3e-11 to 3e-9). PC2 (9.2%) tracks batch (p=5e-28).

## Errors found
1. The original DE model (scripts/GSE44770_analysis_pipeline.R) was `~0 + group` with no covariates, despite the group differences above.
2. Effect sizes were labelled log2FC in the writeup; they are log10 ratios. The writeup also says age was NA for all samples; age:ch2 is complete.
3. In the old results, 11 probes for XIST, RPS4Y1, DDX3Y, KDM5D, UTY and EIF1AY (matched by Entrez ID) were significant (FDR <= 0.038). Five of them (XIST x2, RPS4Y1, DDX3Y, KDM5D) passed the |logFC| > 0.3 cutoff and were in the 542. This reflects the sex imbalance (52% female in AD, 19% in controls). Only these six genes were checked.
4. 26,145 of 39,280 probes had FDR < 0.05 in the old model, so the 542-probe list was set by the |logFC| > 0.3 cutoff, not by statistics.
5. Classifier (AD_classifier.py): its 200 input probes appear to have been chosen with the labels on all 230 samples before the CV (overlap with the top 200 by |logFC| in the old DE: 92; no code on this machine builds AD_expression_for_ML.csv, so the rule is not confirmed). The ROC curve and its AUC (section 7 of the script) are scored on the training data. No covariate control.
6. Writeup inconsistencies: 542 vs 477 DEG counts, a duplicated figure under Figure 4, PPI counts (451 vs 394), Table 2 ordering, GEOquery version, one reference to verify.
- Not found: voom does not appear in either local git repository or in scripts/GSE44770_analysis_pipeline.R, which uses lmFit/eBayes (other local files not searched). The writeup already names the platform as Rosetta/Merck.

## Changes made
- scripts/02_corrected_de.R: limma with age, sex, PMI, pH, RIN, batch. Old results kept in old_results_v1/ (CSVs are untracked via .gitignore).
- scripts/03_export_for_ml.R, scripts/04_cv_leakage_free.py: CV with feature selection (F-test, 200 probes) inside each training fold.

## Corrected results
- Unadjusted rerun reproduces the old results exactly (max |logFC| difference 1.3e-15): 542 DEGs at FDR < 0.05 and |log10FC| > 0.3.
- Adjusted model: 16,870 probes at FDR < 0.05; 71 at |log10FC| > 0.3 (all within the old 542); 643 at |log10FC| > 0.176; 66 in an age >= 55 sensitivity run. All 50 old top probes stay at FDR < 0.05 (median adjusted/old logFC ratio 0.65). The 11 sex-linked probes become null (FDR 0.25 to 0.94).
- Classifier, repeated 5-fold CV x 10, ROC AUC (SD across folds is not a confidence interval): expression only 0.962 (SD 0.024); covariates only 0.921 (0.037); expression + covariates 0.963 (0.025); old design (selection before CV) 0.967 (0.022). Single cohort, no external validation. The models differ (random forest vs logistic regression). Covariates alone reach 0.92, so the expression AUC cannot be taken as evidence of disease-specific signal; this design cannot separate disease from confounding.

## Still pending
- KEGG, PPI and CMap results were built on the unadjusted DEG list and are unvalidated.
- Classifier on covariate-adjusted expression; external validation; permutation test.
- Preservation (dry ice / LNV) sensitivity model; correcting the writeup.
