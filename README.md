# Alzheimer's Disease Drug Discovery Pipeline

## Overview

A computational pipeline for Alzheimer's disease drug discovery using
multi-omics integration and drug repurposing approaches. Current work
uses the GSE44770 dataset (230 prefrontal cortex samples, 129 AD / 101
control) to identify differentially expressed genes and enriched
pathways as a foundation for downstream drug repurposing candidate
prioritization.

## Status

This is an active, in-progress research project. The differential
expression and pathway enrichment stage (below) is complete and
reproducible. Multi-omics integration, drug repurposing screening, and
ML-based candidate prioritization are in development.

## Pipeline Stages

1. Data ingestion & preprocessing — done
2. Differential expression analysis (limma) — done
3. Gene annotation & KEGG pathway enrichment — done
4. STRING PPI network & hub gene identification — done
5. Multi-omics integration — in progress
6. Drug repurposing screen & ML candidate prioritization — in progress

## Tech Stack

- R 4.6.0 (GEOquery, limma, clusterProfiler, STRINGdb, pheatmap)
- Python 3.10+ (planned: pandas, numpy, scikit-learn, rdkit)
- Git for version control

## Project Structure

- `scripts/` → R analysis pipeline (data acquisition through PPI network analysis)
- `results/` → pipeline output figures (e.g. top-50 DEG heatmap)
- `docs/` → mini-project writeup with methods, results, and discussion

## Key Result So Far

Differential expression analysis identified 542 significantly
differentially expressed genes (|log2FC| > 0.3, adj. p < 0.05) between
AD and control samples, with 477 annotated to Entrez gene IDs. KEGG
enrichment of these DEGs recovered pathways consistent with known AD
biology, including complement and coagulation cascades and GABAergic
synaptic signaling. See `docs/GSE44770_PhD_MiniProject.md` for full
methods and results, and `results/GSE44770_heatmap.png` for
unsupervised clustering of the top 50 DEGs.

## Author

Olamiposi Owolabi | Computational Biology | Alzheimer's Drug Discovery

## Status: AI-assisted, under audit
Parts of this pipeline were written with AI tools and were not independently checked at first. An audit is under way; see `docs/01_data_audit.md` for what has been verified, what was found wrong, and what is still pending. Results in `docs/GSE44770_PhD_MiniProject.md` come from an unadjusted model and should not be cited as findings.
