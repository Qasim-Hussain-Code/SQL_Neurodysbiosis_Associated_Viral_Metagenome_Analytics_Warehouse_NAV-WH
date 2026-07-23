[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./license)

# Neurodysbiosis-Associated Viral Metagenome Analytics Warehouse (NAV-WH)

A PostgreSQL star-schema analytics warehouse designed to store, query, and compare viral metagenome operational taxonomic unit (vOTU) read-recruitment metrics from comparative neurodysbiosis viromics experiments.

## Context

This warehouse is designed to ingest data from the **Modular Viromics Pipeline (MVP)** applied to whole metagenome sequencing data from [Sharon et al. 2019](https://doi.org/10.1016/j.cell.2019.05.004), *Cell* 177(3):600-618.

**Study design:** Fecal microbiota from human donors with Autism Spectrum Disorder (ASD) and typically developing (TD) controls were transplanted into germ-free mice. The gut viral metagenome (virome) of colonized mice is profiled to identify condition-associated bacteriophage populations.

| Sample | Cohort | ENA Accession | Donor Type |
|--------|--------|---------------|------------|
| Neurodysbiosis_Cohort | ASD donor mouse | ERR3144321 | ASD |
| Neurotypical_Control | TD donor mouse | ERR3144319 | TD |

## Architecture

```
                    ┌─────────────────────┐
                    │   dim_samples       │
                    │   (cohort metadata) │
                    └────────┬────────────┘
                             │
┌──────────────────┐    ┌────┴───────────────┐    ┌──────────────────────┐
│ dim_quality_     │    │ fact_votu_metrics   │    │ dim_viral_taxonomy   │
│ benchmarks       ├────┤ (CoverM coverage,  ├────┤ (geNomad scores,     │
│ (CheckV quality, │    │  RPKM, TPM, reads) │    │  taxonomy, hallmarks)│
│  completeness)   │    └────────────────────┘    └──────────┬───────────┘
└──────────────────┘                                         │
                                                  ┌──────────┴───────────┐
                                                  │ dim_gene_annotations │
                                                  │ (geNomad/PHROGs/PFAM │
                                                  │  functional genes)   │
                                                  └──────────────────────┘
```

### Tables

| Table | Type | Description |
|-------|------|-------------|
| `dim_samples` | Dimension | Sample identity, cohort assignment, ENA accession |
| `dim_quality_benchmarks` | Dimension | CheckV quality tiers, completeness, gene counts |
| `dim_viral_taxonomy` | Dimension | geNomad virus scores, ICTV taxonomy, genome/host type |
| `dim_gene_annotations` | Dimension | Per-gene functional annotations from geNomad |
| `fact_votu_metrics` | Fact | CoverM read-recruitment metrics (depth, coverage, RPKM, TPM) |
| `mv_high_confidence_votus` | Materialized View | Pre-filtered publication-quality vOTUs |

## Directory Structure

```
NAV-WH/
├── migrations/
│   └── 01_db_init_migration.sql     # DDL: tables, indexes, constraints
├── queries/
│   ├── 01_cohort_abundance_ranking.sql   # Top-N vOTUs per cohort by RPKM
│   ├── 02_json_annotation_parsing.sql    # Functional gene summaries + hallmarks
│   └── 03_high_confidence_votus_view.sql # Mat. view + differential abundance
├── scripts/
│   ├── load_production_data.py      # ETL loader (MVP pipeline → PostgreSQL)
│   └── seed_mock_data.py            # Mock data generator (TBD)
├── benchmarks/
│   └── query_performance_log.txt    # EXPLAIN ANALYZE results
├── .env                             # Database credentials (gitignored)
├── .gitignore
└── README.md
```

## Setup

### Prerequisites

- PostgreSQL 14+
- Python 3.8+ with `pandas` and `psycopg2-binary`
- Completed MVP pipeline run in `~/projects/comparative_neurodysbiosis_viromics_pipeline/`

### 1. Create the database

```bash
createdb nav_wh_db
```

### 2. Run the migration

```bash
psql -d nav_wh_db -f migrations/01_db_init_migration.sql
```

### 3. Set credentials

```bash
export NAV_DB_PASSWORD="your_secure_password"
# Or edit .env and source it:
source .env
```

### 4. Load pipeline data

```bash
# Dry run (no changes committed):
python scripts/load_production_data.py --dry-run

# Production load:
python scripts/load_production_data.py
```

## Analytical Queries

### Query 01: Cohort Abundance Ranking

Returns the top-10 most abundant vOTUs per cohort, filtered by ≥75% genome coverage:

```bash
psql -d nav_wh_db -f queries/01_cohort_abundance_ranking.sql
```

### Query 02: Functional Gene Analysis

Summarizes geNomad functional annotations with virus hallmark gene detection:

```bash
psql -d nav_wh_db -f queries/02_json_annotation_parsing.sql
```

### Query 03: Differential Abundance + High-Confidence View

Creates a materialized view of publication-quality vOTUs and compares RPKM between cohorts with enrichment classification:

```bash
psql -d nav_wh_db -f queries/03_high_confidence_votus_view.sql
```

## Data Sources (MVP Pipeline Mapping)

| Warehouse Component | MVP Source File |
|---------------------|----------------|
| `dim_samples` | `metadata.txt` |
| `dim_quality_benchmarks` | `03_CLUSTERING/MVP_03_*_Representative_*_Quality_Summary.tsv` |
| `dim_viral_taxonomy` | `03_CLUSTERING/MVP_03_*_Representative_*_Quality_Summary.tsv` |
| `fact_votu_metrics` | `04_READ_MAPPING/{sample}/{sample}_CoverM.tsv` |
| `dim_gene_annotations` | `06_FUNCTIONAL_ANNOTATION/MVP_06_*_Gene_Annotation_GENOMAD.tsv` |

## License

This project is licensed under the [MIT License](./license). You are free to use, modify, and distribute this software for academic, research, or commercial purposes. See the [license](./license) file for full terms.
