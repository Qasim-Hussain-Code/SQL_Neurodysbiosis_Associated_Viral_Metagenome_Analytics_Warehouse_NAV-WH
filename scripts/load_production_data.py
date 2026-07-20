#!/usr/bin/env python3
"""
NAV-WH Production ETL Loader
=============================
Extracts data from the Modular Viromics Pipeline (MVP) output directories
and loads it into the NAV-WH PostgreSQL analytics warehouse.

Source pipeline: comparative_neurodysbiosis_viromics_pipeline
Study: Sharon et al. 2019, Cell 177(3):600-618
Data: ERP113632 (ENA)

Usage:
    python load_production_data.py [--pipeline-root PATH] [--dry-run]
"""

import os
import sys
import json
import argparse
import pandas as pd

try:
    import psycopg2
    from psycopg2.extras import Json
except ImportError:
    print("Error: psycopg2 not installed. Run: pip install psycopg2-binary")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Database connection
# ---------------------------------------------------------------------------
DB_PARAMS = {
    "dbname":   os.getenv("NAV_DB_NAME",     "nav_wh_db"),
    "user":     os.getenv("NAV_DB_USER",     "qasim"),
    "password": os.getenv("NAV_DB_PASSWORD", ""),
    "host":     os.getenv("NAV_DB_HOST",     "localhost"),
    "port":     os.getenv("NAV_DB_PORT",     "5432"),
}

# ENA accession mapping for the Sharon et al. 2019 study
ACCESSION_MAP = {
    "Neurodysbiosis_Cohort":  "ERR3144321",
    "Neurotypical_Control":   "ERR3144319",
}


def connect_db():
    """Establish a database connection."""
    try:
        conn = psycopg2.connect(**DB_PARAMS)
        conn.autocommit = False
        return conn
    except Exception as e:
        print(f"[ERROR] Database connection failed: {e}")
        print("Tip: Ensure NAV_DB_PASSWORD is set or peer authentication is enabled.")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Phase 1: Load dim_samples from metadata.txt
# ---------------------------------------------------------------------------
def load_samples(cursor, pipeline_root):
    """Load sample dimension from the pipeline metadata.txt manifest."""
    manifest = os.path.join(pipeline_root, "metadata.txt")
    if not os.path.exists(manifest):
        print(f"[ERROR] Manifest not found: {manifest}")
        return {}

    df = pd.read_csv(manifest, sep="\t")
    sample_id_map = {}

    for _, row in df.iterrows():
        sample_name = row["Sample"]
        cohort_type = sample_name  # column already contains the cohort label
        accession = ACCESSION_MAP.get(sample_name, None)

        cursor.execute("""
            INSERT INTO dim_samples (sample_name, cohort_type, read_accession)
            VALUES (%s, %s, %s)
            ON CONFLICT (sample_name)
            DO UPDATE SET cohort_type = EXCLUDED.cohort_type,
                          read_accession = EXCLUDED.read_accession
            RETURNING sample_id;
        """, (sample_name, cohort_type, accession))

        sample_id_map[sample_name] = cursor.fetchone()[0]
        print(f"  [dim_samples] {sample_name} -> sample_id={sample_id_map[sample_name]}")

    return sample_id_map


# ---------------------------------------------------------------------------
# Phase 2: Load dim_quality_benchmarks + dim_viral_taxonomy from Module 03
# ---------------------------------------------------------------------------
def load_votu_quality_and_taxonomy(cursor, pipeline_root):
    """
    Load quality benchmarks and viral taxonomy dimensions from the
    Module 03 representative vOTU quality summary table. This is the
    single authoritative source for vOTU identity after clustering.
    """
    votu_file = os.path.join(
        pipeline_root, "03_CLUSTERING",
        "MVP_03_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_"
        "Representative_Virus_Proviruses_Quality_Summary.tsv"
    )

    if not os.path.exists(votu_file):
        print(f"[ERROR] Representative vOTU file not found: {votu_file}")
        return {}, {}

    df = pd.read_csv(votu_file, sep="\t")
    quality_id_map = {}   # contig_id -> quality_id
    votu_id_map = {}      # contig_id -> votu_id

    for _, row in df.iterrows():
        contig_id   = str(row["Representative_Sequence"])
        members     = str(row["Sequences"])
        sample      = str(row["Sample"])
        length      = int(row["virus_length"])
        provirus    = str(row["provirus"]).strip().lower() == "yes"
        gene_count  = int(row["gene_count"])       if pd.notna(row["gene_count"])  else 0
        viral_genes = int(row["viral_genes"])      if pd.notna(row["viral_genes"]) else 0
        host_genes  = int(row["host_genes"])       if pd.notna(row["host_genes"])  else 0
        quality     = str(row["checkv_quality"])    if pd.notna(row["checkv_quality"])    else "Not-determined"
        miuvig      = str(row.get("miuvig_quality", ""))  if pd.notna(row.get("miuvig_quality")) else None
        completeness = float(row["completeness"])  if pd.notna(row["completeness"]) else 0.0
        comp_method = str(row["completeness_method"]) if pd.notna(row["completeness_method"]) else None
        kmer_freq   = float(row["kmer_freq"])      if pd.notna(row["kmer_freq"])   else None
        coordinates = str(row["coordinates"])      if pd.notna(row["coordinates"]) and str(row["coordinates"]).strip() else None

        virus_score = float(row["virus_score"])    if pd.notna(row["virus_score"]) else 0.0
        genetic_code = int(row["genetic_code"])    if pd.notna(row["genetic_code"]) else 11
        n_hallmarks = int(row["n_hallmarks"])      if pd.notna(row["n_hallmarks"]) else 0
        marker_enr  = float(row["marker_enrichment"]) if pd.notna(row["marker_enrichment"]) else None
        taxonomy    = str(row["taxonomy"])          if pd.notna(row["taxonomy"])    else "Unclassified"
        genome_type = str(row["Genome type"])       if pd.notna(row.get("Genome type")) else "Unknown"
        host_type   = str(row["Host type"])         if pd.notna(row.get("Host type"))   else "Unknown"

        # --- Insert quality benchmark ---
        cursor.execute("""
            INSERT INTO dim_quality_benchmarks
                (contig_id, checkv_quality, miuvig_quality, completeness,
                 completeness_method, virus_length_bp, provirus,
                 provirus_coordinates, gene_count, viral_genes,
                 host_genes, kmer_freq)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (contig_id)
            DO UPDATE SET checkv_quality = EXCLUDED.checkv_quality,
                          completeness   = EXCLUDED.completeness
            RETURNING quality_id;
        """, (contig_id, quality, miuvig, completeness, comp_method,
              length, provirus, coordinates, gene_count, viral_genes,
              host_genes, kmer_freq))

        quality_id_map[contig_id] = cursor.fetchone()[0]

        # --- Insert viral taxonomy ---
        cursor.execute("""
            INSERT INTO dim_viral_taxonomy
                (contig_id, cluster_members, representative_sample,
                 virus_score, genetic_code, n_hallmarks,
                 marker_enrichment, predicted_taxonomy,
                 genome_type, host_type)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (contig_id)
            DO UPDATE SET virus_score       = EXCLUDED.virus_score,
                          predicted_taxonomy = EXCLUDED.predicted_taxonomy
            RETURNING votu_id;
        """, (contig_id, members, sample, virus_score, genetic_code,
              n_hallmarks, marker_enr, taxonomy, genome_type, host_type))

        votu_id_map[contig_id] = cursor.fetchone()[0]

    print(f"  [dim_quality_benchmarks] Loaded {len(quality_id_map)} vOTU quality records")
    print(f"  [dim_viral_taxonomy]     Loaded {len(votu_id_map)} vOTU taxonomy records")
    return quality_id_map, votu_id_map


# ---------------------------------------------------------------------------
# Phase 3: Load fact_votu_metrics from Module 04 CoverM output
# ---------------------------------------------------------------------------
def load_read_mapping_metrics(cursor, pipeline_root, sample_id_map, quality_id_map, votu_id_map):
    """
    Load read-recruitment quantitative metrics from CoverM TSV output
    for each sample in the pipeline.
    """
    total_loaded = 0

    for sample_name, sample_id in sample_id_map.items():
        coverm_file = os.path.join(
            pipeline_root, "04_READ_MAPPING", sample_name,
            f"{sample_name}_CoverM.tsv"
        )

        if not os.path.exists(coverm_file):
            print(f"  [WARN] CoverM file not found for {sample_name}: {coverm_file}")
            continue

        df = pd.read_csv(coverm_file, sep="\t")
        sample_loaded = 0

        for _, row in df.iterrows():
            contig_id = str(row["Contig"])

            if contig_id not in votu_id_map:
                continue

            votu_id    = votu_id_map[contig_id]
            quality_id = quality_id_map.get(contig_id, None)

            if quality_id is None:
                continue

            cursor.execute("""
                INSERT INTO fact_votu_metrics
                    (sample_id, votu_id, quality_id,
                     mean_depth, trimmed_mean_depth, covered_bases,
                     covered_fraction, coverage_variance, contig_length,
                     raw_read_count, reads_per_base, rpkm, tpm)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (sample_id, votu_id)
                DO UPDATE SET mean_depth       = EXCLUDED.mean_depth,
                              rpkm             = EXCLUDED.rpkm,
                              raw_read_count   = EXCLUDED.raw_read_count,
                              covered_fraction = EXCLUDED.covered_fraction;
            """, (
                sample_id, votu_id, quality_id,
                float(row["Mean"]),
                float(row["Trimmed Mean"]),
                int(row["Covered Bases"]),
                float(row["Covered Fraction"]),
                float(row["Variance"]),
                int(row["Length"]),
                int(row["Read Count"]),
                float(row["Reads per base"]),
                float(row["RPKM"]),
                float(row["TPM"]),
            ))
            sample_loaded += 1

        print(f"  [fact_votu_metrics] {sample_name}: {sample_loaded} vOTU metrics loaded")
        total_loaded += sample_loaded

    return total_loaded


# ---------------------------------------------------------------------------
# Phase 4: Load dim_gene_annotations from Module 06 geNomad annotation
# ---------------------------------------------------------------------------
def load_gene_annotations(cursor, pipeline_root, votu_id_map):
    """
    Load functional gene annotations from the geNomad annotation table
    for representative vOTU sequences.
    """
    annotation_file = os.path.join(
        pipeline_root, "06_FUNCTIONAL_ANNOTATION",
        "MVP_06_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_"
        "Representative_Virus_Proviruses_Gene_Annotation_GENOMAD.tsv"
    )

    if not os.path.exists(annotation_file):
        print(f"  [WARN] Annotation file not found: {annotation_file}")
        return 0

    df = pd.read_csv(annotation_file, sep="\t")
    loaded = 0

    for _, row in df.iterrows():
        contig_id = str(row["Contig_name"]) if pd.notna(row.get("Contig_name")) else None
        if contig_id is None or contig_id not in votu_id_map:
            continue

        gene_id_str = str(row["Viral_gene_ID"]) if pd.notna(row.get("Viral_gene_ID")) else None
        if gene_id_str is None:
            continue

        cursor.execute("""
            INSERT INTO dim_gene_annotations
                (viral_gene_id, contig_id, gene_number,
                 gene_start, gene_end, gene_length, strand,
                 genomad_annotation, genomad_accession, genomad_score,
                 genomad_evalue, gc_content, rbs_motif,
                 genomad_marker, virus_hallmark, amr_annotation)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
        """, (
            gene_id_str,
            contig_id,
            int(row["Gene_number"])    if pd.notna(row.get("Gene_number"))    else None,
            int(row["start"])          if pd.notna(row.get("start"))          else None,
            int(row["end"])            if pd.notna(row.get("end"))            else None,
            int(row["length"])         if pd.notna(row.get("length"))         else None,
            int(row["strand"])         if pd.notna(row.get("strand"))         else None,
            str(row["GENOMAD_Annotation"])          if pd.notna(row.get("GENOMAD_Annotation"))          else "Unknown",
            str(row["GENOMAD_Annotation_accessions"]) if pd.notna(row.get("GENOMAD_Annotation_accessions")) else None,
            float(row["GENOMAD_Score"])              if pd.notna(row.get("GENOMAD_Score"))              else None,
            float(row["GENOMAD_evalue"])             if pd.notna(row.get("GENOMAD_evalue"))             else None,
            float(row["GENOMAD_gc_content"])         if pd.notna(row.get("GENOMAD_gc_content"))         else None,
            str(row["GENOMAD_rbs_motif"])            if pd.notna(row.get("GENOMAD_rbs_motif"))          else None,
            str(row["GENOMAD_marker"])               if pd.notna(row.get("GENOMAD_marker"))             else None,
            str(row["GENOMAD_virus_hallmark"])       if pd.notna(row.get("GENOMAD_virus_hallmark"))     else None,
            str(row["GENOMAD_Annotation_amr"])       if pd.notna(row.get("GENOMAD_Annotation_amr"))     else None,
        ))
        loaded += 1

    print(f"  [dim_gene_annotations] Loaded {loaded} gene annotation records")
    return loaded


# ---------------------------------------------------------------------------
# Main ETL orchestrator
# ---------------------------------------------------------------------------
def run_etl(pipeline_root, dry_run=False):
    """Execute the full ETL pipeline."""
    print("=" * 72)
    print(" NAV-WH Production ETL Loader")
    print(f" Pipeline root: {pipeline_root}")
    print(f" Dry run: {dry_run}")
    print("=" * 72)

    if not os.path.isdir(pipeline_root):
        print(f"[ERROR] Pipeline root does not exist: {pipeline_root}")
        sys.exit(1)

    conn = connect_db()
    cursor = conn.cursor()

    try:
        print("\n[Phase 1] Loading sample dimension...")
        sample_id_map = load_samples(cursor, pipeline_root)

        print("\n[Phase 2] Loading vOTU quality + taxonomy dimensions...")
        quality_id_map, votu_id_map = load_votu_quality_and_taxonomy(cursor, pipeline_root)

        print("\n[Phase 3] Loading read-recruitment fact metrics...")
        n_metrics = load_read_mapping_metrics(
            cursor, pipeline_root, sample_id_map, quality_id_map, votu_id_map
        )

        print("\n[Phase 4] Loading functional gene annotations...")
        n_genes = load_gene_annotations(cursor, pipeline_root, votu_id_map)

        if dry_run:
            conn.rollback()
            print("\n[DRY RUN] All changes rolled back.")
        else:
            conn.commit()
            print(f"\n[SUCCESS] ETL complete:")
            print(f"  Samples:     {len(sample_id_map)}")
            print(f"  vOTUs:       {len(votu_id_map)}")
            print(f"  Metrics:     {n_metrics}")
            print(f"  Annotations: {n_genes}")

    except Exception as e:
        conn.rollback()
        print(f"\n[ERROR] ETL failed, transaction rolled back: {e}")
        raise
    finally:
        cursor.close()
        conn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NAV-WH Production ETL Loader")
    parser.add_argument(
        "--pipeline-root",
        default=os.path.expanduser(
            "~/projects/comparative_neurodysbiosis_viromics_pipeline"
        ),
        help="Path to the MVP pipeline output root directory",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run without committing changes (rolls back at the end)",
    )
    args = parser.parse_args()
    run_etl(args.pipeline_root, args.dry_run)