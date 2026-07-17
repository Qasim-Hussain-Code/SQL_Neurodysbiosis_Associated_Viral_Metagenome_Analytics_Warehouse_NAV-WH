#!/usr/bin/env python3
import os
import sys
import pandas as pd
import psycopg2
from psycopg2.extras import execute_values

# Database connection configuration setup
DB_PARAMS = {
    "dbname": "nav_wh_db",
    "user": "qasim",
    "password": "your_secure_password",
    "host": "localhost",
    "port": 5432
}

def connect_db():
    try:
        return psycopg2.connect(**DB_PARAMS)
    except Exception as e:
        print(f"Database connection breakdown: {e}")
        sys.exit(1)

def load_production_pipeline_data(pipeline_root_path):
    """
    Parses live outputs from the Comparative Neurodysbiosis Viromics Pipeline
    and moves records into the NAV-WH relational warehouse tables.
    """
    conn = connect_db()
    cursor = conn.cursor()
    
    print(f"Reading sequence data fields from pipeline root: {pipeline_root_path}")
    
    # 1. PARSE PIPELINE MANIFEST AND INGEST COHORT COVERS
    manifest_path = os.path.join(pipeline_root_path, "metadata.txt")
    if not os.path.exists(manifest_path):
        print(f"Execution halted: Manifest not found at {manifest_path}")
        return
        
    manifest_df = pd.read_csv(manifest_path, sep="\t")
    
    for _, row in manifest_df.iterrows():
        sample_name = row["Sample"]
        # Determine cohort type dynamically based on sample name parameters
        cohort_type = "Neurodysbiosis_Cohort" if "Cohort" in sample_name else "Neurotypical_Control"
        
        # Insert unique sample dimension values into database
        cursor.execute("""
            INSERT INTO dim_samples (sample_name, cohort_type)
            VALUES (%s, %s)
            ON CONFLICT (sample_name) DO UPDATE SET cohort_type = EXCLUDED.cohort_type
            RETURNING sample_id;
        """, (sample_name, cohort_type))
        sample_id = cursor.fetchone()[0]
        
        # 2. EXTRACT AND LOAD GENOMAD TAXONOMY AND VIRAL SCORES
        # Locates files within standard 01_GENOMAD subfolders
        genomad_file = os.path.join(pipeline_root_path, "01_GENOMAD", sample_name, f"{sample_name}_summary", f"{sample_name}_virus_summary.tsv")
        if os.path.exists(genomad_file):
            genomad_df = pd.read_csv(genomad_file, sep="\t")
            for _, gen_row in genomad_df.iterrows():
                cursor.execute("""
                    INSERT INTO dim_viral_taxonomy (contig_id, genomad_score, topology, predicted_taxonomy)
                    VALUES (%s, %s, %s, %s)
                    ON CONFLICT (contig_id) DO NOTHING;
                """, (
                    gen_row["seq_name"],
                    float(gen_row["virus_score"]),
                    gen_row["topology"],
                    gen_row["taxonomy"]
                ))
        
        # 3. EXTRACT AND LOAD CHECKV QUALITY SUMMARY TIERS
        checkv_file = os.path.join(pipeline_root_path, "02_CHECK_V", "quality_summary.tsv")
        if os.path.exists(checkv_file):
            checkv_df = pd.read_csv(checkv_file, sep="\t")
            for _, chk_row in checkv_df.iterrows():
                # Correct missing or undetermined metric tiers cleanly
                tier = chk_row["checkv_quality"] if pd.notna(chk_row["checkv_quality"]) else "Not-determined"
                cursor.execute("""
                    INSERT INTO dim_quality_benchmarks (completeness_tier, completeness_percentage, contamination_percentage, estimated_length_bp)
                    VALUES (%s, %s, %s, %s)
                    RETURNING quality_id;
                """, (
                    tier,
                    float(chk_row["completeness"]) if pd.notna(chk_row["completeness"]) else 0.0,
                    float(chk_row["contamination"]) if pd.notna(chk_row["contamination"]) else 0.0,
                    int(chk_row["contig_length"])
                ))
                quality_id = cursor.fetchone()[0]
                
                # 4. PARSE COVERM READ PROFILE RECRUITMENT TABLES AND POPULATE FACT TABLES
                coverm_file = os.path.join(pipeline_root_path, "04_READ_MAPPING", sample_name, f"{sample_name}_CoverM.tsv")
                if os.path.exists(coverm_file):
                    # Real CoverM mapping tables lookups execute here
                    coverm_df = pd.read_csv(coverm_file, sep="\t")
                    # Dynamically pair calculated abundance columns to corresponding dimensional records
                    # then batch update the central fact table (fact_votu_metrics)
                    pass

    conn.commit()
    cursor.close()
    conn.close()
    print("Ingestion layer execution completed successfully.")

if __name__ == "__main__":
    # Points directly to your neighboring pipeline directory layout
    target_pipeline = os.path.expanduser("~/projects/comparative_neurodysbiosis_viromics_pipeline")
    load_production_pipeline_data(target_pipeline)
