#!/usr/bin/env python3
import os
import sys
import pandas as pd
import psycopg2

# Fetch credentials from the local terminal environment
DB_PARAMS = {
    "dbname": os.getenv("NAV_DB_NAME", "nav_wh_db"),
    "user": os.getenv("NAV_DB_USER", "qasim"),
    "password": os.getenv("NAV_DB_PASSWORD", ""),
    "host": os.getenv("NAV_DB_HOST", "localhost"),
    "port": os.getenv("NAV_DB_PORT", "5432")
}

def connect_db():
    try:
        return psycopg2.connect(**DB_PARAMS)
    except Exception as e:
        print(f"Database connection error: {e}")
        print("Tip: Ensure you have exported NAV_DB_PASSWORD or that peer access is allowed.")
        sys.exit(1)

def load_production_pipeline_data(pipeline_root_path):
    conn = connect_db()
    cursor = conn.cursor()
    print(f"Starting ETL extraction from pipeline root: {pipeline_root_path}")
    
    manifest_path = os.path.join(pipeline_root_path, "metadata.txt")
    if not os.path.exists(manifest_path):
        print(f"Aborting: Manifest file not found at {manifest_path}")
        return
        
    manifest_df = pd.read_csv(manifest_path, sep="\t")
    
    for _, row in manifest_df.iterrows():
        sample_name = row["Sample"]
        cohort_type = "Neurodysbiosis_Cohort" if "Cohort" in sample_name else "Neurotypical_Control"
        print(f"Processing sample metrics for: {sample_name} ({cohort_type})")
        
        cursor.execute("""
            INSERT INTO dim_samples (sample_name, cohort_type)
            VALUES (%s, %s)
            ON CONFLICT (sample_name) DO UPDATE SET cohort_type = EXCLUDED.cohort_type
            RETURNING sample_id;
        """, (sample_name, cohort_type))
        sample_id = cursor.fetchone()[0]
        
        genomad_file = os.path.join(pipeline_root_path, "01_GENOMAD", sample_name, f"{sample_name}_summary", f"{sample_name}_virus_summary.tsv")
        if os.path.exists(genomad_file):
            genomad_df = pd.read_csv(genomad_file, sep="\t")
            for _, gen_row in genomad_df.iterrows():
                cursor.execute("""
                    INSERT INTO dim_viral_taxonomy (contig_id, genomad_score, topology, predicted_taxonomy)
                    VALUES (%s, %s, %s, %s)
                    ON CONFLICT (contig_id) DO UPDATE SET genomad_score = EXCLUDED.genomad_score;
                """, (gen_row["seq_name"], float(gen_row["virus_score"]), gen_row["topology"], gen_row["taxonomy"]))
        
        checkv_file = os.path.join(pipeline_root_path, "02_CHECK_V", "quality_summary.tsv")
        quality_map = {}
        if os.path.exists(checkv_file):
            checkv_df = pd.read_csv(checkv_file, sep="\t")
            for _, chk_row in checkv_df.iterrows():
                contig_id = chk_row["contig_id"]
                tier = chk_row["checkv_quality"] if pd.notna(chk_row["checkv_quality"]) else "Not-determined"
                cursor.execute("""
                    INSERT INTO dim_quality_benchmarks (completeness_tier, completeness_percentage, contamination_percentage, estimated_length_bp)
                    VALUES (%s, %s, %s, %s)
                    RETURNING quality_id;
                """, (tier, float(chk_row["completeness"]) if pd.notna(chk_row["completeness"]) else 0.0, float(chk_row["contamination"]) if pd.notna(chk_row["contamination"]) else 0.0, int(chk_row["contig_length"])))
                quality_map[contig_id] = cursor.fetchone()[0]
        
        coverm_file = os.path.join(pipeline_root_path, "04_READ_MAPPING", sample_name, f"{sample_name}_CoverM.tsv")
        if os.path.exists(coverm_file):
            coverm_df = pd.read_csv(coverm_file, sep="\t")
            contig_col = coverm_df.columns[0]
            read_count_col = [c for c in coverm_df.columns if "Read Count" in c or "count" in c.lower()]
            cov_col = [c for c in coverm_df.columns if "Coverage" in c or "covered" in c.lower()]
            rel_abund_col = [c for c in coverm_df.columns if "Abundance" in c or "relative" in c.lower()]
            
            for _, cov_row in coverm_df.iterrows():
                contig_id = cov_row[contig_col]
                cursor.execute("SELECT votu_id FROM dim_viral_taxonomy WHERE contig_id = %s;", (contig_id,))
                votu_res = cursor.fetchone()
                if not votu_res:
                    continue
                votu_id = votu_res[0]
                quality_id = quality_map.get(contig_id, 1)
                
                reads = int(cov_row[read_count_col[0]]) if read_count_col else 0
                coverage = float(cov_row[cov_col[0]]) if cov_col else 0.0
                abundance = float(cov_row[rel_abund_col[0]]) if rel_abund_col else 0.0
                
                cursor.execute("""
                    INSERT INTO fact_votu_metrics (sample_id, votu_id, quality_id, raw_read_count, horizontal_coverage, relative_abundance)
                    VALUES (%s, %s, %s, %s, %s, %s);
                """, (sample_id, votu_id, quality_id, reads, coverage, abundance))
                
    conn.commit()
    cursor.close()
    conn.close()
    print("Warehouse data pipeline synchronization complete.")

if __name__ == "__main__":
    target_pipeline = os.path.expanduser("~/projects/comparative_neurodysbiosis_viromics_pipeline")
    load_production_pipeline_data(target_pipeline)