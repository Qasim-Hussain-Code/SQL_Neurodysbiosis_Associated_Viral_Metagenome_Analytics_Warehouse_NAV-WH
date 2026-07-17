-- ============================================================================
-- PROJECT: Neurodysbiosis-Associated Viral Metagenome Analytics Warehouse (NAV-WH)
-- ARCHITECTURE: Star Schema optimized for comparative read-recruitment profiling
-- COMPATIBILITY: PostgreSQL 14+ 
-- ============================================================================

CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- 1. DIMENSION: SAMPLES & COHORT METADATA
CREATE TABLE dim_samples (
    sample_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    sample_uuid UUID DEFAULT uuid_generate_v4() UNIQUE,
    sample_name VARCHAR(100) NOT NULL UNIQUE,
    cohort_type VARCHAR(50) NOT NULL, -- 'Neurodysbiosis_Cohort' or 'Neurotypical_Control'
    sequencing_platform VARCHAR(50) DEFAULT 'Illumina NovaSeq',
    extraction_kit VARCHAR(100),
    geographic_location VARCHAR(100),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    
    CONSTRAINT chk_cohort_type CHECK (cohort_type IN ('Neurodysbiosis_Cohort', 'Neurotypical_Control'))
);

-- 2. DIMENSION: CHECKV QUALITY BENCHMARKS
CREATE TABLE dim_quality_benchmarks (
    quality_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    checkv_version VARCHAR(20) DEFAULT 'v1.0.1',
    completeness_tier VARCHAR(50) NOT NULL, -- 'Complete', 'High-quality', etc.
    completeness_percentage NUMERIC(5, 2) NOT NULL,
    contamination_percentage NUMERIC(5, 2) NOT NULL,
    estimated_length_bp INT NOT NULL,
    provirus_flag BOOLEAN DEFAULT FALSE,
    host_genes_count INT DEFAULT 0,
    viral_genes_count INT DEFAULT 0,
    
    CONSTRAINT chk_completeness_range CHECK (completeness_percentage BETWEEN 0.00 AND 100.00),
    CONSTRAINT chk_contamination_range CHECK (contamination_percentage BETWEEN 0.00 AND 100.00),
    CONSTRAINT chk_quality_tier CHECK (completeness_tier IN ('Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined'))
);

-- 3. DIMENSION: VIRAL TAXONOMY & GENOMAD CLASSIFICATION
CREATE TABLE dim_viral_taxonomy (
    votu_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    contig_id VARCHAR(255) NOT NULL UNIQUE,
    genomad_score NUMERIC(4, 3) NOT NULL,
    topology VARCHAR(20) DEFAULT 'Linear',
    predicted_taxonomy VARCHAR(255) DEFAULT 'Unclassified Viridiplantae/Dark_Matter',
    genetic_code INT DEFAULT 11,
    strand_orientation CHAR(1) DEFAULT '+',
    raw_genomad_annotations JSONB,
    
    CONSTRAINT chk_genomad_score CHECK (genomad_score BETWEEN 0.000 AND 1.000),
    CONSTRAINT chk_strand CHECK (strand_orientation IN ('+', '-', '*'))
);

-- 4. FACT TABLE: HIGH-THROUGHPUT vOTU QUANTITATIVE METRICS
CREATE TABLE fact_votu_metrics (
    metric_id BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    sample_id INT NOT NULL,
    votu_id INT NOT NULL,
    quality_id INT NOT NULL,
    raw_read_count INT NOT NULL DEFAULT 0,
    mapped_read_depth NUMERIC(12, 2) NOT NULL DEFAULT 0.00,
    horizontal_coverage NUMERIC(5, 2) NOT NULL DEFAULT 0.00,
    relative_abundance DOUBLE PRECISION NOT NULL DEFAULT 0.000000,
    
    CONSTRAINT fk_fact_sample FOREIGN KEY (sample_id) REFERENCES dim_samples(sample_id) ON DELETE RESTRICT,
    CONSTRAINT fk_fact_votu FOREIGN KEY (votu_id) REFERENCES dim_viral_taxonomy(votu_id) ON DELETE RESTRICT,
    CONSTRAINT fk_fact_quality FOREIGN KEY (quality_id) REFERENCES dim_quality_benchmarks(quality_id) ON DELETE RESTRICT,
        
    CONSTRAINT chk_horizontal_cov CHECK (horizontal_coverage BETWEEN 0.00 AND 100.00),
    CONSTRAINT chk_raw_reads CHECK (raw_read_count >= 0)
);

-- 5. PERFORMANCE OPTIMIZATION LAYER (INDEXES)
CREATE INDEX idx_fact_sample_votu ON fact_votu_metrics (sample_id, votu_id);
CREATE INDEX idx_fact_high_coverage ON fact_votu_metrics (horizontal_coverage) WHERE horizontal_coverage >= 75.00;
CREATE INDEX idx_gin_genomad_annotations ON dim_viral_taxonomy USING gin (raw_genomad_annotations);
