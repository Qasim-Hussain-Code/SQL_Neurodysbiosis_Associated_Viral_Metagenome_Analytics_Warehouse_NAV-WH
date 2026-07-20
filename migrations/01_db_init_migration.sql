-- ============================================================================
-- PROJECT: Neurodysbiosis-Associated Viral Metagenome Analytics Warehouse (NAV-WH)
-- SOURCE:  Modular Viromics Pipeline (MVP) applied to Sharon et al. 2019 data
-- ARCHITECTURE: Star Schema optimized for comparative viromics profiling
-- COMPATIBILITY: PostgreSQL 14+
-- ============================================================================

CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- ============================================================================
-- 1. DIMENSION: SAMPLES & COHORT METADATA
-- Source: metadata.txt from the viromics pipeline
-- ============================================================================
CREATE TABLE dim_samples (
    sample_id       INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    sample_uuid     UUID DEFAULT uuid_generate_v4() UNIQUE,
    sample_name     VARCHAR(100) NOT NULL UNIQUE,
    cohort_type     VARCHAR(50)  NOT NULL,
    read_accession  VARCHAR(20),                              -- ENA accession (e.g. ERR3144321)
    sequencing_platform VARCHAR(50) DEFAULT 'Illumina HiSeq',
    study_accession VARCHAR(20)  DEFAULT 'ERP113632',         -- Sharon et al. 2019
    created_at      TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT chk_cohort_type CHECK (
        cohort_type IN ('Neurodysbiosis_Cohort', 'Neurotypical_Control')
    )
);

-- ============================================================================
-- 2. DIMENSION: CHECKV QUALITY BENCHMARKS
-- Source: MVP Module 03 representative vOTU quality summary
-- Columns mapped from: checkv_quality, completeness, completeness_method,
--                       virus_length, provirus, gene_count, viral_genes,
--                       host_genes, kmer_freq, miuvig_quality
-- ============================================================================
CREATE TABLE dim_quality_benchmarks (
    quality_id            INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    contig_id             VARCHAR(255) NOT NULL UNIQUE,
    checkv_quality        VARCHAR(50)  NOT NULL,
    miuvig_quality        VARCHAR(50),
    completeness          NUMERIC(6, 2) DEFAULT 0.00,
    completeness_method   VARCHAR(80),
    virus_length_bp       INT          NOT NULL,
    provirus              BOOLEAN      DEFAULT FALSE,
    provirus_coordinates  VARCHAR(50),
    gene_count            INT          DEFAULT 0,
    viral_genes           INT          DEFAULT 0,
    host_genes            INT          DEFAULT 0,
    kmer_freq             NUMERIC(5, 2),

    CONSTRAINT chk_completeness_range CHECK (completeness BETWEEN 0.00 AND 100.00),
    CONSTRAINT chk_quality_tier CHECK (
        checkv_quality IN ('Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined')
    )
);

-- ============================================================================
-- 3. DIMENSION: VIRAL TAXONOMY & GENOMAD CLASSIFICATION
-- Source: MVP Module 03 representative vOTU quality summary + Module 06 annotations
-- Columns mapped from: virus_score, taxonomy, Genome type, Host type,
--                       n_hallmarks, marker_enrichment, genetic_code
-- ============================================================================
CREATE TABLE dim_viral_taxonomy (
    votu_id              INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    contig_id            VARCHAR(255) NOT NULL UNIQUE,
    cluster_members      TEXT,                                   -- comma-separated member IDs
    representative_sample VARCHAR(100),                          -- sample that contributed the representative
    virus_score          NUMERIC(5, 4) NOT NULL,
    genetic_code         INT          DEFAULT 11,
    n_hallmarks          INT          DEFAULT 0,
    marker_enrichment    NUMERIC(7, 4),
    predicted_taxonomy   VARCHAR(512) DEFAULT 'Unclassified',    -- full taxonomy string
    genome_type          VARCHAR(20)  DEFAULT 'dsDNA',           -- dsDNA, ssDNA, RNA, etc.
    host_type            VARCHAR(20)  DEFAULT 'Prokaryote',      -- Prokaryote or Eukaryote
    raw_genomad_annotations JSONB,

    CONSTRAINT chk_virus_score CHECK (virus_score BETWEEN 0.0000 AND 1.0000)
);

-- ============================================================================
-- 4. DIMENSION: FUNCTIONAL GENE ANNOTATIONS
-- Source: MVP Module 06 geNomad/PHROGs/PFAM annotation tables
-- ============================================================================
CREATE TABLE dim_gene_annotations (
    gene_id              INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    viral_gene_id        VARCHAR(255) NOT NULL,
    contig_id            VARCHAR(255) NOT NULL,
    gene_number          INT,
    gene_start           INT,
    gene_end             INT,
    gene_length          INT,
    strand               INT,                                    -- 1 or -1
    genomad_annotation   TEXT         DEFAULT 'Unknown',
    genomad_accession    VARCHAR(100),
    genomad_score        NUMERIC(10, 3),
    genomad_evalue       DOUBLE PRECISION,
    gc_content           NUMERIC(5, 3),
    rbs_motif            VARCHAR(30),
    genomad_marker       VARCHAR(100),
    virus_hallmark       VARCHAR(50),
    amr_annotation       VARCHAR(255),

    CONSTRAINT fk_gene_contig FOREIGN KEY (contig_id)
        REFERENCES dim_viral_taxonomy(contig_id) ON DELETE CASCADE
);

-- ============================================================================
-- 5. FACT TABLE: vOTU READ-RECRUITMENT QUANTITATIVE METRICS
-- Source: MVP Module 04 CoverM output per sample
-- Columns mapped from: Mean, Trimmed Mean, Covered Bases, Covered Fraction,
--                       Variance, Read Count, Reads per base, RPKM, TPM
-- ============================================================================
CREATE TABLE fact_votu_metrics (
    metric_id        BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    sample_id        INT    NOT NULL,
    votu_id          INT    NOT NULL,
    quality_id       INT    NOT NULL,

    -- Coverage metrics (from CoverM)
    mean_depth           NUMERIC(12, 4) DEFAULT 0.0000,
    trimmed_mean_depth   NUMERIC(12, 4) DEFAULT 0.0000,
    covered_bases        INT            DEFAULT 0,
    covered_fraction     NUMERIC(7, 6)  DEFAULT 0.000000,       -- 0.0 to 1.0 scale
    coverage_variance    NUMERIC(12, 4) DEFAULT 0.0000,
    contig_length        INT            DEFAULT 0,

    -- Abundance metrics
    raw_read_count       INT            DEFAULT 0,
    reads_per_base       NUMERIC(12, 8) DEFAULT 0.00000000,
    rpkm                 NUMERIC(12, 4) DEFAULT 0.0000,
    tpm                  NUMERIC(14, 4) DEFAULT 0.0000,

    CONSTRAINT fk_fact_sample  FOREIGN KEY (sample_id)  REFERENCES dim_samples(sample_id)            ON DELETE RESTRICT,
    CONSTRAINT fk_fact_votu    FOREIGN KEY (votu_id)     REFERENCES dim_viral_taxonomy(votu_id)       ON DELETE RESTRICT,
    CONSTRAINT fk_fact_quality FOREIGN KEY (quality_id)  REFERENCES dim_quality_benchmarks(quality_id) ON DELETE RESTRICT,

    CONSTRAINT chk_covered_fraction CHECK (covered_fraction BETWEEN 0.000000 AND 1.000001),
    CONSTRAINT chk_raw_reads        CHECK (raw_read_count >= 0),
    CONSTRAINT uq_sample_votu       UNIQUE (sample_id, votu_id)
);

-- ============================================================================
-- 6. PERFORMANCE OPTIMIZATION LAYER
-- ============================================================================

-- Fast lookup: metrics by sample + vOTU
CREATE INDEX idx_fact_sample_votu     ON fact_votu_metrics (sample_id, votu_id);

-- Partial index: only vOTUs with significant coverage (>= 75% of genome covered)
CREATE INDEX idx_fact_high_coverage   ON fact_votu_metrics (covered_fraction) WHERE covered_fraction >= 0.75;

-- RPKM-based abundance ranking
CREATE INDEX idx_fact_rpkm            ON fact_votu_metrics (rpkm DESC);

-- JSONB index for flexible geNomad annotation queries
CREATE INDEX idx_gin_genomad          ON dim_viral_taxonomy USING gin (raw_genomad_annotations);

-- Taxonomy text search
CREATE INDEX idx_taxonomy_text        ON dim_viral_taxonomy (predicted_taxonomy);

-- Gene annotation lookup by contig
CREATE INDEX idx_gene_contig          ON dim_gene_annotations (contig_id);

-- Quality tier filtering
CREATE INDEX idx_quality_tier         ON dim_quality_benchmarks (checkv_quality);
