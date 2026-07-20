-- ============================================================================
-- QUERY 03: High-Confidence vOTU Materialized View + Differential Abundance
-- Purpose:  Pre-compute a curated set of publication-quality vOTUs and
--           compare their abundance between cohorts.
-- ============================================================================

-- 3a. Materialized view: vOTUs with Medium-quality or better (MIUViG standard)
CREATE MATERIALIZED VIEW mv_high_confidence_votus AS
SELECT
    t.votu_id,
    t.contig_id,
    t.cluster_members,
    t.virus_score,
    t.predicted_taxonomy,
    t.genome_type,
    t.n_hallmarks,
    q.checkv_quality,
    q.completeness,
    q.virus_length_bp,
    q.provirus,
    q.gene_count,
    q.viral_genes,
    q.host_genes
FROM dim_viral_taxonomy t
JOIN dim_quality_benchmarks q ON t.contig_id = q.contig_id
WHERE q.checkv_quality IN ('Complete', 'High-quality', 'Medium-quality')
WITH DATA;

CREATE UNIQUE INDEX idx_mv_hc_votu_id ON mv_high_confidence_votus (votu_id);


-- 3b. Differential abundance comparison: Neurodysbiosis vs Neurotypical
WITH cohort_rpkm AS (
    SELECT
        t.contig_id,
        t.predicted_taxonomy,
        q.checkv_quality,
        q.virus_length_bp,
        q.provirus,
        MAX(f.rpkm) FILTER (WHERE s.cohort_type = 'Neurodysbiosis_Cohort')  AS rpkm_neurodysbiosis,
        MAX(f.rpkm) FILTER (WHERE s.cohort_type = 'Neurotypical_Control')   AS rpkm_neurotypical,
        MAX(f.covered_fraction) FILTER (WHERE s.cohort_type = 'Neurodysbiosis_Cohort') AS cov_neurodysbiosis,
        MAX(f.covered_fraction) FILTER (WHERE s.cohort_type = 'Neurotypical_Control')  AS cov_neurotypical
    FROM fact_votu_metrics f
    JOIN dim_samples s        ON f.sample_id = s.sample_id
    JOIN dim_viral_taxonomy t ON f.votu_id   = t.votu_id
    JOIN dim_quality_benchmarks q ON f.quality_id = q.quality_id
    GROUP BY t.contig_id, t.predicted_taxonomy, q.checkv_quality, q.virus_length_bp, q.provirus
)
SELECT
    contig_id,
    checkv_quality,
    virus_length_bp,
    provirus,
    COALESCE(rpkm_neurodysbiosis, 0) AS rpkm_nd,
    COALESCE(rpkm_neurotypical, 0)   AS rpkm_nt,
    CASE
        WHEN COALESCE(rpkm_neurotypical, 0) = 0 AND COALESCE(rpkm_neurodysbiosis, 0) > 0
            THEN 'ND_exclusive'
        WHEN COALESCE(rpkm_neurodysbiosis, 0) = 0 AND COALESCE(rpkm_neurotypical, 0) > 0
            THEN 'NT_exclusive'
        WHEN rpkm_neurodysbiosis / NULLIF(rpkm_neurotypical, 0) >= 2
            THEN 'ND_enriched'
        WHEN rpkm_neurotypical / NULLIF(rpkm_neurodysbiosis, 0) >= 2
            THEN 'NT_enriched'
        ELSE 'Shared'
    END AS enrichment_pattern,
    ROUND(COALESCE(cov_neurodysbiosis, 0)::NUMERIC * 100, 1) AS coverage_pct_nd,
    ROUND(COALESCE(cov_neurotypical, 0)::NUMERIC * 100, 1)   AS coverage_pct_nt,
    predicted_taxonomy
FROM cohort_rpkm
WHERE COALESCE(rpkm_neurodysbiosis, 0) + COALESCE(rpkm_neurotypical, 0) > 0
ORDER BY
    CASE
        WHEN COALESCE(rpkm_neurotypical, 0) = 0 THEN 1
        WHEN COALESCE(rpkm_neurodysbiosis, 0) = 0 THEN 2
        ELSE 3
    END,
    GREATEST(COALESCE(rpkm_neurodysbiosis, 0), COALESCE(rpkm_neurotypical, 0)) DESC;
