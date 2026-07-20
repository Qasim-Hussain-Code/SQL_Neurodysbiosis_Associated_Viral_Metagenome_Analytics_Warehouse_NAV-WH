-- ============================================================================
-- QUERY 01: Comparative Cohort vOTU Abundance Ranking
-- Purpose:  Identify the top-N most abundant vOTUs per cohort using RPKM,
--           filtered by a minimum genome coverage threshold (75%).
-- ============================================================================

WITH ranked_virome AS (
    SELECT
        s.cohort_type,
        t.votu_id,
        t.contig_id,
        t.predicted_taxonomy,
        q.checkv_quality,
        q.completeness,
        q.virus_length_bp,
        f.rpkm,
        f.raw_read_count,
        f.covered_fraction,
        DENSE_RANK() OVER (
            PARTITION BY s.cohort_type
            ORDER BY f.rpkm DESC
        ) AS abundance_rank
    FROM fact_votu_metrics f
    JOIN dim_samples s        ON f.sample_id = s.sample_id
    JOIN dim_viral_taxonomy t ON f.votu_id   = t.votu_id
    JOIN dim_quality_benchmarks q ON f.quality_id = q.quality_id
    WHERE f.covered_fraction >= 0.75
      AND f.rpkm > 0
)
SELECT
    cohort_type,
    abundance_rank,
    contig_id,
    checkv_quality,
    completeness,
    virus_length_bp,
    rpkm,
    raw_read_count,
    ROUND(covered_fraction::NUMERIC * 100, 1) AS coverage_pct,
    predicted_taxonomy
FROM ranked_virome
WHERE abundance_rank <= 10
ORDER BY cohort_type, abundance_rank;
