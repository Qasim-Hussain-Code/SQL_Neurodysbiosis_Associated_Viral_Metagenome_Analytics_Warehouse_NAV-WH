WITH RankedVirome AS (
    SELECT 
        s.cohort_type,
        t.votu_id,
        t.contig_id,
        f.relative_abundance,
        DENSE_RANK() OVER(PARTITION BY s.cohort_type ORDER BY f.relative_abundance DESC) as abundance_rank
    FROM fact_votu_metrics f
    JOIN dim_samples s ON f.sample_id = s.sample_id
    JOIN dim_viral_taxonomy t ON f.votu_id = t.votu_id
    WHERE f.horizontal_coverage >= 75.00
)
SELECT cohort_type, votu_id, contig_id, relative_abundance
FROM RankedVirome
WHERE abundance_rank <= 5;
