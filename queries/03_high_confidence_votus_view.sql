CREATE MATERIALIZED VIEW mv_high_confidence_votus AS
SELECT t.votu_id, t.contig_id, q.completeness_tier, q.contamination_percentage
FROM dim_viral_taxonomy t
JOIN fact_votu_metrics f ON t.votu_id = f.votu_id
JOIN dim_quality_benchmarks q ON f.quality_id = q.quality_id
WHERE q.completeness_tier IN ('Complete', 'High-quality') AND q.contamination_percentage < 5.00
WITH DATA;

CREATE UNIQUE INDEX idx_mv_votu_id ON mv_high_confidence_votus (votu_id);
