-- ============================================================================
-- QUERY 02: Functional Gene Annotation Summary with Hallmark Detection
-- Purpose:  Extract and summarize geNomad functional annotations, highlighting
--           virus hallmark genes and their distribution across vOTUs.
-- ============================================================================

-- 2a. Hallmark gene inventory per vOTU
SELECT
    g.contig_id,
    t.predicted_taxonomy,
    q.checkv_quality,
    q.completeness,
    COUNT(*)                                          AS total_genes,
    COUNT(*) FILTER (WHERE g.virus_hallmark = 'Virus_hallmark') AS hallmark_genes,
    STRING_AGG(
        DISTINCT g.genomad_annotation, '; '
    ) FILTER (WHERE g.virus_hallmark = 'Virus_hallmark') AS hallmark_functions
FROM dim_gene_annotations g
JOIN dim_viral_taxonomy t      ON g.contig_id = t.contig_id
JOIN dim_quality_benchmarks q  ON g.contig_id = q.contig_id
GROUP BY g.contig_id, t.predicted_taxonomy, q.checkv_quality, q.completeness
ORDER BY hallmark_genes DESC, total_genes DESC;


-- 2b. Top functional categories across all vOTUs (excluding unknowns)
SELECT
    g.genomad_annotation,
    COUNT(*)               AS gene_count,
    COUNT(DISTINCT g.contig_id) AS votu_count,
    BOOL_OR(g.virus_hallmark = 'Virus_hallmark') AS is_hallmark
FROM dim_gene_annotations g
WHERE g.genomad_annotation NOT IN ('Unknown', 'NA', '')
GROUP BY g.genomad_annotation
ORDER BY gene_count DESC
LIMIT 25;


-- 2c. JSONB annotation extraction (when raw_genomad_annotations is populated)
SELECT
    contig_id,
    virus_score,
    genome_type,
    host_type,
    predicted_taxonomy,
    raw_genomad_annotations ->> 'topology'                           AS topology,
    (raw_genomad_annotations -> 'classification' ->> 'fdr')::NUMERIC AS fdr
FROM dim_viral_taxonomy
WHERE raw_genomad_annotations IS NOT NULL;
