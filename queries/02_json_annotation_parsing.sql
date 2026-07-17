SELECT 
    contig_id,
    genomad_score,
    raw_genomad_annotations->>'topology' AS extracted_topology,
    CAST(raw_genomad_annotations->'classification'->>'virus_score' AS DECIMAL(4,3)) AS internal_virus_score
FROM dim_viral_taxonomy
WHERE raw_genomad_annotations IS NOT NULL;
