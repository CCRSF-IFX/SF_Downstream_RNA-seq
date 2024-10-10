rule edger:
    input:
        raw_counts = config['raw_counts'],
        coldata = config['coldata'],
        contrasts = config['contrasts']
    output:
        touch(expand(f"{RESULTS_DIR}/DEG_EdgeR_{{contrast}}.txt", contrast=CONTRASTS)),
        touch(expand(f"{RESULTS_DIR}/Filtered_DEG_EdgeR_{{contrast}}.txt", contrast=CONTRASTS)),
        touch(expand(f"{RESULTS_DIR}/Smearplot_EdgeR_{{contrast}}_mqc.png", contrast=CONTRASTS)),
        touch(expand(f"{RESULTS_DIR}/edgeR_heatmaps_samplebysample_mqc_{{contrast}}.png", contrast=CONTRASTS)),
        touch(expand(f"{RESULTS_DIR}/edgeR_prcomp_mqc_{{contrast}}.png", contrast=CONTRASTS)),
        touch(expand(f"{RESULTS_DIR}/edgeR_pca_mqc_{{contrast}}.png", contrast=CONTRASTS)),
        touch(f"{RESULTS_DIR}/libdistrib_mqc.png"),
        touch(f"{RESULTS_DIR}/MDS_bcv_mqc.png"),
        touch(f"{RESULTS_DIR}/MDS_logFC_mqc.png"),
        touch(f"{RESULTS_DIR}/BCVplot_mqc.png"),
        touch(f"{RESULTS_DIR}/HistEdgeRnormFilter_mqc.png"),
        touch(f"{RESULTS_DIR}/edgeR_normalized_counts.txt"),
        touch(f"{RESULTS_DIR}/edgeR_normalized_counts_log.txt")
    conda:
        "envs/R.yaml"
    log:
        f"{RESULTS_DIR}/edger.log"
    shell:
        """
        Rscript /workflow/rules/scripts/edger.R {input.raw_counts} {input.coldata} {input.contrasts} {RESULTS_DIR} > {log} 2>&1
        """

