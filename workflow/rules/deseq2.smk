rule deseq2:
    input:
        raw_counts = config['raw_counts'],
        coldata = config['coldata'],
        contrasts = config['contrasts']
    output:
        f"{RESULTS_DIR}/Deseq2_normalized_counts.txt",
        f"{RESULTS_DIR}/HistDesq2normFilter_mqc.png",
        expand(f"{RESULTS_DIR}/MAplot_{{contrast}}_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/DEG_{{contrast}}.txt", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/DEG_{{contrast}}_filtered.txt", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/heatmap_{{contrast}}_mqc.png", contrast=CONTRASTS)
    conda:
        "envs/R.yaml"
    log:
        f"{RESULTS_DIR}/deseq2.log"
    shell:
        """
        Rscript /workflow/rules/scripts/deseq2.R {input.raw_counts} {input.coldata} {input.contrasts} {RESULTS_DIR} > {log} 2>&1
        """

