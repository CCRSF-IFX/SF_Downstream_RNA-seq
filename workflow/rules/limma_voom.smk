rule limma_voom:
    input:
        raw_counts = config['raw_counts'],
        coldata = config['coldata'],
        contrasts = config['contrasts']
    output:
        expand(f"{RESULTS_DIR}/limma_voom_results_{{contrast}}.txt", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/limma_voom_filtered_results_{{contrast}}.txt", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/MAplot_limma_voom_{{contrast}}_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/heatmap_limma_voom_{{contrast}}_mqc.png", contrast=CONTRASTS)
    conda:
        "envs/R.yaml"
    log:
        f"{RESULTS_DIR}/limma_voom.log"
    shell:
        """
        Rscript /workflow/rules/scripts/limma_voom.R {input.raw_counts} {input.coldata} {input.contrasts} {RESULTS_DIR} > {log} 2>&1
        """

