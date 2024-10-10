rule venn_diagram:
    input:
        deseq2=f"{RESULTS_DIR}/DEG_{{contrast}}_filtered.txt",
        edger=f"{RESULTS_DIR}/Filtered_DEG_EdgeR_{{contrast}}.txt",
        limma_voom=f"{RESULTS_DIR}/limma_voom_filtered_results_{{contrast}}.txt"
    output:
        venn_diagram=f"{RESULTS_DIR}/venn_diagram_{{contrast}}_mqc.png"
    log:
        f"{RESULTS_DIR}/venn_diagram_{{contrast}}.log"
    shell:
        """
        Rscript /workflow/rules/scripts/venn_diagram.R \
        {input.deseq2} {input.edger} {input.limma_voom} {output.venn_diagram} {log} > {log} 2>&1
        """

