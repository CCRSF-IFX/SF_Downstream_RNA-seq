rule multiqc:
    input:
        # Including only the necessary PNG files for MultiQC
        expand(f"{RESULTS_DIR}/MAplot_{{contrast}}_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/heatmap_{{contrast}}_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/Smearplot_EdgeR_{{contrast}}_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/edgeR_heatmaps_samplebysample_mqc_{{contrast}}.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/edgeR_prcomp_mqc_{{contrast}}.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/edgeR_pca_mqc_{{contrast}}.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/MAplot_limma_voom_{{contrast}}_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/heatmap_limma_voom_{{contrast}}_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/venn_diagram_{{contrast}}_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_dotplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_emapplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_cnetplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_ridgeplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_gseaplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_edger_{{contrast}}_dotplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_edger_{{contrast}}_emapplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_edger_{{contrast}}_cnetplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_edger_{{contrast}}_ridgeplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_edger_{{contrast}}_gseaplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_limma_voom_{{contrast}}_dotplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_limma_voom_{{contrast}}_emapplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_limma_voom_{{contrast}}_cnetplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_limma_voom_{{contrast}}_ridgeplot_mqc.png", contrast=CONTRASTS),
        expand(f"{RESULTS_DIR}/pathway_analysis_limma_voom_{{contrast}}_gseaplot_mqc.png", contrast=CONTRASTS),
        f"{RESULTS_DIR}/HistDesq2normFilter_mqc.png",
        f"{RESULTS_DIR}/libdistrib_mqc.png",
        f"{RESULTS_DIR}/MDS_bcv_mqc.png",
        f"{RESULTS_DIR}/MDS_logFC_mqc.png",
        f"{RESULTS_DIR}/BCVplot_mqc.png",
        f"{RESULTS_DIR}/HistEdgeRnormFilter_mqc.png"
    output:
        f"{RESULTS_DIR}/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc {RESULTS_DIR} -o {RESULTS_DIR} > {RESULTS_DIR}/multiqc.log 2>&1
        """

