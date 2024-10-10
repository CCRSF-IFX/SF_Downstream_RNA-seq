rule pathway_analysis:
    input:
        deseq2_file = f"{RESULTS_DIR}/DEG_{{contrast}}_filtered.txt"
    output:
        dotplot = f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_dotplot_mqc.png",
        emapplot = f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_emapplot_mqc.png",
        cnetplot = f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_cnetplot_mqc.png",
        ridgeplot = f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_ridgeplot_mqc.png",
        gseaplot = f"{RESULTS_DIR}/pathway_analysis_deseq2_{{contrast}}_gseaplot_mqc.png",
        kegg_dotplot = f"{RESULTS_DIR}/pathway_analysis_kegg_{{contrast}}_dotplot_mqc.png",
        kegg_emapplot = f"{RESULTS_DIR}/pathway_analysis_kegg_{{contrast}}_emapplot_mqc.png",
        kegg_cnetplot = f"{RESULTS_DIR}/pathway_analysis_kegg_{{contrast}}_cnetplot_mqc.png",
        kegg_ridgeplot = f"{RESULTS_DIR}/pathway_analysis_kegg_{{contrast}}_ridgeplot_mqc.png",
        kegg_gseaplot = f"{RESULTS_DIR}/pathway_analysis_kegg_{{contrast}}_gseaplot_mqc.png",
        #pathview_output = f"{RESULTS_DIR}/{config['kegg_organism']}04130.{{contrast}}_pathview.png"
    params:
        genome = config["genome"],
        organism_db = config["organism_db"],
        kegg_organism = config["kegg_organism"],
        pathway_id = config["pathway_id"]
    log:
        f"{RESULTS_DIR}/pathway_analysis_{{contrast}}.log"
    shell:
        """
        Rscript /workflow/rules/scripts/pathway.R \
            {input.deseq2_file} \
            {wildcards.contrast} \
            {RESULTS_DIR} \
            {params.genome} \
            {params.organism_db} \
            {params.kegg_organism} \
            {params.pathway_id} > {log} 2>&1
        """

