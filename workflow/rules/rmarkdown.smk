rule rmarkdown_report:
    input:
        expand("{results}/DEG_{contrast}_filtered.txt", results=config["results"], contrast=CONTRASTS),
        expand("{results}/Filtered_DEG_EdgeR_{contrast}.txt", results=config["results"], contrast=CONTRASTS),
        expand("{results}/limma_voom_filtered_results_{contrast}.txt", results=config["results"], contrast=CONTRASTS),
        expand("{results}/venn_diagram_{contrast}_mqc.png", results=config["results"], contrast=CONTRASTS),
    output:
        "{results}/RNA_seq_analysis_report.html"
    conda:
        "envs/R.yaml"
    params:
        config_path=config["config_path"],
        results_dir=config["results"],
        script_path="/workflow/rules/scripts/rmarkdown_report.rmd"  # Correct relative path to the Rmd file
    shell:
        """
        Rscript -e "rmarkdown::render('{params.script_path}', params=list(config_path='{params.config_path}'), intermediates_dir='{params.results_dir}', output_file='RNA_seq_analysis_report.html', output_dir='{params.results_dir}', quiet=FALSE)"
        """

