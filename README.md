# Sequencing Facility RNA-Seq Downstream Analysis pipeline
Pipeline will intake raw counts file and run differential expression analysis using three different tools (Deseq2, EdgeR, and Limma_Voom) and performs pathway analysis using Deseq2 output.

[link of jpeg lucid workflow]

![SF_Fastq-QC](https://github.com/CCRSF-IFX/SF_Fastq-QC/blob/main/resources/Fastq-QC_workflow.jpeg CHANGE LINK)

## Analysis Tools:
[Deseq2](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html): From Deseq2, users can expect a comprehensive analysis of differential gene expression using DESeq2. The script takes as input raw gene counts, sample metadata, and contrasts between experimental conditions. It normalizes the gene counts, performs differential expression analysis, and outputs normalized counts, lists of differentially expressed genes (DEGs), and corresponding MA plots. The DEGs are filtered based on a log fold change threshold of 2 for upregulation and 0.5 for downregulation, with a significance threshold of p-value < 0.05. Additionally, heatmaps for the top 30 differentially expressed genes are generated to visually represent expression changes across the relevant samples. The code ensures that all outputs, including tables and plots, are saved to a specified directory for easy retrieval and downstream analysis.

[EdgeR]( https://bioconductor.org/packages/release/bioc/html/edgeR.html): The edgeR code performs differential gene expression analysis by normalizing raw RNA-seq counts using TMM, accounting for library size differences, and estimating gene-specific dispersions. Users can expect the code to detect significantly differentially expressed genes between experimental groups using an exact test. It filters results based on an FDR threshold of 0.05 and log fold change thresholds of ≥ 2 or ≤ 0.5 (corresponding to a log2 fold change of ±1). The code generates visualizations such as MDS plots, MA plots, and heatmaps to help users explore gene expression patterns across samples.

[Limma_Voom](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html): The limma-voom code applies the voom transformation to RNA-seq data, stabilizing variance across gene expression levels and making the data suitable for linear modeling. Users can expect the code to fit linear models to the transformed data and perform empirical Bayes moderation to improve statistical accuracy. Contrasts between experimental groups are analyzed for differential expression, with results filtered based on an adjusted p-value threshold of 0.05 and a log fold change threshold of ≥ 2 or ≤ -2 (corresponding to a log2 fold change of ±1). The code generates MA plots and heatmaps to visualize the top differentially expressed genes, providing a detailed summary of expression changes.

[Venn_Diagram](https://r-graph-gallery.com/14-venn-diagramm): This Venn diagram script takes filtered gene lists as input from multiple differential expression analyses, specifically from DESeq2, edgeR, and limma-voom. These gene lists are used to generate a Venn diagram, which visually represents the overlap of differentially expressed genes identified by each method. The script processes the input files containing filtered genes, ensures that the required GeneID column is present, and creates the diagram. The resulting Venn diagram is saved as a PNG file, providing a clear visualization of shared and unique genes between the analyses.  

[Pathway_Analysis]( https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/): Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether a pre-defined set of genes (ex: those beloging to a specific GO term or KEGG pathway) shows statistically significant, concordant differences between two biological states.

## Usage

### Step 1: Obtain a copy of this workflow
Clone the Repository: Clone the new repository to your local machine, choosing the directory where you want to perform data analysis. Instructions for cloning can be found [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).

### Step 2: Configure workflow

1) Have input RawCountFile_rsemgenes.txt available in test folder.  
-Example of input RawCountFile_rsemgenes.txt
gene_id	Sample1	Sample2	Sample3	Sample4	Sample5	Sample6
ENSG00000000003.14_TSPAN6	1500	1600	1450	1300	1250	1400
ENSG00000000005.6_TNMD	10	12	9	15	14	18
ENSG00000000419.12_DPM1	3000	3100	2900	3200	3300	3400
ENSG00000000457.14_SCYL3	500	520	480	600	590	610
.
.
.
.
.
ENSG00000005955.13_ACTN2	500	530	460	550	540	570
ENSG00000006062.15_RPS6	600	650	580	700	680	720

2) Modify the contrasts.txt file in config folder. 
-Add sample condition contrast pairs. Example: 
parent control

3) Modify the sampleinfo.txt file in config folder.
-Add sample details keeping the same format such as:
samplename	condition	batch 
sample1	parent		b1
sample2	parent		b2
sample3	parent		b2
sample4	control		b1
sample5	control		b2
sample6	control		b2

4) In the /config/config.yaml:
-Add full path to the raw_counts data i.e. /test/RawCountFile_rsemgenes.txt
-Add full path to results directory i.e. /results
-Add full path to contrasts file i.e /config/contrasts.txt
-Select the reference genome of choice. For example: “mouse”
-Select the organism_db accordingly. For example:
-Select the kegg_organism. For example: "mmu"

5) In the /workflow/Snakefile:
-Give full path to configfile in Snakefile
-Give full path to rule files in Snakefile 

6) In the /workflow/rules/common.smk 
-Give full path to configfile

6) Ensure full paths to the Rscripts in snakemake (.smk) files under 'shell' or "params' for the following:
-deseq2.smk
-pathway.smk
-edger.smk
-venn.smk
-rmarkdown.smk
-limma_voom.smk

7) Ensure full path to R libraries as .libPaths in /workflow/rules/scripts/pathway.R

### Step 3: Load the snakemake version 8 or above 

`module load snakemake/8.4.8`

### Step 4: Create a conda environment

`$NAME=my_environment_name`

`conda create -n $NAME`

### Step 5: Execute the Workflow

Activate the Conda Environment:

`conda activate $NAME`

Install `mamba`

`conda install -c conda-forge mamba`

Perform a dry-run to validate your setup:

`snakemake --use-conda -np`

Local Execution: Execute the workflow on your local machine using $N cores:

`snakemake --use-conda --cores $N --latency-wait 55 --force`

Here, $N represents the number of cores you wish to allocate for the workflow.

### Step 6: Investigate results

After successful execution, output analysis files along with rmarkdown html report will be generated in results folder. 

Move the results folder to a different path and the workflow will be ready to run another dataset. 

`conda deactivate`


## For any questions, contact
CCRSF_IFX@nih.gov

