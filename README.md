# Fiber-seq in plants

<br>

## Fiber-seq analysis workflow

This workflow is used to run the **Fiber-seq** analysis workflow in different plants (such as ***Zea mays***, ***Arabidopsis thaliana*** etc.) to study chromatin structure, accessibility, and DNA methylation modifications.  

When running on a **SLURM** cluster, choose the corresponding script file according to the species.

🔹 Example: **Zea mays** B73 

Run the following command in the terminal:

        sbatch Fiber-seq_Zeamays_B73.sh data sample.bam

This generates a series of files, mainly including:

1. Open chromatin signals: `*fire.bw`, `*percent.bw`

2. Nucleosome signals: `*nuc.bw`, `*linker.bw`

3. 5mC methylation: `*CPG*.bw`

4. raw m6A : `*cram`,`*bam`

5. High-confidence peaks: `*fire-v0.1.1-wide-peaks.bed`
   
6. Low-confidence peaks: `*fire-v0.1.1-01-fire-wide-peaks.bed`

<br>
<br>

## Fiber-seq differential peaks analysis workflow

### 1. Run `FIRE_Get_differential_peaks.sh`:
   

        sh FIRE_Get_differential_peaks.sh 'sample1 sample2' 
         

  Requires the Peak_Diff_V1.R file and should be run in the base environment.

  **Each p-value will generate two files at 0.05 intervals, containing BED files of the upregulated and downregulated peaks in the two samples.**

<br>
<br>

### 2. Run `FIRE_Differentail_Peak_visualization.R` (use `--help` to view the help information)  

It is generally recommended to use this script when there are more than three samples.

This will generate clustering heatmaps and boxplots based on z-scores.

<br>
<br>

## ATAC-seq analysis workflow

When running on a **SLURM** cluster, choose the corresponding script file according to the species.

🔹 Example: **Arabidopsis thaliana**

When you want to process the files `sample_R1.fq.gz` and `sample_R2.fq.gz`

Run the following command in the terminal:
    
        sbatch ATAC-seq_Arabidopsis_Col.sh date sample

A folder named `date_sample` will be generated, containing the corresponding open chromatin region `.bw` files for the sample.

<br>
<br>

## ATAC-seq callpeak workflow

To perform peak calling for the above ATAC-seq data, you can run the following command:

        sbatch ATAC-seq_CallPeak_Arabidopsis_Col.sh date_sample sample
        
This will generate a folder named `date_sample_analysis`, which contains files with the peak counts.

<br>
<br>

## RNA-seq analysis workflow

When running on a **SLURM** cluster, choose the corresponding script file according to the species.

When you want to process the files `sample_R1.fq.gz` and `sample_R2.fq.gz`

🔹 Example: **Zea mays** B73 
    
        sbatch RNA-seq_Zeamays_B73.sh date sample

This generates the .bw files of gene expression for the corresponding sample.

<br>
<br>

## Genome assembly workflow

Activate the Assembly virtual environment

This workflow does not involve Hi-C data, so a reference genome (Control.fa) is required to assist with chromosome assembly.

Run the script:

        sh Genome_assembly_pipeline.sh data Control.fa sample.bam 

This generates the genome assembled from the input BAM file, as well as visualization plots comparing the newly assembled genome with the reference genome.

<br>
<br>

## Overlap

### 1. `Get_IDOverlap.R` can be used to obtain the overlap between two gene ID files

        Rscript Get_IDOverlap.R A.bed B.bed

        
### 2. `Venn_Pvalue.R` can be used to obtain the overlap between two BED files. (use `--help` to view the help information)

Options:

        -a A_BED, --A_bed=A_BED
                Path to first BED file (A)

        -b B_BED, --B_bed=B_BED
                Path to second BED file (B)

        -g GENOME, --genome=GENOME
                Path to genome length file

        -o OUTPUT, --output=OUTPUT
                Output PDF filename [default: Venn_with_Pvalue.pdf]

        -h, --help
                Show this help message and exit

<br>
<br>

## All software versions used in this workflow have been uploaded as Conda YAML files in the `envs` folder.

<br>
<br>

## Acknowledgements

We would like to thank all members of our research group for their helpful discussions and technical support during this project. We also acknowledge the developers of *Fibertools* for providing powerful tools that enabled the integration of long-read epigenetic and genetic analyses in our study. The following publication was referenced in this work:

Jha, A., Bohaczuk, S. C., Mao, Y., Ranchalis, J., Mallory, B. J., Min, A. T., Hamm, M. O., Swanson, E., Dubocanin, D., Finkbeiner, C., Li, T., Whittington, D., Noble, W. S., Stergachis, A. B., & Vollger, M. R. (2024). DNA-m6A calling and integrated long-read epigenetic and genetic analysis with *Fibertools*. *Genome Research*. https://doi.org/10.1101/gr.279095.124






