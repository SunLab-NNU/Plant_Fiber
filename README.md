# Plant Fiber Project

## All Conda virtual environment YAML files have been uploaded envs folder.

## Fiber-seq Analysis Workflow

This project is used to run the **Fiber-seq** analysis workflow in different plants (such as maize B73, Arabidopsis, peanut, etc.) to study chromatin structure, accessibility, and DNA methylation modifications.  

When running on a **SLURM** cluster, choose the corresponding script file according to the species.

🔹 Example: Maize B73 Command

Run the following command in the terminal:


        sbatch B73_V5.sh 251005 yumifibercj.bc2065.bc2065.HiFi.bam

This generates a series of files, mainly including:

1. Open chromatin regions: `*fire.bw`, `*percent.bw`

2. Nucleosome regions: `*nuc.bw`, `*linker.bw`

3. 5mC methylation: `*CPG*.bw`

4. raw m6A : `*cram`,`*bam`

5. High-confidence peaks: `*fire-v0.1.1-wide-peaks.bed`
   
6. Low-confidence peaks: `*fire-v0.1.1-01-fire-wide-peaks.bed`

## Fiber-seq Differential Analysis Workflow

Find Arabidopsis Fiber FIRE Peaks

1. Run `FIRE_Peak_UpStream.sh`, e.g.:
   

        sh FIRE_Peak_UpStream.sh 'sample1 sample2' 
         

     (Requires the Peak_Diff_V1.R file and should be run in the base environment.)

    #### Each p-value will generate two files, with intervals of 0.05 :

   BED files of the upregulated and downregulated peaks in the two samples

3. Run `Fiber_peak_diffbind_V1.R` (use `--help` to view the help information)  

    It is generally recommended to use this script when there are more than three samples.

    This will generate clustering heatmaps and boxplots based on z-scores.



## ATAC-seq Analysis Workflow

When running on a **SLURM** cluster, choose the corresponding script file according to the species.

🔹 Example: Arabidopsis Command

When you want to process the files `sample_R1.fq.gz` and `sample_R2.fq.gz`

Run the following command in the terminal:
    
        sbatch ATAC-seq_Arabidopsis_Col.sh date sample

A folder named `date_sample` will be generated, containing the corresponding open chromatin region `.bw` files for the sample.


## ATAC-seq CallPeak Workflow

To perform peak calling for the above ATAC-seq data, you can run the following command:

        sbatch ATAC-seq_CallPeak_Arabidopsis_Col.sh date_sample sample
        
This will generate a folder named `date_sample_analysis`, which contains files with the peak counts.


## RNA-seq Analysis Workflow

When running on a **SLURM** cluster, choose the corresponding script file according to the species.

When you want to process the files `sample_R1.fq.gz` and `sample_R2.fq.gz`

🔹 Example: Maize B73 Command
    
        sbatch RNA-seq_Zeamays_B73.sh date sample

This generates the .bw files of gene expression for the corresponding sample.



## Genome Assembly Workflow

Activate the Assembly virtual environment

Run the script:

   
        sh Genome_assembly_pipeline.sh data Control.fa sample.bam 

This generates the genome assembled from the input BAM file, as well as visualization plots comparing the newly assembled genome with the reference genome.

## overlap

`Get_IDOverlap.R` can be used to obtain the overlap between two gene ID files

`Venn_Pvalue.R` can be used to obtain the overlap between two BED files.




