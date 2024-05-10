# This repository contains the data used to generage figure 3H as well as supplemental information about the virus sequences derived from trubinates of immunocompromised mice with long-term SARS-CoV-2 beta variant infections and late-onset antiviral treatment. 

## <ins>**rava_output** contains the following data:</ins>
1. RAVA output directories, which contain variant analysis results with respect to the inoculum consensus genome. (https://github.com/greninger-lab/lava/tree/Rava_Slippage-Patch).

		- The RAVA output directories are separated by NGS library replicate (1 or 2) and timepoint (14, 21, or 28 dpi).
   		- There are also results for inoculum specimens at passage 3 and passage 4 (metagenomic and amplicon preparation).
		- These directories contain interactive html plots.
   
4. The inoculum consensus fasta sequence and gff annotation files.  

## <ins>**bin** contains the following scripts:</ins>
1. **figure3H_NGS_final.R**, which will produce the figure 3H.

		Input:
			1. The "visualization.csv" files in RAVA output directories. 
				- Each mouse turbinate specimen was sequenced twice using amplicon-tiling for SARS-CoV-2.
				- The "visualization" tables contain annotated information in tabular format for all
   				  strand-bias filtered variants with >=1% allele frequency.
				- The tables were derived from vcf files and are also the basis of the RAVA html plots.
   
		Outputs:
			1. figure3H_NGS_table.csv, the list of variants considered in the analysis.
   				- All variants appeared in both NGS replicates at allele frequency >=1%, depth >= 100.
   				- Furthermore, mean allele frequency >=5% and mean depth >= 100.
	   
			2. figure3H_NGS.pdf, figure 3H in pdf format.
	   
			3. stats_mean-var-per-animal_D21.csv, a list of p-values from One-way ANOVA analysis of
   			   mean variants per animal between treatment groups at 21 dpi.
	   
			4. stats_mean-var-per-animal_Vehicle.csv, a list of p-values from One-way ANOVA analysis of
   			   mean variants per animal between timepoints (14, 21, 28 dpi) for vehicle-treated animals.
	   
			5. stats_mean-var-per-animal_Nirmatrelvir.csv, a list of p-values from One-way ANOVA analysis of
   			   mean variants per animal between timepoints (21, 28 dpi) for Nirmatrelvir/Ritonavir-treated animals.
   	
3. **extract_fastqs_from_taylor.sh**, which extracts the soft-clipped reads from the bam files produced by the TAYLOR pipeline (https://github.com/greninger-lab/covid_swift_pipeline) and calls an R script (remove_soft_clips.R) that converts the bams to hard-clipped fastq files (i.e. a fastq file with reads that align to Wuhan-Hu-1 and do not contain amplicon primer sequences). 
4. **remove_soft_clips.R**, which reads the CIGAR strings from a sam file, convert soft-clips to hard-clips, and writes a new hard-clipped sam file. 

## <ins>**supplemental_tables** contains the following files:</ins>
1. **GISAID_Results_final.csv**, a table of variants from 4'-FlU and Nirmatrelvir/Ritonavir-treated animals at 21 dpi with mean allele frequency >=20% and mean depth >=100. It lists the variants in GISAID format with respect to the inoculum consensus genome and to Wuhan-Hu-1 (GISAID entry hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124), along with the number and percentage of genomes in the GISAID database that contain the SNVs as of May 2024.
2. **inoculum-consensus_vs_Wuhan-Hu-1.csv**, a table of nucelotide and amino acid changes in the inoculum consensus genome with respect to the Wuhan-Hu-1 genome (GenBank NC_045512.2).
3. **sample_metadata.csv**, a table of the BioProject, BioSample, and SRA accessions, along with treatment group, timepoint, animal number, NGS technical replicate, and "Passage" name used in RAVA pipeline.

## <ins>**Steps to reproduce this analysis**</ins>
<br>Note: I use nextflow tower with AWS S3, but you can modify the nextflow commands according to your machine's specifications.</br>
1. Download raw reads from NCBI SRA Archive
2. Run TAYLOR pipeline (https://github.com/greninger-lab/covid_swift_pipeline) on the raw reads:
		
  		nextflow run [path to covid_swift_pipeline/main.nf] \
			 --INPUT [path to raw reads from SRA archive] \
			 --OUTDIR [path to desired output directory] \
			 --MIN_LEN 50 \
			 -with-docker 'ubuntu:18.04' \
			 -with-report \
			 -c [path to the config file with your nextflow tower credentials] \
			 -profile cloud_bigger \
			 -with-tower \
			 --PRIMERS sarscov2_swift_v2
3. Run extract_fastqs_from_taylor.sh on the sorted bam files in the TAYLOR output directory:

   		./extract_fastqs_from_taylor.sh [path to bam files from TAYLOR output] [path to NC_045512.2.fasta, including filename]
	
5. Run RAVA on the extracted fastqs, using inoculum consensus fasta and gff files. For the RAVA metadata file, use the "Passage" names provided in **supplemental_tables/sample_metadata.csv** in this repository:

		nextflow run greninger-lab/lava \
			 -r Rava_Slippage-Patch \
			 --FASTA [full path to inoculum_consensus.fasta, including filename] \
			 --GFF [full path to inoculum_consensus.gff, including filename] \
			 --METADATA [full path to your metadata csv file, including filename] \
			 --OUTDIR [full path to your output directory] \
			 -with-docker 'ubuntu:18.04' \
			 -c [path to the config file with your nextflow tower credentials] \
			 -profile Cloud \
			 -with-tower
8. Run figure3H_NGS_final.R on the RAVA output. You can either edit the hard-coded rava parent directory in the R script, or provide it as the first command-line argument:

   		Rscript figure3H_NGS_final.R [full path to directory that contains all the RAVA output directories]

## **Alternatively, simply download this repository and run figure3H_NGS_final.R using the local path to this repository as the RAVA parent directory. You can either edit the hard-coded RAVA parent directory in the script, or provide it as the first command-line argument.**
