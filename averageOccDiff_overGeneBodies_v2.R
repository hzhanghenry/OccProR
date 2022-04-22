#!/usr/bin/env Rscript
library("optparse")

options = list(
	make_option(c("-c", "--control_bam_files"), 
		type="character", default=NULL,
		help="Bam file names of control samples, separated by comma 
		[e.g.: -c control_1.bam,control_2.bam,control_3.bam]"),
	make_option(c("-t", "--treat_bam_files"), 
		type="character", default=NULL,
		help="Bam file names of treatment samples, separated by comma 
		[e.g.: -c treatment_1.bam,treatment_2.bam,treatment_3.bam]"),
	make_option("--control_label", type="character",
		default="Control", help="Control label [default = %default]"),
	make_option("--treat_label", type="character", default="Treatment",
		help="Treatment label [default = %default]"),
	make_option(c("-l", "--min_length"), type="integer", default=50, 
		help="The smallest DNA fragment to be considered 
		[default = %default]"),
	make_option(c("-L", "--max_length"), type="integer", default=300, 
		help="The largest DNA fragment to be considered 
		[default = %default]"),
	make_option(c("-p", "--presorted_gene_list"), type="character", 
		default=NULL,
		help="Presorted list of genes (csv file)"),
	make_option("--color_scale_diff", type="double", default=1,
		help="Maximum color scale in difference heat map 
		[default = %default]"),
	make_option("--annotation_file", type="character", 
		default="sacCer3_annotations.csv",
		help="Name of file (csv format) containing the annotations 
		[default = %default]"),
	make_option("--normalization", type="character", default="razvan",
		help="Normalization methods, include 'razvan' and 'scaling'\
		[default = %default]"),
	make_option("--control_scale_factor", type="numeric", default=1,
		help="Scale factor for control samples. [default = %default]"),
	make_option("--treat_scale_factor", type="numeric", default=1,
		help="Scale factor for treatment samples. [default = %default]"),
	make_option("--output_dir", type="character", 
		default="Avg_Occ_overGeneBodies",
		help="Scale factor for treatment samples. [default = %default]")		
) 


opt_parser = OptionParser(option_list=options);
opt = parse_args(opt_parser);

if (is.null(opt$control_bam_files) | is.null(opt$treat_bam_files)){
  print_help(opt_parser);
  stop("Missing one or two set of bam file names.", call.=FALSE);
}

#*	Load the necessary R packages
#*	-----------------------------------------------
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(caTools)
  library(colorRamps)
  library(Rsamtools)
});

source("bamR_source_March_2021.R");




#* 	I.	Initialization (parsering of input parameteres)
#*
#*	=====================================================================

#*	Labels
#*	-------------------------------------------------------
control_label <- opt$control_label;
treat_label <- opt$treat_label;

#*	fragment length selection parameters:
#*	-------------------------------------------------------
min_length = opt$min_length ;
max_length = opt$max_length ;
control_scale_factor <- opt$control_scale_factor;
treat_scale_factor <- opt$treat_scale_factor;
normalize_method <- opt$normalization;

if(!normalize_method %in% c("razvan", "scaling"))
	stop("Please use 'razvan' or 'scaling' for normalize");


#*	Input data and default parameters
#*	-------------------------------------------------------
control_bam <- strsplit(opt$control_bam_files, ",")[[1]];
treat_bam <- strsplit(opt$treat_bam_files, ",")[[1]];

gene_list_file <- opt$presorted_gene_list;
annotation_file <- opt$annotation_file;
color_scale_diff <- opt$color_scale_diff;


#*	Default objects 
#*	-------------------------------------------------------
output_dir <- opt$output_dir;
sacCer3_chr <- getSacCer3ChromosomeRanges();
sacCer3_rDNA <- getSacCer3rRNARegion();
bam_param <- getScanPairedEndUnstrandedBamParameters();


#*	Sample name for output files
#*	-------------------------------------------------------
plot_title <- paste0(treat_label, "-", control_label);

if (! is.null(opt$presorted_gene_list)){
  sample_name <- paste(treat_label, "-", control_label, 
	".presorted", sep="");
} else {
  sample_name <- paste0(treat_label, "-", control_label);
}
file_name_base <- paste0(output_dir, 
	"/Avg_Occ_diff_over_gene_bodies.", sample_name, ".", 
	min_length, "_", max_length, "_", normalize_method);


#*	End of part I.
#*	*********************************************************************




#*	II.	Read bam files, remove the reads from rDNA, and keep the reads 
#*		with desired length
#*	=====================================================================

control_reads <- scanUnstrandedPairedEndBamFiles(control_bam, bam_param);
treat_reads <- scanUnstrandedPairedEndBamFiles(treat_bam, bam_param);

control_reads <- removeUnQualifiedReads(control_reads, sacCer3_rDNA,
					"chrM", min_length, max_length);
treat_reads <- removeUnQualifiedReads(treat_reads,sacCer3_rDNA,
					"chrM", min_length, max_length)

#*	End of part II.
#*	*********************************************************************




#*	III.	Calculate and normalize occupancy (coverage) 
#*	=====================================================================

scale_ctrl_Occ <- getReadsOccupancyAlongChromosome(control_reads,
					control_scale_factor, normalize=normalize_method);
scale_treat_Occ <- getReadsOccupancyAlongChromosome(treat_reads,
					treat_scale_factor, normalize=normalize_method);

#*	End of part III.
#*	*********************************************************************




#*	IV.	Compute average occupancy and the difference between control 
#*		and treatment samples.
#*	=====================================================================

sacCer3_annot <- loadAnnotationData(annotation_file, ",", 
		gene_list_file, sacCer3_chr);

avg_Occ_control <- averageOccOverGeneBody(Occupancy=scale_ctrl_Occ, 
				Reference=sacCer3_annot);
				
avg_Occ_treatment <- averageOccOverGeneBody(Occupancy=scale_treat_Occ, 
				Reference=sacCer3_annot);

if(normalize_method == "razvan") {
	avg_Occ_control[avg_Occ_control < 0] <- 0 ;
	avg_Occ_treatment[avg_Occ_treatment < 0] <- 0;
} else {
	avg_Occ_control <- log2(avg_Occ_control + 0.01) ;
	avg_Occ_treatment <- log2(avg_Occ_treatment + 0.01);
}
avg_Occ_diff <- avg_Occ_treatment - avg_Occ_control;

avg_Occ <- avg_Occ_diff;
smoothed_avg_Occ <- cbind(runmean(avg_Occ, 21, alg="C", endrule="mean"), 
	runmean(avg_Occ, 21, alg="C", endrule="mean"));
max_scale_diff <- color_scale_diff;

			
#* 	Save the average occupancy
#*	-----------------------------------------------------------
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE);

save(avg_Occ_diff, file=paste0(file_name_base, ".RData"));
write.csv(data.frame(Gene=names(avg_Occ_diff), 
	Average_Occ_diff=avg_Occ_diff), 
	file=paste0(file_name_base, ".csv"), row.names=FALSE)

#*	End of part IV.
#*	*********************************************************************




#*	V.	Heatmap plot 
#*	=====================================================================

pdf_file <- paste0(file_name_base, ".pdf");
pdf(pdf_file, width=5, height=8);

#*	heamap
#*	---------------------------------------------------------------
layout(matrix(c(1,2), ncol=2), widths=c(3,2));
par(mar=c(5,5,2,2))
image(c(1,2), 1:nrow(smoothed_avg_Occ), t(smoothed_avg_Occ), 
	col=matlab.like(102), ylim=c(nrow(smoothed_avg_Occ)+0.5,0.5), 
    breaks=c(min(-0.001-max_scale_diff, min(smoothed_avg_Occ)), 
		seq(-max_scale_diff, max_scale_diff, length.out = 101), 
			max(0.001+max_scale_diff,max(smoothed_avg_Occ))), 
     axes=FALSE, xlab="", useRaster=TRUE, cex.lab=1.4,
	 ylab="Average occupancy difference over gene bodies"
);
title(main=plot_title, cex.main=1);
box();

#*	color key scale
#*	---------------------------------------------------------------
par(mar=c(25,0,2,4))
image.scale(col=matlab.like(100), horiz=FALSE, xlab="", yaxt="n", 
	breaks=seq(-max_scale_diff, max_scale_diff, length.out=101));
axis(4, at=seq(-max_scale_diff, max_scale_diff, length.out=5), las=2)
box()
garbage <- dev.off();

#*	End of part V.
#*	*********************************************************************




#*	End of the script. 
#*	Last updated on March 29, 2021
#*	========================================================================
