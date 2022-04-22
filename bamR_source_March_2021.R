#*
#*		R Source for MNaseSeq Occupancy Calculation
#*
#*		Version 0.9
#*		Last updated on January 31, 2020
#*		Reference source: https://github.com/rchereji/bamR
#*
#*		Functions adopted/modified:
#*
#*		1.	getScanPairedEndUnstrandedBamParameters()
#*		2.	scanUnstrandedPairedEndBamFiles()
#*		3.	getSacCer3rRNARegion()
#*		4.	removeUnQualifiedReads()
#*		5.	getReadsOccupancyAlongChromosome()
#*		6.	getSacCer3ChromosomeRanges()
#*		7.	loadAnnotationData()
#*		8.	averageOccOverGeneBody()
#*		9.	image.scale()
#*
#*	_____________________________________________________________________
#*	<R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>



require(Rsamtools)


#*	=====================================================================
#*
#*		1.	getScanPairedEndUnstrandedBamParameters()
#*
#*		Purpose:	Get default parameters for scanBam(). This is for 
#*					scanning bam files with unstranded pair-end sequence 
#*					data and only will keep mapped read pairs which are 
#*					mapped on forward strand only.
#*			
#*		Argument:	None
#*		Return:		ScanBamParam object
#*
#*
getScanPairedEndUnstrandedBamParameters <- function()
{
	all_fields <- c("rname", "pos", "isize");
	param <- ScanBamParam(what=all_fields, 
				flag=scanBamFlag(isPaired=TRUE, isProperPair=TRUE, 
					isUnmappedQuery=FALSE, hasUnmappedMate=FALSE, 
					isMinusStrand=FALSE, isMateMinusStrand=TRUE,
					isNotPassingQualityControls=FALSE)
			 );
			 
	return (param);
}


#*	=====================================================================
#*
#*		2.	scanUnstrandedPairedEndBamFiles()
#*	
#*		Purpose:	Scan bam file(s) to GRanges Object. 
#*	
#*		Arguments:
#*	
#*			bam_files:	Name (and path) of one or more bam file(s). Bam 
#*						file(s) must be from unstranded and pair-end 
#*						sequence data.
#*			params:		A ScanBamParam object returned from calling to
#*						getScanPairedEndUnstrandedBamParameters().
#*
#*		Return:		GRanges list with all qualified reads and strand of 
#*					GRanges will be all '*'.
#*

scanUnstrandedPairedEndBamFiles <- function(bam_files=NULL, params=NULL)
{
	if(is.null(bam_files) | is.null(params))
		stop("Missing functon argument(s)");
		
	num_of_files <- length(bam_files);
	reads <- GRanges();	
	
	for (a_bam in 1:num_of_files)
	{
		bam <- scanBam(bam_files[a_bam], param=params)
  
		#*	Keep only the proper reads, with the length > 0
		gooddReads <- which(bam[[1]]$isize > 0);
		reads0 <- GRanges(seqnames=Rle(bam[[1]]$rname[gooddReads]),
					ranges = IRanges(start=bam[[1]]$pos[gooddReads], 
						width=bam[[1]]$isize[gooddReads]),
					strand = "*")
		rm(bam);
		reads <- c(reads, reads0);
	}

	return (reads);
}


#*	=====================================================================
#*
#*		3.	getSacCer3rRNARegion()
#*
#*		Purpose:	Get Saccharomyces cerevisiae rDNA gene range.
#*			
#*		Argument:	None
#*		Return:		A GRanges object with sacCer rDNA range.
#*
getSacCer3rRNARegion <- function()
{
	rDNA_region <- GRanges(seqnames = "chrXII",
					ranges = IRanges(start=451000, end=469000),
					strand = "*");
	return (rDNA_region);
}


#*	=====================================================================
#*
#*		4.	removeUnQualifiedReads()
#*
#*		Purpose:	Filter GRanges of reads scanned from bam file(s)
#*					to remove reads from yeast rDNA, mitochondrial
#*					DNA, and the reads out of desired lengths.
#*
#*		Arguments:
#*
#*			reads:			GRange list of reads to be filtered
#*			rDNA_region:	A GRanges object with rDNA gene range.
#*			chrM_name:		Character vector, mitochondrial chromosome 
#*							name, defaul: "chrM"
#*			min_length:		Positive integer, minimum length of reads
#*							to keep
#*			max_length:		Positive integer, maximum length of reads
#*							to keep
#*
#*		Return:		GRanges list of filtered reads
#*
removeUnQualifiedReads <- function(reads=NULL, rDNA_region=NULL,
			 chrM_name="chrM", min_length=50, max_length=300)
{
	if(is.null(reads) | is.null(rDNA_region))
		stop("GRanges for reads and rDNA region must be provided");
		
	#*	Discard the reads from rDNA, chrXII:451000-469000
	#*	-------------------------------------------------------

	rDNA_index <- overlapsAny(reads, rDNA_region, ignore.strand=TRUE);
	reads <- reads[!rDNA_index];

	#*	Discard the reads from the yeast mitochondrial DNA
	#*	-------------------------------------------------------
	if(length(grep("chrM", seqlevels(reads))) > 0) {
		reads <- reads[seqnames(reads) != chrM_name];
		reads <- dropSeqlevels(reads, chrM_name);
	}

	#*	Eliminate the reads that are shorter or longer than
	#*	desired length range
	#*	----------------------------------------------------------------
	read_length <- width(reads)
	is_good <- ((read_length >= min_length) & (read_length <= max_length));
	reads <- reads[is_good];
	
	return (reads);
}


#*	=====================================================================
#*
#*		5.	getReadsOccupancyAlongChromosome()
#*
#*		Purpose:	Calculate coverage of reads and normalize/scale it 
#*					if requested.
#*
#*		Arguments:
#*
#*			reads:			GRange list of reads to be filtered
#*			scale_factor:	Positive numeric, scaling factor for OCC.
#*							default: 1 and can be a customized such
#*							as caluclated from CPM method. Default is
#*							0 for using Razvan's original methods.
#*
#*		Return:		A RLE list presenting coverage of reads along each 
#*					chromosome
#*
#*		Last edited on February 12, 2020

getReadsOccupancyAlongChromosome <- function(reads, scale_factor=1, 
	normalize=c("razvan", "scaling", "ratio_to_target", "fold_to_raw"))
{
	if(is.null(reads)) stop("Reads are missing");
	
	raw_occ <- coverage(reads);
	chr_label <- seqlevels(raw_occ);
	num_of_chr <- length(chr_label);

	normalize <- tolower(normalize);
	coverage_weight = list();
	for(a_chr in chr_label) {	
		if(normalize == "razvan") {
			coverage_weight[[a_chr]]  <- 1/mean(raw_occ[[a_chr]]);
		} else if (normalize == "scaling") {
			coverage_weight[[a_chr]] <- scale_factor;
		} else {  stop("Unknown normalize method!"); }
	}

	scaled_occ <- coverage(reads, weight=coverage_weight);

	return (scaled_occ);
}


#*	=====================================================================
#*	
#*		6.	getSacCer3ChromosomeRanges()
#*
#*		Purpose:	Construct GRanges for all full length of  
#*					Saccharomy cescerevisiae chromosomes.
#*
#*		Argument: 	None
#*		Return:		GRanges object of full length of Saccharomy  
#*					cescerevisiae chromosomes.
#*
getSacCer3ChromosomeRanges <- function()
{
	#*	Construct GRanges for the entire chromosomes...
	#*	---------------------------------------------------------
	chr_length <- c(230218,  813184, 316620, 1531933,  576874, 
					270161, 1090940, 562643,  439888,  745751,
					666816, 1078177, 924431,  784333, 1091291,
					948066)
	names(chr_length) <- c('chrI','chrII','chrIII','chrIV','chrV',
			'chrVI','chrVII','chrVIII','chrIX','chrX','chrXI',
			'chrXII','chrXIII','chrXIV','chrXV','chrXVI');
			
	whole_chr <- GRanges(seqnames=names(chr_length),
                   IRanges(start=rep(1, length(chr_length)),
                           end=chr_length));
	
	return (whole_chr);
}


#*	=====================================================================
#*	
#*		7.	loadAnnotationData()
#*
#*		Purpose:	Read annotation file such as a CSV file derived from 
#*					MNaseSeq profiling showing Dyas ORF information.
#*
#*		Arguments:
#*
#*			annotation_file:	Character vector, name (and path) of 
#*								the annotation file to read. The file
#*								should be in csv or tab-delimited format.
#*			gene_list_file:		Character vector, name (and path) of file
#*								with gene list to filter annotations.The 
#*								format should be same as annotation file.
#*			delimitor:			Character, either ',' or '\t'.
#*			chr_ranges:			GRanges object with full length of all
#*								chromosomes to remove regions fall to
#*								outside of chromosome ranges.
#*
#*		Return:	GRanges object with all reference regions.
#*
loadAnnotationData <- function(annotation_file=NULL, delimitor=",",
			gene_list_file=NULL, chr_ranges=NULL)
{
	if(is.null(annotation_file) || is.null(chr_ranges))
		stop("Missing one or more arguments.");
	if(delimitor != "," && delimitor != "\t")
		stop("Unsupported file delimitor.");
	
	transcripts <- read.delim(annotation_file, sep=delimitor,
			header=TRUE, stringsAsFactors=FALSE);
	rownames(transcripts) <- transcripts$ORF;

	if (! is.null(gene_list_file)){
		presorted_genes = read.delim(gene_list_file, header=FALSE, 
				sep=delimitor, stringsAsFactors=FALSE);
		presorted_genes <- presorted_genes$V1
		index <- match(presorted_genes, rownames(transcripts));
		transcripts <- transcripts[index,];
	}

	ORF_Start <- transcripts$ORF_Start
	ORF_End <- transcripts$ORF_End
	Reference_chr <- transcripts$Chr
	Ref_strand <- transcripts$Strand
	ORF <- transcripts$ORF

	Watson = (Ref_strand == 1)
	left_edge = ORF_Start
	right_edge = ORF_End

	left_edge_crick = ORF_End
	right_edge_rrick = ORF_Start

	left_edge[!Watson] <- left_edge_crick[!Watson];
	right_edge[!Watson] <- right_edge_rrick[!Watson];

	annoations <- GRanges(seqnames=Reference_chr,
                IRanges(start=left_edge,end=right_edge),
                strand=Ref_strand,
                ORF=ORF)

	hits <- as.matrix(findOverlaps(annoations, chr_ranges, 
					type="within", ignore.strand=TRUE))
	index_to_keep <- hits[,"queryHits"]
	annoations <- annoations[index_to_keep];
	
	return (annoations);
}


#*	===========================================================================
#*	
#*		8.	averageOccOverGeneBody()
#*
#*		Purpose:	Calculate average coverage over gene body regions.
#*	
#*		Arguments:
#*
#*			Occupancy:           - The profile that needs to be aligned, 
#*							e.g. coverage/occupancy profile
#*			Reference:	GRanges with the windows centered on the 
#*								reference points
#*
#*	Output:
#*  	AvgOcc            - a vector with the average occupancy for each 
#*							gene body

averageOccOverGeneBody <- function(Occupancy, Reference)
{
	#*	Create Views with all the ReferenceGRanges
	#*	---------------------------------------------------
	chr_name <- unique(as.character(seqnames(Reference)));
	chr_views = Views(Occupancy[chr_name], as(Reference, 
					"IntegerRangesList")[chr_name])
	avg_profile_list <- viewMeans(chr_views);
	average_profile <- unlist(avg_profile_list);
  
	#*	Get the index of ReferenceGRanges, which were reorganized  
	#*	by as(ReferenceGRanges, "IntegerRangesList")
	#*	-------------------------------------------------------
	list_index <- split(1:length(Reference), 
			as.factor(seqnames(Reference)));
	index <- do.call("c", list_index);
  
	names(average_profile) <- index;
	average_profile <- average_profile[order(index)];
	names(average_profile) <- Reference$ORF;
	
	return(average_profile);
}


#*	===========================================================================
#*	
#*		9.	image.scale()
#*
#*		Purpose:	Drawing the image color scale 
#*
#*		Arguments:
#*
#*		Return: None.

image.scale <- function(z, zlim, col = heat.colors(12),
                    breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...)
{
	if(!missing(breaks)){
		if(length(breaks) != (length(col)+1))
		{stop("must have one more break than colour")}
	}
	if(missing(breaks) & !missing(zlim)){
		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
	}
	if(missing(breaks) & missing(zlim)){
		zlim <- range(z, na.rm=TRUE)
	
		#*	adds a bit to the range in both directions
		#*
		zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)
		zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
	}
	
	poly <- vector(mode="list", length(col));
	for(i in seq(poly)){
		poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
	}
	
	xaxt <- ifelse(horiz, "s", "n");
	yaxt <- ifelse(horiz, "n", "s");
	
	if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)};
	if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)};
	if(missing(xlim)) xlim=XLIM;
	if(missing(ylim)) ylim=YLIM;
	
	plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, 
				xaxs="i", yaxs="i", ...)  
	for(i in seq(poly)){
		if(horiz) polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA);
		if(!horiz) polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA);
	}
}

