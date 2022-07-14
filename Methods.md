# Table of Contents
- Data analysis
	- Software tools
	- Data files
	- Specific examples of how to apply the workflow
		- G4 enrichment at ageing clock CpGs
		- G4 enrichment at protein binding sites
		- Obtaining global CpG sites

# Data analysis
This file contains a general overview of how to implement the colocalisation analysis from the paper. Furthermore, there are more specific examples of how the data shown was produced. Data and software tools used are given below for reproducing the results. 
Utility functions created for this workflow are within `functions.R`. 
All genomic coordinates were converted to the hg19 (human) or mm10 (mouse) genome build before analysis using the `liftOver` function of the `rtracklayer` package using the appropriate chain files.
For better reproducibility, no exact paths or file names are given for the general workflow. Instead `$` is used in combination with informative placeholder names. 
## Software tools

- R (4.0.3)
- dplyr (1.0.5)
- (rjson) (0.2.20)
- GenomicRanges (1.42.0)
- regioneR (1.22.0)
- IRanges (2.24.1)
- Biostrings (2.58.0)
- tidyverse (1.3.0)
- scales (1.1.1)
- rtracklayer (1.49.5)
- stringr (1.4.0)

## Data files
Geo repositories are given in parentheses. 
**Required**
- G4 data:
	- G4-seq: 
		- Human G4-seq (K+) (GSM3003539)
		- Mice G4-seq (K+) (GSM3003547)
	- BG4-ChIP: 
		- NHEK (GSE76688)
		- HaCaT (GSE76688)
		- K562 (GSE107690)
- CpG Aging clocks:
	- Horvath: [Supplementary Data](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115)
	- Levine: [Supplementary Data](https://www.aging-us.com/article/101414/text)
	- Stubbs / Meer / Petkovich: [Supplementary Data](https://elifesciences.org/articles/40675)
- ChIP-seq TET, DNMT
	- TET1 (mouse ESCs) (GSE24841)
	- TET2 (mouse ESCs) (GSE41720)
	- TET1 (human ESCs) (GSE99346)
	- DNMT1 (K562) (GSE92213)
	- DNMT1 (HepG2) (GSE170872)
	- DNMT3b (HepG2) (GSE95953)


## Examples of how to apply the colocalisation workflow.
### G4 enrichment at ageing clock CpGs. 
Following this example will output a table containing the type of data that was used to generate the published plots. Use the specified datasets to reproduce the exact results. The first example shows how to get the G4 enrichment at ageing clock CpGs relative to random sites and global CpG sites in the genome. 

1. Import functions and libraries. 
```
library(dplyr)
library(rjson)
source("functions.R")
```

2. Load G4 data and catenate both strands.
```
G4_locs <- load_G4_data(G4_file_plus = $G4_file_strand_1, G4_file_minus=$G4_file_strand_2)
```

3. Load CpG data, converting from `.csv` to `.bed` format.
```
CpG_locs <- load_CpG_data(CpG_file=$CpGs_ageing_clock_data_file)
```

4. Lift over genomic coordinates using the appropriate chain file. This will make different genome assemblies compatible. 
```
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=$chain_file)
```

5. Get number of all CpGs in the genome and their overlaps with the G4 dataset. The output file will contain data on how many global CpGs overlap with G4s. This will be used later by other functions to obtain the enrichment value of G4s relative to global CpGs instead of relative to random regions in the genome. The `$all_CpGs_file.bed` containing the coordinates of all CpGs in the genome needs to be generated beforehand (see examples below).  The file `$genome_file` is a text file characterising valid regions of the genome (an [example](link_to_genome_file) is included in this repository). This step has a long execution time.
```
results_all_CpGs <- global_CpG_overlap(all_CpGs_file=$all_CpGs_file.bed,
	G4_locs=G4_locs,genome_file=$genome_file, reduce=F, window_sizes=window_sizes)
write.table(results_all_CpGs, file = "out/all_CpGs_props/all_CpGs_props_all_G4s_mouse.csv", 
	sep=";", row.names=F)
```

6. Run the colocalisation analysis for mutliple window sizes and obtain results on statistical tests. The function lets you specify a query and search set (make sure to shuffle the smaller dataset for increased performance by setting parameter `shuffle="q"` for shuffling the query or `shuffle="s"` for search set). `window sizes` can be an arbitrary vector of different window sizes to test. `random_trials` specifies the number of trials of random shuffling.  
```
results <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file=$genome_file,
	window_sizes=window_sizes, all_CpGs_props=results_all_CpGs, random_trials = 30, 
	shuffle = "q")
```

7. Write result statistics to file.
```
write.table(results$stats, file=$output_file, sep=";", row.names = F)
```

### G4 enrichment at protein binding sites
This is the adapted workflow of how ChIP-seq data of epigenetic regulators (e.g. TET1) was colocalised with G4s to obtain potential enrichment of G4s within protein binding sites. This example can be adapted to the other ChIP-seq datasets (TET1/2 mouse, DNMT1, DNMT3b) to obtain the other published results. 

1. Load ChIP-seq data and bring in the correct format. Do some filtering to remove X, Y and mitochondrial chromosomes.
```
TET_data <- read.csv($ChIP-seq_data_file_for_TET, sep="\t", header=F)
names(TET_data) <- c("chromosome", "start", "end", ".", "strand")
TET_data <- TET_data %>% select(chromosome, start, end) %>% filter(chromosome %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))
```

2. Load G4s (autosomes only). Filter out mitochondrial G4s. 
```
G4_locs <- load_G4_data(G4_file_plus=$G4_file_strand_1, G4_file_minus=$G4_file_strand_2, autosomes=T)
G4_locs <- G4_locs %>% filter(!(chromosome=="chrM"))
```

3. Overlap TET peaks and G4s. Use `analyse_overlap()` function if not applying window sizes (difference to CpG case). 
```
G4s_in_TETs <- analyse_overlap(query=G4_locs, search_set=TET_data, genome_file=$genome_file, shuffle="s", random_trials = 30)
```

4. Construct output dataframe to write to file.
```
stats_G4s <- list("direction"="G4s_in_TETs", "num_overlap_subset"=G4s_in_TETs$num_overlap_subset, "num_overlap_bp"=G4s_in_TETs$num_overlap_bases, "p_value_MC_subset"=G4s_in_TETs$p_value_MC_subset, "FE_MC_subset"=G4s_in_TETs$FE_MC_subset, "p_value_MC_bp"=G4s_in_TETs$p_value_MC_bp, "FE_MC_bp"=G4s_in_TETs$FE_MC_bp)
stats_G4s <- data.frame(stats_G4s)
```

5. Write results to file.
```
write.table(stats_G4s, file="out/TET_analysis/mouse_TET1_C_all_G4_stats.csv", sep=";", row.names = F)
```

### Obtaining global CpG sites 
An example of how the sites of all CpG sites across the genome were obtained (here from mouse genome assembly mm10). 
```
library(tidyverse)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
setwd("D:/proj_epigen/aging_clocks")

# load genome
genome <- BSgenome.Mmusculus.UCSC.mm10
chromosomes <-seqnames(genome)[0:19]
pattern <- "CG"

# get number of CpGs in all chromosomes and find genomic coordinates
CpG_occurence <- vector(mode='list', length=23)
names(CpG_occurence) <- c(chromosomes, 'total')
total <- 0
CpG_coordinates <- tibble()
sequences <- tibble()
seqs_chr <- vector(mode='list', length=22)

#apply window to CpG locations and save in extended CpG coordinates tibble
extension_bases <- 500
first <- 1

# loop through chromosomes to get CpGs for each one
for (chr_name in chromosomes) {
	current_chr <- genome[[chr_name]]
	#count CpG occurence and save in named list
	n <- countPattern(pattern, current_chr)
	CpG_occurence[chr_name] <- n
	total <- total + n
	
	#match pattern to get genomic coordinates of each CpG
	seqs <- matchPattern(pattern, current_chr) #gets XStringView object of coordinates
	seqs <- as(seqs, "IRanges") #convert object to IRanges object
	
	# get start and end coordinates, save in temporary tibble, then combine with
	# data from other chromosomes
	starts <- start(seqs) - extension_bases
	ends <- end(seqs) + extension_bases
	cur_chr_tibble <- tibble(chr=chr_name, start=starts, end=ends)
	CpG_coordinates <- rbind(CpG_coordinates, cur_chr_tibble)
	a(sequence, identifiers, "sequences/all_CpG_mouse_seq.fasta", open="a", as.string=T)
}
CpG_occurence[[23]] <- total

# merge coordinates before saving
r <- makeGRangesFromDataFrame(CpG_coordinates, keep.extra.columns = T)
r <- reduce(r)
starts <- start(ranges(r))
ends <- end(ranges(r))
chroms <- as.data.frame(seqnames(r))
merged_coor <- cbind(chroms, starts, ends)
write.table(merged_coor, file="data/CpG_lists/all_CpGs_hg19_ws_50_merged.bed", sep="\t", col.names = F, row.names = F, quote = F)
```


