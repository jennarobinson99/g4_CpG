# This file contains all functions needed for overlap analysis. 


overlap <- function(query=NULL, search_set=NULL, reduce=F){
  ### A and B should be Ranges objects which can be reduced and overlapped (subset overlap, bedtools -wa mode)
  library(GenomicRanges)
  ranges_A <- makeGRangesFromDataFrame(query, keep.extra.columns = T)
  ranges_B <- makeGRangesFromDataFrame(search_set, keep.extra.columns = T)
  
  # overlap based on bases 
  overlap_bases <- GenomicRanges::intersect(ranges_A, ranges_B, ignore.strand=T)
  # subset overlap mode: giving subset of query that is overlapping
  overlap_ranges <- subsetByOverlaps(ranges_A, overlap_bases)
  
  if (reduce){
    overlap_ranges <- GenomicRanges::reduce(overlap_ranges) # merges potential duplicates
  }
  
  # convert to data frame, return resulting dataframe
  overlap_coor <- GRanges_to_df(ranges=overlap_ranges)
  overlap_bases <- GRanges_to_df(ranges=overlap_bases)
  
  return(list("overlapping_bases"=overlap_bases, "overlap_subset"=overlap_coor))
}


overlap_bases <- function(query=NULL, search_set=NULL){
  ### Overlaps A and B and extracts all genomic regions present in both 
  
  library(GenomicRanges)
  ranges_A <- makeGRangesFromDataFrame(query, keep.extra.columns = T)
  ranges_B <- makeGRangesFromDataFrame(search_set, keep.extra.columns = T)
  
  intersect <- intersect(ranges_A, ranges_B)
  intersect_df <- GRanges_to_df(ranges=intersect)
  return(intersect_df)
}


analyse_overlap <- function(query = NULL, search_set = NULL, genome_file=NULL, random_trials=5, all_CpGs_props=NULL, shuffle="q", mask=NULL) {
  ### Function overlaps two files and analyses the overlap, gives out bp and CpG-based p-values, enrichment
  
  library(regioneR)
  library(IRanges)
  library(Biostrings)
  library(stats)
  
  # Overlap the file and get query entries that are in the search set
  overlap_results <- overlap(query=query, search_set = search_set)
  overlap_bases <- overlap_results[["overlapping_bases"]]
  num_overlap_bases <- sum(overlap_bases$end - overlap_bases$start + 1)
  overlap_coor <- overlap_results[["overlap_subset"]]
  overlap_coor <- query %>% filter(start %in% overlap_coor$start)
  num_overlaps <- dim(overlap_coor)[1]
  if(num_overlaps != 0){
    overlap_ranges <- makeGRangesFromDataFrame(overlap_coor)
  }
  
  # Load and process genome file
  genome <- read.csv(genome_file, header = F,sep="\t", col.names = c("chromosome", "length", "unknown"))
  genome <- cbind.data.frame(genome$chromosome, 1, genome$length)
  names(genome) <- c("chromosome", "start", "end")
  
  # Random shuffling overlap (Monte Carlo method)
  rand_overlaps <- vector()
  rand_overlap_bases <- vector()
  for(i in 1:random_trials){
    set.seed(i)
    print("Shuffling trial: ")
    print(i)
    if(shuffle=="q"){
      shuffled_coor <- GRanges_to_df(randomizeRegions(query, genome=genome, mask = mask))
      shuffled_overlap <- overlap(query=shuffled_coor, search_set = search_set)
      rand_overlaps[i] <- dim(shuffled_overlap[["overlap_subset"]])[1]
      rand_overlap_bases[i] <- sum(shuffled_overlap[["overlapping_bases"]]$end -
                                     shuffled_overlap[["overlapping_bases"]]$start + 1)
    } else if(shuffle=="s"){
      shuffled_coor <- GRanges_to_df(randomizeRegions(search_set, genome=genome, mask=mask))
      shuffled_overlap <- overlap(query=query, search_set = shuffled_coor)
      rand_overlaps[i] <- dim(shuffled_overlap[["overlap_subset"]])[1]
      rand_overlap_bases[i] <- sum(shuffled_overlap[["overlapping_bases"]]$end -
                                     shuffled_overlap[["overlapping_bases"]]$start + 1)
    } else{
      print("Invalid argument for shuffle.")
    }
  }
  # calculate p-value based on MC method CpG_based
  p_value_MC_subset <- p_value_MC(num_overlaps, rand_overlaps)
  # Calculate p-value based on base overlap
  p_value_MC_bp <- p_value_MC(num_overlap_bases, rand_overlap_bases)
  
  
  # Calculate fold enrichment based on MC method for subset overlap
  average_rand_overlap <- mean(rand_overlaps)
  if(num_overlaps == 0){
    print("No overlaps found.")
    enrichment_MC_subset <- 0
  }
  else {
    enrichment_MC_subset <- num_overlaps / (average_rand_overlap + 1)
  }
  
  # Calculate fold enrichment based on MC method for base overlap
  average_rand_base_overlap <- mean(rand_overlap_bases)
  if(num_overlap_bases == 0){
    print("No overlaps found for base overlap.")
    enrichment_MC_bp <- 0
  }
  else {
    enrichment_MC_bp <- num_overlap_bases / (average_rand_base_overlap + 1)
  }
  
  # Calculate p-value and enrichment with analytical method: against the whole genome according to number of basepairs (using hypergeometric dist)
  if(!is.null(all_CpGs_props)){
    num_all_bases <- all_CpGs_props[["num_bases"]]
    num_all_bases_overlap <- all_CpGs_props[["num_overlaps_bases"]]
    num_query_bases <- sum(query$end - query$start + 1)
    num_overlap_bases <- sum(overlap_bases$end - overlap_bases$start + 1)
    p_value_hyper_bp <- phyper(q = num_overlap_bases, m = num_all_bases_overlap, n = (num_all_bases-num_all_bases_overlap), k = num_query_bases, lower.tail = F)
    expected_overlap_bases <- num_all_bases_overlap/num_all_bases * num_query_bases
    enrichment_hyper_bp <- num_overlap_bases / expected_overlap_bases
  } else{
    p_value_hyper_bp <- NA
    enrichment_hyper_bp <- NA
  }
  
  # Calculate p-value and enrichment with analytical method: against number of CpGs genome-wide according to number of CpGs (NOT basepairs)
  if(!is.null(all_CpGs_props)){
    num_all_CpGs <- all_CpGs_props[[1]]
    num_all_CpGs_overlap <- all_CpGs_props[[2]]
    num_AC_CpGs <- dim(query)[1]
    num_AC_CpGs_overlap <- dim(overlap_coor)[1]
    p_value_hyper_CpGs <- phyper(q=num_AC_CpGs_overlap, m=num_all_CpGs_overlap, 
                           n=(num_all_CpGs - num_all_CpGs_overlap), 
                           k=num_AC_CpGs, lower.tail=F)
    expected_CpGs <- num_all_CpGs_overlap/num_all_CpGs * num_AC_CpGs
    enrichment_hyper_CpGs<- num_AC_CpGs_overlap/expected_CpGs
  } else {
    p_value_hyper_CpGs <- NA
    enrichment_hyper_CpGs <- NA
  }
  
  # save infos on how many positive / negative correlation; CGI context; chromatin accessibility
  if("coefficient" %in% names(overlap_coor)){
    num_pos <- dim(overlap_coor %>% filter(coefficient>0))[1] 
    num_neg <- dim(overlap_coor %>% filter(coefficient<0))[1]
    percentage_pos <- round(num_pos/dim(overlap_coor)[1], digits=4)*100
    percentage_neg <- round(num_neg/dim(overlap_coor)[1], digits=4)*100
  }else{
    num_pos <- "N/A"
    percentage_pos <- "N/A"
  }
  if("CGI_context" %in% names(overlap_coor)){
    num_in_CGI <- dim(overlap_coor %>% filter(CGI_context=="inside"))[1]
    num_outside_CGI <- dim(overlap_coor %>% filter(CGI_context=="outside"))[1]
    percentage_in_CGI <- round(num_in_CGI/dim(overlap_coor)[1], digits=4)*100
    percentage_outside_CGI <- round(num_outside_CGI/dim(overlap_coor)[1], digits=4)*100
  } else {
    num_in_CGI <- "N/A"
    percentage_in_CGI <- "N/A"
  }
  if("chromatin" %in% names(overlap_coor)){
    num_open <- dim(overlap_coor %>% filter(chromatin=="open"))[1]
    num_closed <- dim(overlap_coor %>% filter(chromatin=="closed"))[1]
    percentage_open <- round(num_open/dim(overlap_coor)[1], digits=4) * 100
    percentage_closed <- round(num_closed/dim(overlap_coor)[1], digits=4) * 100
  } else{
    num_open <- "N/A"
    percentage_open <- "N/A"
  }
  return(list("overlap_coordinates" = overlap_coor, "num_overlap_subset" = num_overlaps,
              "num_overlap_bases" = num_overlap_bases, 
              "p_value_MC_subset" = p_value_MC_subset, "FE_MC_subset"=enrichment_MC_subset,
              "p_value_MC_bp" = p_value_MC_bp, "FE_MC_bp" = enrichment_MC_bp,
              "p_value_hyper_CpGs" = p_value_hyper_CpGs, "FE_hyper_CpGs" = enrichment_hyper_CpGs,
              "p_value_hyper_bp" = p_value_hyper_bp, "FE_hyper_bp" = enrichment_hyper_bp,
              "num_pos" = num_pos, "percentage_pos" = percentage_pos, 
              "num_in_CGI" = num_in_CGI,"percentage_in_CGI" = percentage_in_CGI, 
              "num_open" = num_open, "percentage_open" = percentage_open))
}


analyse_window_size <- function(query=NULL, search_set=NULL, genome_file=NULL, window_sizes=NULL, all_CpGs_props=NULL, shuffle="q", random_trials=5, mask=NULL){
  ### Function loops through all given window sizes and calculates statistical metrics about distribution and enrichment for all. Outputs dataframe containing all results and produces graphs. 
  
  library(dplyr)
  i <- 1
  overlap_coor <- list()
  num_overlap_subset <- numeric()
  num_overlap_bases <- numeric()
  p_value_MC_subset <- numeric()
  FE_MC_subset <- numeric()
  p_value_MC_bp <- numeric()
  FE_MC_bp <- numeric()
  p_value_hyper_CpGs <- numeric()
  FE_hyper_CpGs <- numeric()
  p_value_hyper_bp <- numeric()
  FE_hyper_bp <- numeric()
  num_pos <- numeric()
  num_neg <- numeric()
  percentage_neg <- numeric()
  percentage_pos <- numeric()
  num_in_CGI <- numeric()
  num_outside_CGI <- numeric()
  percentage_outside_CGI <- numeric()
  percentage_in_CGI <- numeric()
  num_open <- numeric()
  percentage_open <- numeric()
  
  for(window_size in window_sizes){
    print("Analysing window size: ")
    print(window_size)
    #apply window size
    query_ext <- query %>% mutate(start=start-round(window_size/2), end=end+round(window_size/2)-1)
    # select correct all_CpGs_props entry
    ws <- window_size
    if(!is.null(all_CpGs_props)){
      curr_all_CpGs_props <- all_CpGs_props %>% filter(window_size==ws)
    }
    #get results for current window size
    temp_results <- analyse_overlap(query=query_ext, search_set = search_set, 
                                    genome_file=genome_file, all_CpGs_props=curr_all_CpGs_props, 
                                    shuffle=shuffle, random_trials = random_trials, mask=mask)
    
    # save coordinates of overlapping G4s and CpGs and augment with coefficient from input
    overlap_coor[i] <- temp_results["overlap_coordinates"]
    # save statistics of results
    num_overlap_subset[i] <- temp_results[["num_overlap_subset"]]
    num_overlap_bases[i] <- temp_results[["num_overlap_bases"]]
    p_value_MC_subset[i] <- temp_results[["p_value_MC_subset"]]
    FE_MC_subset[i] <- temp_results[["FE_MC_subset"]]
    p_value_MC_bp[i] <- temp_results[["p_value_MC_bp"]]
    FE_MC_bp[i] <- temp_results[["FE_MC_bp"]]
    p_value_hyper_CpGs[i] <- temp_results[["p_value_hyper_CpGs"]]
    FE_hyper_CpGs[i] <- temp_results[["FE_hyper_CpGs"]]
    p_value_hyper_bp[i] <- temp_results[["p_value_hyper_bp"]]
    FE_hyper_bp[i] <- temp_results[["FE_hyper_bp"]]
    num_pos[i] <- temp_results[["num_pos"]]
    percentage_pos[i] <- temp_results[["percentage_pos"]]
    num_in_CGI[i] <- temp_results[["num_in_CGI"]]
    percentage_in_CGI[i] <- temp_results[["percentage_in_CGI"]]
    num_open[i] <- temp_results[["num_open"]]
    percentage_open[i] <- temp_results[["percentage_open"]]
    i <- i + 1
  }
  
  stats <- cbind.data.frame(window_sizes, as.numeric(num_overlap_subset), as.numeric(num_overlap_bases),
                            as.numeric(p_value_MC_subset), as.numeric(FE_MC_subset), 
                            as.numeric(p_value_MC_bp), as.numeric(FE_MC_bp),
                            as.numeric(p_value_hyper_CpGs),as.numeric(FE_hyper_CpGs),
                            as.numeric(p_value_hyper_bp),as.numeric(FE_hyper_bp),
                            as.numeric(num_pos), as.numeric(percentage_pos), 
                            as.numeric(num_in_CGI), as.numeric(percentage_in_CGI), 
                            as.numeric(num_open), as.numeric(percentage_open))
  names(stats) <- c("window_size", "num_overlap_subset", "num_overlap_bases", "p_value_MC_subset",
                    "FE_MC_subset","p_value_MC_bp", "FE_MC_bp", "p_value_hyper_CpGs", "FE_hyper_CpGs",
                    "p_value_hyper_bp", "FE_hyper_bp", "num_pos", "percentage_pos", 
                    "num_in_CGI", "percentage_in_CGI", "num_open", "percentage_open")
  return(list("overlap_coordinates" = overlap_coor, "stats" = stats))
}


plot_results <- function(overlap_results=NULL, query_name="CpGs", search_set_name="G4s", 
                         window_size=NULL, stats=NULL, 
                         figure_name_enrichment="out/enrichment.pdf",
                         figure_name_p_value="out/p_value.pdf",
                         figure_name_CGI_dist="out/CpG_CGI_context_dist.pdf",
                         figure_name_coeff_dist="out/CpG_coeff_dist.pdf", 
                         figure_name_chromatin_dist="out/CpG_chromatin_dist.pdf", 
                         legend_position_fe = c(0.8,0.9), legend_position_p = c(0.8, 0.5)) {
  ### Function visualises the results of overlap by plotting enrichment vs. window size,
  ### CpG distribution according to CGI context and coefficient
  
  library(tidyverse)
  library(scales)
  if(!is.null(overlap_results)){
    stats <- overlap_results[["stats"]]
    overlap_coor <- overlap_results[["overlap_coordinates"]]
    idx <- match(window_size, stats$window_size)
    overlap_ws <- overlap_coor[[idx]]
  } else if(!is.null(stats)){
    stats <- stats
  }
  
  
  # construct data frame containing all fold enrichment values
  enrichment_df_1 <- cbind.data.frame(stats$window_size, stats$FE_MC_subset, "Counts-based MC test")
  enrichment_df_2 <- cbind.data.frame(stats$window_size, stats$FE_MC_bp, "Base pair-based MC test")
  enrichment_df_3 <- cbind.data.frame(stats$window_size, stats$FE_hyper_CpGs, "CpG-based analytical test")
  names(enrichment_df_1) <- c("window_size", "FE", "Test")
  names(enrichment_df_2) <- c("window_size", "FE", "Test")
  names(enrichment_df_3) <- c("window_size", "FE", "Test")
  enrichment_df <- rbind(enrichment_df_1, enrichment_df_2, enrichment_df_3)
  
  # construct p-value data frame
  p_value_df_1 <- cbind.data.frame(stats$window_size, stats$p_value_MC_subset, "Counts-based MC test")
  p_value_df_2 <- cbind.data.frame(stats$window_size, stats$p_value_MC_bp, "Base pair-based MC test")
  p_value_df_3 <- cbind.data.frame(stats$window_size, stats$p_value_hyper_CpGs, "CpG-based analytical test")
  names(p_value_df_1) <- c("window_size", "p_value", "Test")
  names(p_value_df_2) <- c("window_size", "p_value", "Test")
  names(p_value_df_3) <- c("window_size", "p_value", "Test")
  p_value_df <- rbind(p_value_df_1, p_value_df_2, p_value_df_3)
  p_value_df_zero <- p_value_df %>% filter(p_value < 1e-50) %>% mutate(p_value=0)
  p_value_df <- rbind((p_value_df %>% filter(p_value > 1e-50)), p_value_df_zero) %>% arrange(Test, window_size)
  
  #plot enrichment
  ylabel <- paste("Fold enrichment of", query_name, "in", search_set_name, sep=" ")
  xlabel <- paste("Window size in bp")
  # plot G4s
  ggplot(enrichment_df) +
    geom_line(aes(x = window_size, y = FE, colour=Test)) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_hline(yintercept = 1) +
    ylab(ylabel) + 
    xlab(xlabel) + 
    theme(text = element_text(size=18), legend.position = legend_position_fe)
  ggsave(figure_name_enrichment)
  
  #plot p_value
  ylabel <- paste("P-value of", query_name, "in", search_set_name, sep=" ")
  xlabel <- paste("Window size in bp")
  # plot G4s
  ggplot(p_value_df) +
    geom_line(aes(x = window_size, y = p_value, colour=Test), size=1) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_hline(yintercept = 0.05) +
    ylab(ylabel) + 
    xlab(xlabel) +
    theme(text = element_text(size=18), legend.position = legend_position_p)
  ggsave(figure_name_p_value)
  
  #plot coefficient distribution and print stats
  overlap_pos <- overlap_ws %>% filter(coefficient > 0) %>% mutate(sign="positive")
  overlap_neg <- overlap_ws %>% filter(coefficient < 0) %>% mutate(sign="negative")
  overlap_signed <- rbind(overlap_pos, overlap_neg) %>% arrange(chromosome, start)
  ggplot(overlap_signed) +
    geom_bar(aes(x=chromosome)) + 
    facet_grid(sign~.)
  ggsave(figure_name_coeff_dist) 
  
  #plot CGI context distribution
  if(NA %in% stats$num_in_CGI){
    print("No CGI information given. CpG distribution according to CGI context is not plotted. ")
  }else{
    ggplot(overlap_ws) +
      geom_bar(aes(x=chromosome)) + 
      facet_grid(CGI_context~.)
    ggsave(figure_name_CGI_dist) 
  }
  
  #plot chromatin accessibility distribution
  if(NA %in% stats$num_open){
    print("No chromatin accessibility information given. CpG distribution according to this factor is not plotted. ")
  } else{
    ggplot(overlap_ws) + 
      geom_bar(aes(x=chromosome)) + 
      facet_grid(chromatin~.)
    ggsave(figure_name_chromatin_dist)
  }
  
  # print stats
  print(stats)
}

plot_grouped_results <- function(overlap_results=NULL, query_name="CpGs", search_set_name="G4s", 
                         window_size=NULL, stats_1=NULL, stats_2=NULL, group_1="", group_2="", 
                         figure_name_enrichment="out/enrichment.pdf",
                         figure_name_p_value="out/p_value.pdf",
                         figure_name_CGI_dist="out/CpG_CGI_context_dist.pdf",
                         figure_name_coeff_dist="out/CpG_coeff_dist.pdf", 
                         figure_name_chromatin_dist="out/CpG_chromatin_dist.pdf", 
                         legend_position_fe = c(0.8,0.9), legend_position_p = c(0.8, 0.5)) {
  ### Function visualises the results of overlap by plotting enrichment vs. window size,
  ### CpG distribution according to CGI context and coefficient
  
  library(tidyverse)
  library(scales)
  if(!is.null(overlap_results)){
    stats <- overlap_results[["stats"]]
    overlap_coor <- overlap_results[["overlap_coordinates"]]
    idx <- match(window_size, stats$window_size)
    overlap_ws <- overlap_coor[[idx]]
  } else if(!is.null(stats_1)){
    stats_1 <- mutate(stats_1, group=group_1)
    stats_2 <- mutate(stats_2, group=group_2)
    stats <- rbind(stats_1, stats_2)
  }

  
  # construct data frame containing all fold enrichment values
  enrichment_df_1 <- cbind.data.frame(stats$window_size, stats$FE_MC_subset, stats$group, "Counts MC")
  enrichment_df_2 <- cbind.data.frame(stats$window_size, stats$FE_MC_bp, stats$group,"Base pair MC")
  enrichment_df_3 <- cbind.data.frame(stats$window_size, stats$FE_hyper_CpGs,stats$group, "CpG analytical ")
  names(enrichment_df_1) <- c("window_size", "FE","Group", "Test")
  names(enrichment_df_2) <- c("window_size", "FE","Group", "Test")
  names(enrichment_df_3) <- c("window_size", "FE", "Group","Test")
  enrichment_df <- rbind(enrichment_df_1, enrichment_df_2, enrichment_df_3)
  
  # construct p-value data frame
  p_value_df_1 <- cbind.data.frame(stats$window_size, stats$p_value_MC_subset,stats$group, "Counts MC")
  p_value_df_2 <- cbind.data.frame(stats$window_size, stats$p_value_MC_bp,stats$group, "Base pair MC")
  p_value_df_3 <- cbind.data.frame(stats$window_size, stats$p_value_hyper_CpGs,stats$group, "CpG analytical")
  names(p_value_df_1) <- c("window_size", "p_value","Group", "Test")
  names(p_value_df_2) <- c("window_size", "p_value","Group", "Test")
  names(p_value_df_3) <- c("window_size", "p_value","Group", "Test")
  p_value_df <- rbind(p_value_df_1, p_value_df_2, p_value_df_3)
  p_value_df_zero <- p_value_df %>% filter(p_value < 1e-50) %>% mutate(p_value=0)
  p_value_df <- rbind((p_value_df %>% filter(p_value > 1e-50)), p_value_df_zero) %>% arrange(Test, window_size)
  
  #plot enrichment
  ylabel <- paste("Fold enrichment of", query_name, "in", search_set_name, sep=" ")
  xlabel <- paste("Window size in bp")
  # plot G4s
  ggplot(enrichment_df) +
    geom_line(aes(x = window_size, y = FE, colour=Test, linetype=Group), size=1) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_hline(yintercept = 1) +
    ylab(ylabel) + 
    xlab(xlabel) + 
    theme(text = element_text(size=25), legend.position = legend_position_fe)
  ggsave(figure_name_enrichment)
  
  #plot p_value
  ylabel <- paste("P-value of", query_name, "in", search_set_name, sep=" ")
  xlabel <- paste("Window size in bp")
  # plot G4s
  ggplot(p_value_df) +
    geom_line(aes(x = window_size, y = p_value, colour=Test, linetype=Group), size=1) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_hline(yintercept = 0.05) +
    ylab(ylabel) + 
    xlab(xlabel) +
    theme(text = element_text(size=25), legend.position = legend_position_p)
  ggsave(figure_name_p_value)
}


global_CpG_overlap <- function(all_CpGs_file=NULL, G4_locs=NULL, genome_file=genome_file, random_trials=3, reduce=F, window_sizes=NULL){
  ### Function that loads and overlaps all CpGs in the genome and gives out overlap results
  
  library(GenomicRanges)
  library(dplyr)
  
  # load all CpGs and merge / reduce entries
  all_CpGs <- read.csv(all_CpGs_file, header=F, sep="\t")
  names(all_CpGs) <- c("chromosome", "start", "end")
  
  all_CpGs_props <- data.frame()
  
  for(window_size in window_sizes){
    print(window_size)
    all_CpGs_temp <- all_CpGs %>% mutate(start=start-round(window_size/2), end=end+round(window_size/2))
    num_bases <- sum(all_CpGs_temp$end - all_CpGs_temp$start + 1)
    overlap_results <- overlap(query = all_CpGs_temp, search_set=G4_locs, reduce = reduce)
    overlap_coor <- overlap_results[["overlap_subset"]]
    overlap_bases <- overlap_results[["overlapping_bases"]]
    num_overlap_bases <- sum(overlap_bases$end - overlap_bases$start + 1)
    # overlap CpGs with G4s 
    # (crashes if run!) results <- analyse_overlap(query=all_CpGs, search_set=G4_locs, genome_file=genome_file, random_trials=random_trials)
    all_CpGs_props_temp <- data.frame(dim(all_CpGs)[1], dim(overlap_coor)[1], 
                                      num_bases, num_overlap_bases, window_size)
    all_CpGs_props <- rbind(all_CpGs_props, all_CpGs_props_temp)
  }
  
  names(all_CpGs_props) <- c("num_CpGs", "num_overlaps_CpGs","num_bases", "num_overlaps_bases", "window_size")
  return(all_CpGs_props)
}


load_CpG_data <- function(CpG_file=NULL){
  ### Function loads csv data into dataframe and returns it
  
  # load data
  CpG_locs <- read.csv(CpG_file, header = TRUE, sep = ";")
  # get chromosomes vector in .bed file format
  chr_list <- list(CpG_locs[,"Chr"])
  chromosomes <- lapply(chr_list, chr_string)
  # get start end end positions
  starts <-  CpG_locs[, "MapInfo"]
  ends <-  CpG_locs[, "MapInfo"] + 1
  coefficients <- CpG_locs[,"CoefficientTraining"]
  # create data frame with .bed file appropriate format
  CpG_locs <-  cbind(data.frame(chromosomes), starts, ends, coefficients)
  names(CpG_locs) <- c("chromosome", "start", "end", "coefficient") 
  return(CpG_locs)
}


lift_over <- function(coordinates=NULL, chain_file=NULL) {
  ### Function lifts over genome coordinates from one assembly to another given a chain file
  
  if (chain_file=="-" | chain_file=="" | chain_file=="/" | chain_file=="\\" | is.null(chain_file) ) {
    #load the bed file containing genome coordinates 
    print("No lift over requested, proceeding without it.")
    return(coordinates)
  } else {
    library(rtracklayer)
    #import chain file
    chain <- import.chain(chain_file)
    #coerce dataframe into GRanges
    ranges <- makeGRangesFromDataFrame(coordinates)
    # use R implementation of Lift Over tool to convert ranges to hg38 coordinate system
    lifted_ranges <- liftOver(ranges, chain)
    # get the object in the right format and write to file
    temp_df <- data.frame(iranges=lifted_ranges)
    lifted_coor <- data.frame(cbind.data.frame(temp_df$iranges.seqnames,
                                               temp_df$iranges.start, temp_df$iranges.end))
    if("coefficient" %in% names(coordinates)){
      lifted_coor <- data.frame(cbind.data.frame(lifted_coor, coordinates$coefficient))
      names(lifted_coor) <- c("chromosome", "start", "end", "coefficient")
    } else {
      names(lifted_coor) <- c("chromosome", "start", "end")
    }
    
    return(lifted_coor)
  }
}


load_G4_data <- function(G4_file_plus=NULL, G4_file_minus=NULL, narrow_peak=F, autosomes=T){
  ### Function loads in given G4 data; when narrow_peak is true load from narrow_peak format
  library(dplyr)
  G4_locs <- read.csv(G4_file_plus, header=F, sep="\t")
  if (!is.null(G4_file_minus)){
    print("Two G4 files supplied.")
    G4_minus <- read.csv(G4_file_minus, header=F, sep="\t")
    G4_locs <- rbind.data.frame(G4_locs, G4_minus)
  }
  if (narrow_peak==T){
    keep = c(1, 2, 3, 5)
    G4_locs <- subset(G4_locs, select = keep)
  }
  if (dim(G4_locs)[2]<4){
    G4_locs <- G4_locs %>% mutate(score=NA)
  }
  names(G4_locs) <- c("chromosome", "start", "end", "score")
  G4_locs <- G4_locs %>% arrange(chromosome, start, end) 
  if(autosomes){
    G4_locs <- G4_locs %>% filter(chromosome!="chrX" & chromosome!="chrY" & chromosome!="chrM")
  }
  return(G4_locs)
}


load_ATAC_data <- function(ATAC_file=NULL, autosomes=T){
  ### Function loads given ATAC-seq data file, assumed to be .narrowPeak file
  
  library(dplyr)
  ATAC_data <- read.csv(ATAC_file, sep="\t", header=F)
  ATAC_data <- subset(ATAC_data, select=c(1, 2, 3, 5))
  names(ATAC_data) <- c("chromosome", "start", "end", "score")
  if(autosomes){
    ATAC_data <- ATAC_data %>% filter(chromosome!="chrX" & chromosome!="chrY")
  }
  return(ATAC_data)
}


annotate_CpGs <- function(CpG_coordinates=NULL, CGI_map_file=NULL, ATAC_map_file=NULL){
  ### Function to add columns to CpG_locs dataframe, such as CGI context, ATAC context
  # annotate CGI context
  if(is.null(CGI_map_file)){
    print("No CGI map file supplied. CGI context won't be annotated.")
    CpGs_annotated <- CpG_coordinates
  } else{
    library(dplyr)
    # Load the CGI map 
    CGI_coordinates <- read.csv(CGI_map_file, sep="\t", col.names = c("chromosome", "start", "end", "name", "?", "??"))
    CGI_coordinates <- CGI_coordinates[,c("chromosome", "start", "end")]
    # Overlap CpGs with CGI map to get CpGs in CGIs
    CpGs_in_CGIs <- overlap(query=CpG_coordinates, search_set=CGI_coordinates, reduce=F)
    CpGs_in_CGIs <- CpGs_in_CGIs[["overlap_subset"]]
    CpGs_in_CGIs <- CpG_locs %>% 
      filter(start %in% CpGs_in_CGIs$start & end %in% CpGs_in_CGIs$end) %>%
      arrange(chromosome, start) %>%
      mutate(CGI_context="inside")
    CpGs_outside_CGIs <- CpG_locs %>% 
      filter(!(start %in% CpGs_in_CGIs$start & end %in% CpGs_in_CGIs$end)) %>% 
      arrange(chromosome, start) %>% 
      mutate(CGI_context="outside")
    CpGs_annotated <- rbind(CpGs_in_CGIs, CpGs_outside_CGIs) %>% arrange(chromosome, start)
  }
  #annotate ATAc context
  if(is.null(ATAC_map_file)){
    print("No ATAC map file supplied. ATAC context won't be annotated.")
    CpGs_annotated <- CpGs_annotated
  } else{
    library(dplyr)
    ATAC_coordinates <- read.csv(ATAC_map_file, sep="\t", header=F, col.names=c("chromosome", "start", "end", "name", "strand", ".", "..", "...", "....", "....."))
    ATAC_coordinates <- ATAC_coordinates[,c("chromosome", "start","end")]
    CpGs_accessible <- overlap(query=CpGs_annotated, search_set=ATAC_coordinates, reduce=F)
    CpGs_accessible <- CpGs_accessible[["overlap_subset"]]
    CpGs_accessible <- CpGs_annotated %>% 
      filter(start %in% CpGs_accessible$start & end %in% CpGs_accessible$end) %>% 
      arrange(chromosome, start) %>% 
      mutate(chromatin="open")
    CpGs_inaccessible <- CpGs_annotated %>% 
      filter(!(start %in% CpGs_accessible$start & end %in% CpGs_accessible$end)) %>%
      arrange(chromosome, start) %>%
      mutate(chromatin="closed")
    CpGs_annotated <- rbind(CpGs_accessible, CpGs_inaccessible) %>% arrange(chromosome, start)
    
  }
  return(CpGs_annotated)
}


GRanges_to_df <- function(ranges=NULL){
  ### Function converts GRanges object to data frame
  df <- data.frame(iranges=ranges)
  df <- data.frame(cbind.data.frame(df$iranges.seqnames, df$iranges.start, df$iranges.end))
  names(df) <- c("chr", "start", "end")
  return(df)
}


chr_string <- function(number) {
  ## Function to convert numeric chromosome number to .bed file chromosome string (e.g.:"chr1")
  library(stringr)
  if (is.integer(number)){
    char <- as.character(number)
    chromosome = paste("chr", char, sep = "")
  }
  else if (str_detect(number, "chr")[1]) 
    chromosome = number
  else{
    chromosome=NaN
  }
  return(chromosome)
}


merge_coordinates <- function(df){
  ##Function to merge intervals in a data frame 
  library(GenomicRanges)
  ranges <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
  ranges <- GenomicRanges::reduce(ranges)
  df <- GRanges_to_df(ranges)
  return(df)
}


merge_replicates <- function(rep1, rep2){
  ## Function that merges 2 replicates and extracts only common peaks 
  library(GenomicRanges)
  merged <- overlap(query=rep1, search_set=rep2)
  merged <- merged[["overlapping_bases"]]
  return(merged)
  
}


p_value_MC <- function(num_overlaps, rand_overlaps){
  ## Function calculating the p-value given results from the MC-based method
  average_rand_overlaps <- mean(rand_overlaps)
  std_dev_rand_overlaps <- sd(rand_overlaps)
  z_score_overlap <- (num_overlaps - average_rand_overlaps) / std_dev_rand_overlaps
  p_value <- stats::pnorm(z_score_overlap, lower.tail = F)
  return(p_value)
}


sum_bases <- function(coordinates){
  bases <- sum(coordinates$end - coordinates$start + 1)
  return(bases)
}