if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("readr" %in% rownames(installed.packages()) == FALSE) {install.packages("readr")}
if("tidyr" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyr")}
library (readr, quietly = T)
library (tidyr, quietly = T)
library(dplyr, quietly = T)

options(dplyr.summarise.inform = FALSE)

args = commandArgs(trailingOnly=TRUE)

sampleList = read_tsv(args[1], col_names = "Sample",
                      show_col_types = FALSE)
  analysis = read_tsv(args[2], 
                      show_col_types = FALSE) %>%
    rowwise() %>%
    mutate (total_IS_depth = sum(`5prime`+`3prime`)) %>%
    filter (total_IS_depth >= 4) %>%
    mutate (low_end = min (c(`5prime`,`3prime`))) %>%
    mutate (high_end = max (c(`5prime`,`3prime`))) %>%
    filter (low_end >= 0.1*(high_end)) %>%
    distinct()
  
  
  nonIS_mapping  = data.frame()
  for (i in seq(nrow(sampleList))){
    temp = read_tsv (paste("mapping_files/", sampleList$Sample[i], ".contig.reads_mapped.depth", sep=""),
                     col_names = c("contig","start_pos","max_depth_site_insertion_pos","non_IS_depth"), 
                     show_col_types = FALSE) %>%
      mutate (sample = sampleList$Sample[i]) %>%
      mutate (max_depth_site_insertion_pos = ifelse(start_pos==0, 0, max_depth_site_insertion_pos))
    nonIS_mapping = rbind (nonIS_mapping, temp)
  }
  nonIS_mapping = nonIS_mapping %>% distinct()
  
  analysis_annot = analysis %>% rowwise () %>%
    left_join (nonIS_mapping, by=c("contig" = "contig", "max_depth_site_insertion_pos" = "max_depth_site_insertion_pos",
                                   "sample"= "sample")) %>%
    mutate (depth_percentage = (100*(sum(`5prime`+`3prime`)/(sum(`5prime`+`3prime`)+non_IS_depth))))
  
    ins_itr = analysis_annot %>% select (contig, insertion_pos) %>% distinct () %>%
      ungroup() %>%
      arrange (contig, insertion_pos) %>%
      mutate (lag_insertion_pos = lag(insertion_pos)) %>%
      mutate (adj_insertion_pos = insertion_pos) %>%
      mutate (lag_is_nearby =  (contig == lag(contig) &
                                  between (insertion_pos, (lag_insertion_pos-15), (lag_insertion_pos+15)) &
                                  !is.na (lag_insertion_pos) &
                                  lag(adj_insertion_pos) != adj_insertion_pos))
    
    while(any(ins_itr$lag_is_nearby)) {
      ins_itr = ins_itr %>% mutate (adj_insertion_pos = ifelse (lag_is_nearby,
                                                                lag(adj_insertion_pos), adj_insertion_pos)) %>%
        mutate (lag_is_nearby =  (contig == lag(contig) &
                                    between (insertion_pos, (lag_insertion_pos-15), (lag_insertion_pos+15)) &
                                    !is.na (lag_insertion_pos) &
                                    lag(adj_insertion_pos) != adj_insertion_pos))
    }
    temp = ins_itr %>% select (contig, adj_insertion_pos) %>% distinct ()
    temp$unique_ins_itr = seq(nrow(temp))
    ins_itr = ins_itr %>% left_join(temp)
    
    analysis_annot = analysis_annot %>% left_join(ins_itr, by=c("contig" = "contig", "insertion_pos"="insertion_pos")) %>%
      select (sample, IS,contig, adj_insertion_pos,`5prime`, `3prime`,
              total_IS_depth, "max_depth_site_insertion_pos", max_depth_site_depth, non_IS_depth, depth_percentage) %>%
      rename ("Insertion_Position" = adj_insertion_pos) %>%
      rename ("5prime_IS_depth" = `5prime`) %>%
      rename ("3prime_IS_depth" = `3prime`) 

#this block removes IS elements false hits (likely) and adds a better naming strategy for the IS elements

IS_fam_df  = read_tsv(paste (args[3],"/IS_fam_annot.txt", sep="")) %>%
  rowwise () %>%
  mutate (IS = gsub("_", ":", IS))

analysis_annot = allResults %>% 
  left_join(IS_fam_df)

if (args[4] == "contig"){
  inOrf  = read_tsv ("final_results/IS_hits_in_orfs.txt", 
                     col_names = c("contig","start_pos","max_depth_site_insertion_pos", "ORF"),
                     show_col_types = FALSE) %>%
    select (contig,max_depth_site_insertion_pos,ORF) %>%
    mutate (intragenic = "Yes") %>% distinct(contig, max_depth_site_insertion_pos, .keep_all = TRUE)
  analysis_annot = analysis_annot %>%
    left_join(inOrf, by=c("contig" = "contig", "max_depth_site_insertion_pos" = "max_depth_site_insertion_pos")) %>% 
    mutate (intragenic =ifelse(is.na(intragenic), "No", intragenic)) %>%
    distinct() %>%
    rename ("Max_Depth_Site_Position" = max_depth_site_insertion_pos) %>%
    rename ("IS_Depth_at_Max_Depth_Site" = max_depth_site_depth)
  write_tsv (analysis_annot, "final_results/pseudoR_output.contig.tsv")
} else if (args[4] == "ORF"){
  analysis_annot = analysis_annot %>%
    mutate ("ORF"=contig) %>%
    rename ("Max_Depth_Site_Position" = max_depth_site_insertion_pos) %>%
    rename ("IS_Depth_at_Max_Depth_Site" = max_depth_site_depth)
  write_tsv (analysis_annot, "final_results/pseudoR_output.ORF.tsv")
}
