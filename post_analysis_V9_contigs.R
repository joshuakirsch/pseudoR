if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("readr" %in% rownames(installed.packages()) == FALSE) {install.packages("readr")}
if("tidyr" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyr")}
library (readr, quietly = T)
library (tidyr, quietly = T)
library(dplyr, quietly = T)

options(dplyr.summarise.inform = FALSE)

args = commandArgs(trailingOnly=TRUE)
#args=c("sampleList.txt", "blank", "results")
sampleList = read_tsv(args[1], col_names = "Sample",show_col_types = FALSE)
sam_file_cols = c("qname", "flag", "contig","pos", "mapq", "cigar")
depth_cols = c("contig", "start" ,"end", "mapped_depth", rep("blah", 3))
depth_cols_select = c("contig", "start" ,"end", "mapped_depth")
orf_loc_cols = c("orf", "start", "end")
analysis= data.frame()
IS_ends = data.frame()
contig_output = data.frame()

for (i in seq (nrow(sampleList))) {
  temp_IS = read_tsv(paste ("mapping_files/",sampleList$Sample[i], ".unmapped.contig_mapping.filtered.sam", sep=""), col_names = sam_file_cols,
                     show_col_types = FALSE) %>%
    mutate (sample = sampleList$Sample[i]) %>%
    select (c(sam_file_cols, sample))
    IS_annot = read_csv(paste ("mapping_files/",sampleList$Sample[i], "_pre_mapping_output.csv", sep=""), show_col_types = FALSE) %>%
      rename ("qname"=read_id ) %>% rename ('IS_end'=sseqid) %>%
      rowwise()%>%
      mutate (IS = strsplit(IS_end, "_")[[1]][1]) %>%
      #label which side the read was trimmed on
      mutate (end_trim = ifelse (between(qstart,0,15), "5prime",
                                    ifelse (between(qend, (read_length-15), read_length), "3prime","0")))
    temp_IS = temp_IS %>% left_join(IS_annot, by="qname") %>%
      #make sure that IS insertion is likely inside the ORF
      mutate (insertion_pos = ifelse ((flag==0 & end_trim=="5prime"), pos, 
                                      ifelse((flag==0 & end_trim=="3prime"), pos+trimmed_seq_length, 
                                             ifelse ((flag==16 & end_trim=="5prime"), pos+trimmed_seq_length,
                                                     ifelse ((flag==16 & end_trim=="3prime"), pos,NA)))))
    #sum all IS elements at each insertion position
    temp_IS_n = temp_IS %>% group_by (IS, contig, insertion_pos,sample) %>% summarise (IS_depth = n())
    #generate dataframe with the max depth for each orf
    insertion_df = temp_IS_n %>% ungroup() %>% filter (IS_depth >= 2)
    if (nrow(insertion_df)>0){
    print (paste ("Working on:  ", sampleList$Sample[i], sep=""))
    insertion_df$ins_id = seq(nrow(insertion_df))
    
    #build dataframe with the number each of IS insertions with the separate IS element ends (5' and 3')
    temp_IS_ends = temp_IS %>% group_by (IS, contig, insertion_pos,sample, IS_end) %>%
      summarise (IS_depth = n()) %>%
      separate(col=IS_end, into =c("IS_type", "IS_end"), sep ="_") %>%
      pivot_wider(names_from = IS_end, values_from=IS_depth )
    if(!("5prime" %in% colnames(temp_IS_ends))){
      temp_IS_ends$`5prime` = rep(NA, nrow(temp_IS_ends))
    }
    if(!("3prime" %in% colnames(temp_IS_ends))){
      temp_IS_ends$`3prime` = rep(NA, nrow(temp_IS_ends))
    }
    temp_IS_ends [is.na(temp_IS_ends)] = 0
    temp_IS_ends = temp_IS_ends %>% ungroup() %>% filter (`5prime`> 0 | `3prime` > 0)
    IS_ends = rbind (IS_ends, temp_IS_ends)
    
    
    five_ends = temp_IS_ends %>% filter (`5prime` > 0) %>% select (!`3prime`)

    three_ends  = temp_IS_ends %>% filter (`3prime` > 0) %>% select (!`5prime`)
    if (nrow(five_ends) > 0 & nrow(three_ends)>0){
      five_ends$itr = seq(nrow(five_ends))
    for (j in seq(nrow(five_ends))) {
      temp_ins = five_ends [j,]
      temp_five_total = five_ends %>% filter (contig == temp_ins$contig) %>%
        filter (IS==temp_ins$IS) %>% 
        filter (between(insertion_pos, temp_ins$insertion_pos-20, temp_ins$insertion_pos+20))%>%
        filter (! insertion_pos == temp_ins$insertion_pos)
      temp_three_total = three_ends %>% filter (contig == temp_ins$contig) %>%
        filter (IS==temp_ins$IS) %>% 
        filter (between(insertion_pos, temp_ins$insertion_pos-20, temp_ins$insertion_pos+20))
      if (nrow(temp_three_total) > 0 & nrow(temp_five_total) > 0) {
        max_itr = (rbind (temp_five_total, temp_ins) %>% filter (`5prime` == max (`5prime`)) %>% select (itr) %>% arrange (itr)) [1,]
          if (!(max_itr$itr %in% contig_output$max_itr)){
            
            grouped_hits = rbind (temp_five_total %>% rename ("depth" = `5prime`) %>% select (!itr), 
                                  temp_three_total%>% rename ("depth" = `3prime`), 
                                  temp_ins %>% rename ("depth" = `5prime`)%>% select (!itr))
            max_hit = (grouped_hits %>% filter (depth == max (depth)))[1,]
            temp_ins$max_depth_site_insertion_pos = max_hit$insertion_pos
            temp_ins$max_depth_site_depth = max_hit$depth
            temp_ins$`3prime` = sum (temp_three_total$`3prime`)
            temp_ins$`5prime` = (sum (temp_five_total$`5prime`) + temp_ins$`5prime`)
            temp_ins$max_itr = max_itr$itr
            contig_output = rbind (contig_output, temp_ins)
          }
      } else if (nrow(temp_three_total) > 0 & nrow(temp_five_total) == 0) {
          grouped_hits = rbind (temp_three_total%>% rename ("depth" = `3prime`), 
                                temp_ins %>% rename ("depth" = `5prime`)%>% select (!itr))
          max_hit = (grouped_hits %>% filter (depth == max (depth)))[1,]
          temp_ins$max_depth_site_insertion_pos = max_hit$insertion_pos
          temp_ins$max_depth_site_depth = max_hit$depth
          temp_ins$`3prime` = sum (temp_three_total$`3prime`)
          temp_ins$max_itr = temp_ins$itr
          contig_output = rbind (contig_output, temp_ins)
        }
    }
    }
    }
}    


contig_output = contig_output %>% relocate (sample, IS, contig, insertion_pos, `5prime`, `3prime`, max_depth_site_insertion_pos, max_depth_site_depth)
write_tsv (contig_output, file="final_results/contig_analysis.step1.tsv")
contig_bed_output = contig_output %>% select (contig, max_depth_site_insertion_pos) %>%
  mutate (start_pos = (max_depth_site_insertion_pos-1)) %>%
  mutate (max_depth_site_insertion_pos = ifelse(start_pos<0, 1, max_depth_site_insertion_pos)) %>%
  mutate (start_pos = ifelse(start_pos < 0, 0, start_pos)) %>%
  relocate (contig,start_pos,max_depth_site_insertion_pos)
write_tsv (contig_bed_output, file="final_results/contig_analysis.bed", col_names = FALSE)