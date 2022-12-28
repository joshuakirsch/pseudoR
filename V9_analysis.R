if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("readr" %in% rownames(installed.packages()) == FALSE) {install.packages("readr")}
if("tidyr" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyr")}
library (readr)
library (tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
blast_results=readr::read_tsv (paste("mapping_files/",args[1],"_blast_results.txt", sep=""), col_names = FALSE, show_col_types = FALSE)
colnames (blast_results) = c("read_id", 'sseqid', 'pident','lngth', 'mismatch', 
                             'gapopen','qstart', 'qend', 'sstart','send', 'evalue', 'bitscore')
reads = read_tsv("blast_results_reads.tsv", col_names = c("read_id","read_seq","read_qual"), show_col_types = FALSE)
blast_results = blast_results %>% left_join (reads) %>%
  mutate (read_length = nchar (read_seq)) %>%
  filter ((read_length-lngth) > 34)%>%
  filter (pident > 85) %>%
  filter (lngth > 27) %>%
  rowwise() %>%
  mutate (IR_end = strsplit(sseqid, "_")[[1]][2])%>%
  mutate (subj_max = max(c(sstart,send)))%>%
  mutate (subj_min = min(c(sstart,send))) %>%
  #make sure that the IR is positioned at the end of the read to avoid matches to the inside of the read
  filter (between (qstart,0,15) || between (qend, (read_length-15), read_length)) %>%
   #this block makes sure that the IR is positioned correctly so that the rest of the read is not the internal
  #parts of the IS element. For example, the 3' IRs should be at the start of the read (or end if RC)
  filter ( (IR_end=='3prime' && between (qstart,0,15) && sstart<send) ||
             (IR_end=='5prime' && between (qstart,0,15) && sstart>send) ||
             (IR_end=='5prime' && between (qend, (read_length-15), read_length) && sstart<send) ||
             (IR_end=='3prime' && between (qend, (read_length-15), read_length) && sstart>send)) %>%
  
  #make sure that the part of the IR that matches is correct end. For example, the 5' end of the 5prime IRs
  #is what should be included in the match. This also reduces matches to the inside of the IR.
  
  filter ((IR_end == '5prime' && between (subj_min,0,15))||  
            (IR_end == '3prime' && between (subj_max,135,150))) %>%
  #revise nchar to read_length
  mutate (trimmed_seq = ifelse (between(qstart,0,15), substring(read_seq, (qend+1), read_length), 
                                ifelse (between(qend, (read_length-15), read_length), substring (read_seq, 0, (qstart-1)),"0"))) %>%
  mutate (trimmed_qual = ifelse (between(qstart,0,15), substring(read_qual, (qend+1), read_length), 
                                ifelse (between(qend, (read_length-15), read_length), substring (read_qual, 0, (qstart-1)),"0"))) %>%
  mutate (trimmed_seq_length = nchar(trimmed_seq)) %>%
  group_by (read_id) %>% filter (trimmed_seq_length == max(trimmed_seq_length)) %>%
  filter (lngth == max(lngth)) %>% filter (pident == max(pident)) %>% ungroup () %>%
  distinct (read_id, .keep_all = TRUE)
write_csv(blast_results, paste("mapping_files/",args[1], "_pre_mapping_output.csv",sep=""))
write_tsv (blast_results %>% select (read_id, trimmed_seq, trimmed_qual), paste(args[1],"_trimmed_reads.tsv",sep=""), col_names = FALSE)

