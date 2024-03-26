while getopts "s:t:d:r:c:m:1:2:3:x:n:h" flag
do
    case "${flag}" in
     s) sraList=${OPTARG};;
     t) threads=${OPTARG};;
	d) database=${OPTARG};;
	c) dedupe=${OPTARG};;
	m) mode=${OPTARG};;
	r) reference=${OPTARG};;
 	1) read1_ending=${OPTARG};;
  	2) read2_ending=${OPTARG};;
   	3) singleton_ending=${OPTARG};;
    x) contig_ending=${OPTARG};;
	n) sraList_New=${OPTARG};;

    esac
done

cd results

max_threads=${threads}
wd=$(pwd)

rm ref/contigs*
rm -r blast_results
rm -r unmapped.fa.split

for i in $(cat ../${sraList_New})
do

#quality trim reads
bbduk.sh in1=../${dedupe}/${i}${read1_ending} \
in2=../${dedupe}/${i}${read2_ending} \
out=temp/read1.fq out2=temp/read2.fq \
qtrim=rl threads=${max_threads} overwrite=true

bbduk.sh in=../${dedupe}/${i}${singleton_ending} \
out=temp/readS.fq \
qtrim=rl threads=${max_threads} overwrite=true

#build contigs mapping database
seqkit_var="${i}""${contig_ending}"
seqkit seq ../${reference}/${seqkit_var} -m 1000 > ref/contigs.fa

bowtie2-build ref/contigs.fa ref/contigs --threads ${max_threads} -q

#map to contigs and do IS trimming/mapping
bowtie2 -x ref/contigs -p ${max_threads} -1 temp/read1.fq \
-2 temp/read2.fq \
-U temp/readS.fq -S temp/temp.sam
samtools view -Sb -@ {max_threads} temp/temp.sam > temp/temp1.bam
samtools sort -@ ${max_threads} -m 4G -o mapping_files/${i}.contig.reads_mapped.bam temp/temp1.bam

#samtools fastq -f 4 -N -@ ${max_threads} -s temp/unmapped.s.fq -1 temp/unmapped.1.fq -2 temp/unmapped.2.fq mapping_files/${i}.contig.reads_mapped.bam > blank.txt
#cat temp/unmapped.s.fq temp/unmapped.1.fq temp/unmapped.2.fq > mapping_files/${i}.unmapped.fq
#reformat.sh in=temp/temp.sam out=mapping_files/${i}.unmapped.fq unmappedonly=t overwrite=true primaryonly=t

samtools fastq -f 4 temp/temp.sam > mapping_files/${i}.unmapped.fq

seqkit replace -p .+ -r "seq_{nr}" mapping_files/${i}.unmapped.fq > mapping_files/${i}.unmapped.numbered.fq
seqkit fq2fa mapping_files/${i}.unmapped.numbered.fq > unmapped.fa

mkdir blast_results

seqkit split2 -p ${max_threads} --force unmapped.fa

ls unmapped.fa.split > unmapped_list.txt

cat unmapped_list.txt | xargs -P${max_threads} -n1 -I% blastn -db ${IR_database}/IRs.fa -query unmapped.fa.split/% -out blast_results/%_blast_results.txt -num_threads 1 -outfmt 6

cat blast_results/* > mapping_files/${i}_blast_results.txt


cut -f1 mapping_files/${i}_blast_results.txt | seqkit grep -f - mapping_files/${i}.unmapped.numbered.fq | seqkit fx2tab > blast_results_reads.tsv
timeout 20m Rscript --vanilla ${database}/V9_analysis.R ${i}
seqkit tab2fx ${i}_trimmed_reads.tsv > trimmed_reads.fq
bowtie2 -x ref/contigs --very-sensitive -p ${max_threads} -U trimmed_reads.fq | samtools view -@ 5 -Sb -o mapping_files/${i}.unmapped.contig_mapping.bam
samtools view -F 4 mapping_files/${i}.unmapped.contig_mapping.bam > mapping_files/${i}.unmapped.contig_mapping.filtered.sam

cut -f3 mapping_files/${i}.unmapped.contig_mapping.filtered.sam | sort | uniq > mapping_files/${i}.contig_IS_hits.txt
rm ref/contigs*
rm -r blast_results
rm -r unmapped.fa.split

#map to ORFs and do IS trimming/mapping [use previously found mapping reads and unmapped reads]
#samtools fastq -F 4 -N -@ ${max_threads} -s temp/mapped.s.fq -1 temp/mapped.1.fq -2 temp/mapped.2.fq mapping_files/${i}.contig.reads_mapped.bam > blank.txt
#repair.sh in=temp/mapped.1.fq in2=temp/mapped.2.fq out1=temp/mapped.repaired.1.fq out2=temp/mapped.repaired.2.fq overwrite=true
#reformat.sh in=temp/temp.sam out=temp/mapped.fastq mappedonly=t overwrite=true primaryonly=t
samtools fastq -F 4 temp/temp.sam > temp/mapped.fastq

repair.sh in=temp/mapped.fastq out1=temp/mapped.repaired.1.fq out2=temp/mapped.repaired.2.fq outs=temp/mapped.s.fq overwrite=true


bowtie2 -x ref/orfs.nucl -p ${max_threads} -1 temp/mapped.repaired.1.fq \
-2 temp/mapped.repaired.2.fq \
-U temp/mapped.s.fq -S temp/temp.sam
samtools view -Sb -@ {max_threads} temp/temp.sam > temp/temp1.bam
samtools sort -@ ${max_threads} -m 4G -o mapping_files/${i}.orf.reads_mapped.bam temp/temp1.bam

seqkit tab2fx ${i}_trimmed_reads.tsv > trimmed_reads.fq
bowtie2 -x ref/orfs.nucl --very-sensitive -p ${max_threads} -U trimmed_reads.fq | samtools view -@ 5 -Sb -o mapping_files/${i}.unmapped.orf_mapping.bam
samtools view -F 4 mapping_files/${i}.unmapped.orf_mapping.bam > mapping_files/${i}.unmapped.orf_mapping.filtered.sam

cut -f3 mapping_files/${i}.orf_mapping.filtered.sam | sort | uniq > mapping_files/${i}.contig_IS_hits.txt


done

bash ${database}/post-analysis_V3.sh -s ../${sraList} -d ${database} 
Rscript --vanilla ${database}/Final_Output.multiReference.R ../${sraList} final_results/contig_analysis.step1.tsv ${database} contig
Rscript --vanilla ${database}/Final_Output.multiReference.R ../${sraList} final_results/orf_analysis.step1.tsv ${database} ORF