while getopts "s:t:d:r:c:m:h" flag
do
    case "${flag}" in
        s) sraList=${OPTARG};;
        t) threads=${OPTARG};;
	d) database=${OPTARG};;
	c) dedupe=${OPTARG};;
	m) mode=${OPTARG};;
	r) reference=${OPTARG};;

    esac
done

mkdir results
cd results
mkdir final_results
mkdir mapping_files
mkdir temp
mkdir ref

max_threads=${threads}
wd=$(pwd)
echo $wd

seqkit seq ../${reference} -m 1000 > ref/contigs.fa
bowtie2-build ref/contigs.fa ref/contigs --threads ${max_threads} -q

python ${database}/pprodigal.py -i ref/contigs.fa -p meta -d ref/orfs.nucl.fa -T ${max_threads}
dedupe.sh in=ref/orfs.nucl.fa out=ref/orfs.nucl.dedupe.fa minidentity=90 overwrite=true threads=${max_threads}
seqkit fx2tab -n -i -l ref/orfs.nucl.dedupe.fa > ref/orfs.nucl.dedupe.lengths.txt
bowtie2-build ref/orfs.nucl.dedupe.fa ref/orfs.nucl --threads ${max_threads} -q

seqkit fx2tab -n ref/orfs.nucl.dedupe.fa > temp/orf_temp.txt
cut -f1 -d " " temp/orf_temp.txt > temp/temp_name_1.txt
cut -f3 -d " " temp/orf_temp.txt > temp/temp_name_2.txt
cut -f5 -d " " temp/orf_temp.txt > temp/temp_name_3.txt
cut -f7 -d " " temp/orf_temp.txt > temp/temp_name_4.txt
paste temp/temp_name_1.txt temp/temp_name_2.txt temp/temp_name_3.txt temp/temp_name_4.txt > ref/orfs.ID-loc.txt

seqkit fx2tab -n ref/orfs.nucl.fa > temp/orf_temp.txt
cut -f3 -d " " temp/orf_temp.txt > temp/temp_name_2.txt
cut -f5 -d " " temp/orf_temp.txt > temp/temp_name_3.txt
cut -f1 -d " " temp/orf_temp.txt | rev | cut -f2- -d "_" | rev > temp/temp_name_contig.txt
paste temp/temp_name_contig.txt temp/temp_name_2.txt temp/temp_name_3.txt > ref/contig.orf.bed


for i in $(cat ../${sraList})
do

#quality trim reads
bbduk.sh in1=../${dedupe}/${i}.dedupe_reads.1.fq.gz \
in2=../${dedupe}/${i}.dedupe_reads.2.fq.gz \
out=temp/read1.fq out2=temp/read2.fq \
qtrim=rl threads=${max_threads}

bbduk.sh in=../${dedupe}/${i}.dedupe_reads.s.fq.gz \
out=temp/readS.fq \
qtrim=rl threads=${max_threads}

#map to contigs and do IS trimming/mapping
bowtie2 -x ref/contigs -p ${max_threads} -1 temp/read1.fq \
-2 temp/read2.fq \
-U temp/readS.fq -S temp/temp.sam
samtools view -Sb -@ {max_threads} temp/temp.sam > temp/temp1.bam
samtools sort -@ ${max_threads} -m 4G -o mapping_files/${i}.contig.reads_mapped.bam temp/temp1.bam

samtools fastq -f 4 -N -@ ${max_threads} -s temp/unmapped.s.fq -1 temp/unmapped.1.fq -2 temp/unmapped.2.fq mapping_files/${i}.contig.reads_mapped.bam > blank.txt
cat temp/unmapped.s.fq temp/unmapped.1.fq temp/unmapped.2.fq > mapping_files/${i}.unmapped.fq

seqkit replace -p .+ -r "seq_{nr}" mapping_files/${i}.unmapped.fq > mapping_files/${i}.unmapped.numbered.fq
seqkit fq2fa mapping_files/${i}.unmapped.numbered.fq > unmapped.fa
blastn -db ${IR_database}/IRs.fa -query unmapped.fa -out mapping_files/${i}_blast_results.txt -num_threads ${max_threads} -outfmt 6
cut -f1 mapping_files/${i}_blast_results.txt | seqkit grep -f - mapping_files/${i}.unmapped.numbered.fq | seqkit fx2tab > blast_results_reads.tsv
timeout 20m Rscript --vanilla ${database}/V9_analysis.R ${i}
seqkit tab2fx ${i}_trimmed_reads.tsv > trimmed_reads.fq
bowtie2 -x ref/contigs --very-sensitive -p ${max_threads} -U trimmed_reads.fq | samtools view -@ 5 -Sb -o mapping_files/${i}.unmapped.contig_mapping.bam
samtools view -F 4 mapping_files/${i}.unmapped.contig_mapping.bam > mapping_files/${i}.unmapped.contig_mapping.filtered.sam

cut -f3 mapping_files/${i}.unmapped.contig_mapping.filtered.sam | sort | uniq > mapping_files/${i}.contig_IS_hits.txt

#map to ORFs and do IS trimming/mapping [use previously found mapping reads and unmapped reads]
samtools fastq -F 4 -N -@ ${max_threads} -s temp/mapped.s.fq -1 temp/mapped.1.fq -2 temp/mapped.2.fq mapping_files/${i}.contig.reads_mapped.bam > blank.txt
repair.sh in=temp/mapped.1.fq in2=temp/mapped.2.fq out1=temp/mapped.repaired.1.fq out2=temp/mapped.repaired.2.fq overwrite=true

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