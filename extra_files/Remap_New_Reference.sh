 while getopts "s:t:d:p:r:" flag
do
    case "${flag}" in
        s) sraList=${OPTARG};;
        t) threads=${OPTARG};;
        d) database=${OPTARG};;
        p) past_results=${OPTARG};;
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

cp ../${past_results}/mapping_files/*_pre_mapping_output.csv mapping_files

bowtie2-build ../${reference} ref/contigs --threads ${max_threads} -q
 
python ${database}/pprodigal.py -i ../${reference} -p meta -d ref/orfs.nucl.fa -T ${max_threads}
 
 
seqkit fx2tab -n ref/orfs.nucl.fa > temp/orf_temp.txt
cut -f3 -d " " temp/orf_temp.txt > temp/temp_name_2.txt
cut -f5 -d " " temp/orf_temp.txt > temp/temp_name_3.txt
cut -f1 -d " " temp/orf_temp.txt | rev | cut -f2- -d "_" | rev > temp/temp_name_contig.txt
paste temp/temp_name_contig.txt temp/temp_name_2.txt temp/temp_name_3.txt > ref/contig.orf.bed
 
 for i in $(cat ../${sraList})
do

samtools fastq -F 4 -N -@ ${max_threads} -s temp/mapped.s.fq -1 temp/mapped.1.fq -2 temp/mapped.2.fq ../${past_results}/mapping_files/${i}.contig.reads_mapped.bam > blank.txt
repair.sh in=temp/mapped.1.fq in2=temp/mapped.2.fq out1=temp/mapped.repaired.1.fq out2=temp/mapped.repaired.2.fq overwrite=true

bowtie2 -x ref/contigs -p ${max_threads} -1 temp/mapped.repaired.1.fq \
-2 temp/mapped.repaired.2.fq \
-U temp/mapped.s.fq -S temp/temp.sam
samtools view -Sb -@ {max_threads} temp/temp.sam > temp/temp1.bam
samtools sort -@ ${max_threads} -m 4G -o mapping_files/${i}.contig.reads_mapped.bam temp/temp1.bam

seqkit tab2fx ../${past_results}${i}_trimmed_reads.tsv > trimmed_reads.fq
bowtie2 -x ref/contigs --very-sensitive -p ${max_threads} -U trimmed_reads.fq | samtools view -@ 5 -Sb -o mapping_files/${i}.unmapped.contig_mapping.bam
samtools view -F 4 mapping_files/${i}.unmapped.contig_mapping.bam > mapping_files/${i}.unmapped.contig_mapping.filtered.sam

done

bash ${database}/post-analysis_V3.sh -s ../${sraList} -d ${database}
