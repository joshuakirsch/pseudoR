while getopts "s:d:r:c:f:m:h" flag
do
    case "${flag}" in
        s) sraList=${OPTARG};;
	   d) database=${OPTARG};;
	   r) ref=${OPTARG};;
		f) results=${OPTARG};;
        h)
		echo "This program aggregates IS hits in both contigs and ORFs."
		echo "	-s 	List of samples"
		echo "	-d	Location of helper scripts (usually repo location)"
		echo "	-m	m for one reference per sample and s for a single reference per sample"
        echo"   -r location of results folder (output folder from the main program)"
		echo "	-t	Number of threads"
;;

    esac
done

echo "Begin finding IS Circles"
Rscript --vanilla ${database}/IS-circles_analysis.R ${sraList}
echo "Finished finding IS Circles"
echo "Begin finding IS Insertions in Contigs"
Rscript --vanilla ${database}/post_analysis_V9_contigs.R ${sraList}
echo "Finished finding IS Insertions in Contigs"
echo "Begin finding IS Insertions in ORFs"
Rscript --vanilla ${database}/post_analysis_V9_orfs.R ${sraList}
echo "Finished finding IS Insertions in ORFs"

echo "Begin finding if IS Insertions in contigs are in ORFs"
cut -f1 final_results/contig_analysis.bed > temp/contig_ID.txt
grep -w -f temp/contig_ID.txt ref/contig.orf.bed > temp/orf.bed

seqkit fx2tab -n ref/orfs.nucl.fa > temp/orf_temp.txt
cut -f1 -d " " temp/orf_temp.txt > temp/temp_name_1.txt
cut -f3 -d " " temp/orf_temp.txt > temp/temp_name_2.txt
cut -f5 -d " " temp/orf_temp.txt > temp/temp_name_3.txt
cut -f7 -d " " temp/orf_temp.txt > temp/temp_name_4.txt
cut -f1 -d " " temp/orf_temp.txt | rev | cut -f2- -d "_" | rev > temp/temp_name_contig.txt

paste temp/temp_name_contig.txt temp/temp_name_2.txt temp/temp_name_3.txt temp/temp_name_1.txt> ref/contig2orf.txt

bedtools intersect -a ref/contig2orf.txt -b final_results/contig_analysis.bed > final_results/IS_hits_in_orfs.txt
bedtools intersect -wa -a temp/orf.bed -b final_results/contig_analysis.bed > temp/IS_ids.txt



echo "Finished finding if IS Insertions in contigs are in ORFs"

echo "Begin finding depth of originally mapped reads at points of IS insertions in ORFs"
for i in $(cat ${sraList})
do

samtools index -@ 3 mapping_files/${i}.orf.reads_mapped.bam mapping_files/${i}.orf.reads_mapped.bam.bai
mosdepth -t 3 -x -n --by final_results/orf_analysis.bed temp/temp mapping_files/${i}.orf.reads_mapped.bam
gunzip -c temp/temp.regions.bed.gz > mapping_files/${i}.orf.reads_mapped.depth


done
echo "Finished finding depth of originally mapped reads at points of IS insertions in ORFs"

echo "Begin finding depth of originally mapped reads at points of IS insertions in contigs"
for i in $(cat ${sraList})
do

samtools index -@ 3 mapping_files/${i}.contig.reads_mapped.bam mapping_files/${i}.contig.reads_mapped.bam.bai
mosdepth -t 3 -x -n --by final_results/contig_analysis.bed temp/temp mapping_files/${i}.contig.reads_mapped.bam
gunzip -c temp/temp.regions.bed.gz > mapping_files/${i}.contig.reads_mapped.depth

done
echo "Finished finding depth of originally mapped reads at points of IS insertions in contigs"
