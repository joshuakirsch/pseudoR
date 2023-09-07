while getopts "s:t:d:r:c:m:h" flag
do
    case "${flag}" in
        t) threads=${OPTARG};;
	d) database=${OPTARG};;
	r) results=${OPTARG};;
	h)
		echo "This program begins the process of finding IS insertions in metagenomes by first finding insertions in contigs."
		echo "	-d	Location of kraken2 database"
		echo "	-t	Number of threads"
		echo "  -r      Location of 'results' folder"

;;

    esac
done

#seqkit grep -f final_results/filtered_orfs_ID.txt ref/orfs.nucl.dedupe.fa > final_results/filtered_orfs.fna
kraken2 --db ${database} --threads ${threads} ${results}/final_results/filtered_orfs.fna --use-names > ${results}/final_results/orfs.kraken.out
cut -f2 ${results}/final_results/orfs.kraken.out > ${results}/temp/k1_temp1.txt
cut -f3 ${results}/final_results/orfs.kraken.out | cut -f1 -d " " > ${results}/temp/k1_temp2.txt
paste ${results}/temp/k1_temp1.txt ${results}/temp/k1_temp2.txt > ${results}/final_results/orfs_taxonomy.txt

#seqkit grep -f final_results/filtered_contigs_ID.txt ref/allContigs.fa > final_results/filtered_contigs.fna
kraken2 --db ${database} --threads ${threads} ${results}/final_results/filtered_contigs.fna --use-names > ${results}/final_results/contigs.kraken.out
cut -f2 ${results}/final_results/contigs.kraken.out > ${results}/temp/k1_temp1.txt
cut -f3 ${results}/final_results/contigs.kraken.out | cut -f1 -d " " > ${results}/temp/k1_temp2.txt
paste ${results}/temp/k1_temp1.txt ${results}/temp/k1_temp2.txt > ${results}/final_results/contigs_taxonomy.txt
