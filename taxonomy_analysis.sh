while getopts "s:t:d:r:c:m:h" flag
do
    case "${flag}" in
        t) threads=${OPTARG};;
	d) database=${OPTARG};;
	h)
		echo "This program begins the process of finding IS insertions in metagenomes by first finding insertions in contigs."
		echo "	-d	Location of kraken2 database"
		echo "	-t	Number of threads"
;;

    esac
done

seqkit grep -f final_results/filtered_orfs_ID.txt ref/orfs.nucl.dedupe.fa > final_results/filtered_orfs.fna
kraken2 --db ../${database} --threads ${threads} final_results/filtered_orfs.fna --use-names > final_results/orfs.kraken.out
cut -f2 final_results/orfs.kraken.out > temp/k1_temp1.txt
cut -f3 final_results/orfs.kraken.out | cut -f1 -d " " > temp/k1_temp2.txt
paste temp/k1_temp1.txt temp/k1_temp2.txt > final_results/orfs_taxonomy.txt

seqkit grep -f final_results/filtered_contigs_ID.txt ref/allContigs.fa > final_results/filtered_contigs.fna
kraken2 --db ../${database} --threads ${threads} filtered_contigs.fna --use-names > contigs.kraken.out
cut -f2 contigs.kraken.out > temp/k1_temp1.txt
cut -f3 contigs.kraken.out | cut -f1 -d " " > temp/k1_temp2.txt
paste temp/k1_temp1.txt temp/k1_temp2.txt > final_results/contigs_taxonomy.txt