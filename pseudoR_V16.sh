while getopts "s:t:d:r:c:m:i:h" flag
do
    case "${flag}" in
        s) sraList=${OPTARG};;
        t) threads=${OPTARG};;
	d) database=${OPTARG};;
	c) dedupe=${OPTARG};;
	m) mode=${OPTARG};;
	r) reference=${OPTARG};;
    i) IR_database=${OPTARG};;
	h)
		echo "This program begins the process of finding IS insertions in metagenomes by first finding insertions in contigs."
		echo "	-s 	List of samples"
		echo "	-d	Location of the repo"
		echo "	-c	Folder containing reads"
		echo "	-m	m for one reference per sample and s for a single reference per sample"
		echo "	-r	Folder of assemblies if m is selected or location of single reference contigs if sr is selected"
		echo "	-t	Number of threads"
        echo "  -i  OPTIONAL: Folder location of IS element blast database (default is the same location as '-d'. However, a custom database can be used and will replace the provided IS element database from ISFinder. Please make this database using the IR-Database-Generation_V2.sh script)"
;;

    esac
done

if [ -z "$IR_database" ]
then
      IR_database=${database}
fi
echo "Location of IS Database:  ${IR_database}"

if [ "$mode" = "m" ]
then

source ${database}/pseudoR_multi-reference.sh -s ${sralist} -t ${threads} -d ${database} -c ${dedupe} -r ${reference} -i ${IR_database}

elif [ "${mode}" = "s" ]
then

source ${database}/pseudoR_single-reference.sh -s ${sralist} -t ${threads} -d ${database} -c ${dedupe} -r ${reference} -i ${IR_database}

fi

#cat mapping_files/*.contig_IS_hits.txt | sort| uniq > final_results/IS_insertion_orfs.ID.txt