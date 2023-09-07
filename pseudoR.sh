while getopts "s:t:d:r:c:m:i:1:2:3:x:h" flag
do
    case "${flag}" in
        s) sraList=${OPTARG};;
        t) threads=${OPTARG};;
	d) database=${OPTARG};;
	c) dedupe=${OPTARG};;
	m) mode=${OPTARG};;
	r) reference=${OPTARG};;
    	i) IR_database=${OPTARG};;
	1) read1_ending=${OPTARG};;
  	2) read2_ending=${OPTARG};;
   	3) singleton_ending=${OPTARG};;
    	x) contig_ending=${OPTARG};;
	h)
		echo "This program begins the process of finding IS insertions in metagenomes by first finding insertions in contigs."
		echo "	-s 	List of samples"
		echo "	-d	Location of the repo"
		echo "	-c	Folder containing reads"
		echo "	-m	m for one reference per sample and s for a single reference per sample"
		echo "	-r	Folder of assemblies if m is selected or location of single reference contigs if s is selected"
		echo "	-t	Number of threads"
	 	echo "	-1	ending of read1 files. The sample name MUST proceed this delimiter. For instance, for sample1 the read1 file is sample1.read1.fq.gz. The read1 ending should be -1 'read1.fq.gz'"
   		echo "	-2	ending of read2 files. The sample name MUST proceed this delimiter. For instance, for sample1 the read2 file is sample1.read2.fq.gz. The read1 ending should be -2 'read2.fq.gz'"
     		echo "	-3	ending of single read files. The sample name MUST proceed this delimiter and you MUST provide a single read file. If none exists, make a blank file as placeholder"
       		echo "	-x	ending of contig files. The sample name MUST proceed this delimiter. For instance, for sample1 the assembly file is sample1.contigs.fa. The contig ending should be -x '.contigs.fa'
	 			Only used for multi-reference mode."
	 	echo "  -i	OPTIONAL: Folder location of IS element blast database (default is the same location as '-d'. However, a custom database can be used and will replace the provided IS element database from the ISOSDB. Please make this database using the IR-Database-Generation_V2.sh script)"
   		
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

source ${database}/pseudoR_multi-reference.sh -s ${sralist} -t ${threads} -d ${database} -c ${dedupe} -r ${reference} -i ${IR_database} -1 ${read1_ending} -2 ${read2_ending} -3 ${singleton_ending} -x ${contig_ending}

elif [ "${mode}" = "s" ]
then

source ${database}/pseudoR_single-reference.sh -s ${sralist} -t ${threads} -d ${database} -c ${dedupe} -r ${reference} -i ${IR_database} -1 ${read1_ending} -2 ${read2_ending} -3 ${singleton_ending}

fi

