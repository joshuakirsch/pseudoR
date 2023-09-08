Clone directory using git: `git clone https://github.com/joshuakirsch/pseudoR.git`
Install dependencies with mamba [might work with conda but not sure]: 
`mamba create -n pseudoR bedtools=2.30.0 blast=2.12.0 bowtie2=2.4.5 mosdepth=0.3.3 prodigal=2.6.3 samtools=1.6 seqkit=2.2.0 r r-essentials r-readr r-tidyr r-dplyr`
BBTools must be installed using the instructions from this website: [BBTools] (https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/)
Add the BBTools folder to your $PATH variable before running the pseudoR pipeline
pseudoR pipeline can be run with a single command:
`bash pseudoR.sh`
Running `bash pseudoR.sh` shows the following information:
```
This program begins the process of finding IS insertions in metagenomes by first finding insertions in contigs.
        -s      List of samples
        -d      Location of the repo
        -c      Folder containing reads
        -m      m for one reference per sample and s for a single reference per sample
        -r      Folder of assemblies if m is selected or location of single reference contigs if s is selected
        -t      Number of threads
        -1      ending of read1 files. The sample name MUST proceed this delimiter. For instance, for sample1 the read1 file is sample1.read1.fq.gz. The read1 ending should be -1 'read1.fq.gz'
        -2      ending of read2 files. The sample name MUST proceed this delimiter. For instance, for sample1 the read2 file is sample1.read2.fq.gz. The read1 ending should be -2 'read2.fq.gz'
        -3      ending of single read files. The sample name MUST proceed this delimiter and you MUST provide a single read file. If none exists, make a blank file as placeholder
        -x      ending of contig files. The sample name MUST proceed this delimiter. For instance, for sample1 the assembly file is sample1.contigs.fa. The contig ending should be -x '.contigs.fa' Only used for multi-reference mode.
  	-i    OPTIONAL: Folder location of IS element blast database (default is the same location as '-d'. However, a custom data         base can be used and will replace the provided IS element database from the ISOSDB. Please make this database using the IR-Database-Generation_V2.sh script)


```
I recommend using over 10 cores. There is a lot of mapping and BAM sorting in this program which takes a good bit of time. Usually, I budget 1 hour per sample. I also highly recommend that reads are deduplicated before using in this program.

This program does not create a log file automatically. To obtain a record of the program's readout, run `bash pseudoR.sh 1>>log.file 2>>log.file`

Files in the `results/ref/` folder can be used for downstream analyses if need be:
* orfs.nucl.fa – All predicted ORFs found in contigs greater than 1 kB
* orfs.nucl.dedupe.fa – Deduplicated ORFs from orfs.nucl.fa

There are many output files from this program and they can take up a lot of space. For each sample (in this case sample1), the following files in `results/mapping_files/` are produced:
* `sample1.contig.reads_mapped.bam` – sorted BAM file of reads mapped to whole contigs
* `sample1.unmapped.fq` – FASTQ file of reads that did not map to the contigs
* `sample1.unmapped.numbered.fq`  - FASTQ file with the same reads as in sample1.unmapped.fastq but with the names changed to seq_#, where # is a number
* `sample1_blast_results.txt` – Output of BLASTN from comparing sample1.unmapped.numbered.fq  to the IS termini database
* `sample1_trimmed_reads.tsv` – Unmapped reads trimmed of the IS termini [is found in results/]
* `sample1.unmapped.contig_mapping.bam` – Unsorted BAM file of unmapped and trimmed reads mapped to contigs
* `sample1.unmapped.contig_mapping.filtered.sam` – SAM version of sample1.unmapped.contig_mapping.bam without reads that did not map
* `sample1.orf.reads_mapped.bam` – Sorted BAM file of reads that mapped to the contigs and also mapped to the ORF database
* `sample1.unmapped.orf_mapping.bam` – Unsorted BAM file of unmapped and trimmed reads mapped to contigs
* `sample1.unmapped.orf_mapping.filtered.sam` – SAM version of sample1.unmapped.orf_mapping.bam without reads that did not map
* `sample1.contig.reads_mapped.depth` – Read depth of every position where an IS insertion was predicted in contigs
* `sample1.orf.reads_mapped.depth` – Read depth of every position where an IS insertion was predicted in the ORF database

Useful Output files are 
1. `final_results/pseudoR_output.contig.tsv`
	- This the contig mapping output of the pipeline. Here are what the columns in the tsv are:
		- `sample` = Sample where the insertion was found
		- `contig` = Contig where the insertion was found
		- `Insertion_Position` = Site of highest IS read depth from insertion
		- `IS_element` = IS element forming insertion
		- `5prime_IS_depth` = Depth of reads from 5' terminus of IS element
		- `3prime_IS_depth` = Depth of read from 3' terminus of IS element
		- `total_IS_depth` = Sum of 5prime_IS_depth and 3prime_IS_depth
		- `non_IS_depth` = Non-IS read depth (ie reads that originally mapped to this loci) at the insertion position
		- `depth_percentage`  = Percentage of total_IS_depth/non_IS_depth (ie 100* (total_IS_depth/non_IS_depth))
		- `ORF` = ORF where insertion is found. NA if insertion is intergenic
2. `final_results/pseudoR_output.ORF.tsv`
	- This the contig mapping output of the pipeline. Here are what the columns in the tsv are:
		- `sample` = Sample where the insertion was found
		- `ORF` = ORF where the insertion was found
		- `Insertion_Position` = Site of highest IS read depth from insertion
		- `IS_element` = IS element forming insertion
		- `5prime_IS_depth` = Depth of reads from 5' terminus of IS element
		- `3prime_IS_depth` = Depth of read from 3' terminus of IS element
		- `total_IS_depth` = Sum of 5prime_IS_depth and 3prime_IS_depth
		- `non_IS_depth` = Non-IS read depth (ie reads that originally mapped to this loci) at the insertion position
		- `depth_percentage`  = Percentage of total_IS_depth/non_IS_depth (ie 100* (total_IS_depth/non_IS_depth))
