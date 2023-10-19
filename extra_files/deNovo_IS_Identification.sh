cat assemblies/*.fa > assemblies/allContigs.fa
isescan.py --seqfile assemblies/allContigs.fa --output ISEScan_Output --nthread 20

#filter ISEScan hits to only include complete hits over 500bp
awk 'BEGIN {FS=","; OFS="_"} {if ($22=="c" && $6>500) print $1,$4,$5,$18}' ISEScan_Output/assemblies/allContigs.fa.csv > ISEScan_Output/filtered_IS_elements.ID.txt 
#remove weird description section of IS element headings in fna file from ISEscan
seqkit replace -p "\s.+" ISEScan_Output/assemblies/allContigs.fa.is.fna > ISEScan_Output/assemblies/allContigs.fa.is.noDescriptions.fna
#grab predicted IS elements
seqkit grep -f ISEScan_Output/filtered_IS_elements.ID.txt ISEScan_Output/assemblies/allContigs.fa.is.noDescriptions.fna > ISEScan_Output/filtered_IS_elements.fna

#get IRs from IS elements (this block has been checked)
mkdir temp
mkdir IR_database
cd-hit-est -i ISEScan_Output/filtered_IS_elements.fna -o temp/IS_nucl.fa -T 4 -c 0.9
perl clstr2txt.pl temp/IS_nucl.fa.cdhit.clstr > IR_database/IS.clstr.txt
seqkit fx2tab temp/IS_nucl.fa > temp/IS_nucl.txt
cut -f1 temp/IS_nucl.txt | cut -f1 -d "_" | awk '{print $0 "_5prime"}' > temp/name_5.txt
cut -f1 temp/IS_nucl.txt | cut -f1 -d "_" | awk '{print $0 "_3prime"}' > temp/name_3.txt
cut -f2 temp/IS_nucl.txt > temp/seqs.txt
cut -c -150 temp/seqs.txt > temp/seqs_5_150bp.txt
paste temp/name_5.txt temp/seqs_5_150bp.txt > temp/IRs.5.txt
rev temp/seqs.txt | cut -c -150 | rev > temp/seqs_3_150bp.txt
paste temp/name_3.txt temp/seqs_3_150bp.txt > temp/IRs.3.txt
cat temp/IRs.5.txt temp/IRs.3.txt > temp/IRs.txt
seqkit tab2fx temp/IRs.txt > IR_database/IRs.fa
makeblastdb -in IR_database/IRs.fa -dbtype 'nucl'
