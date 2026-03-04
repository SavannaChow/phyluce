# -c 會在輸出中增加一致性資訊，--secondary=no 確保盡量減少多重比對的干擾

minimap2 -x sr -c --secondary=no target_contigs.fasta query_probes.fasta > output.paf

spades.py -t 32 -1 ./GPY-Ah4_R1_trimmed.fastq -2 ./GPY-Ah4_R2_trimmed.fastq -o ./GPY_Ah4


DG4_Ah_R1.fastq
DG4_Ah_R2.fastq


get the raw path
make a "01_Trimmed" folder for fastp output

get the sample ID R1.fastq and R2.fastq pair
ex: 
DG4_Ah_R1.fastq
DG4_Ah_R2.fastq
the sample ID is the DG4_Ah

command:
fastp 
--overrepresentation_analysis \
--detect_adapter_for_pe \
--thread #threads \
-i DG4_Ah_R1.fq 
-I DG4_Ah_R2.fq
-o DG4_Ah_trimmed_R1.fq 
-O DG4_Ah_trimmed_R2.fq

put trimmed files R1 & R2 in "sample ID" folder under trimmed
put the trimming report in sample ID sampleID_fastq_reports in the sub folder

then do the assembly

make all the assemble in a folder called 02_Spades_assembly
read fastq from trimmed folder
do:

spades.py 
-t 8
-careful
-1 trimmed_R1.fastq 
-2 trimmed_R2.fastq 
-o sampleID

偵測有多少threads, 然後跑多個spades.py，每個py跑8個，所以可以跑$threads/8個process

Spades will creat a folder from the sampleID, and put everything asselbly in the folder

跑完assembly後，做alignment with UCE，要建立03_Alignment folder

minimap2 -x sr -c --secondary=no target_contigs.fasta query_probes.fasta > sampleID.paf

target_contigs會在02_Spades_assembly/sampleID/contigs.fasta
query_probes要使用者輸入路徑，或是讀取目前目錄下是否有fasta讓人選
每個樣本的paf輸出記得放到子資料夾

04_prepare_uce_hits_from_paf

05_Filtering_uce_hits_trace
跑phyluce_assembly_filter_uce_hits_trace

06_extract_uce_fasta_from_hits




SampleID_pop_species.fastq
DG4_Ah_R1.fastq
DG4_Ah_R2.fastq

DG4_hyaD_Ah.fastq
GPY5_hyaC_Ah.fastq