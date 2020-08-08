module load SMRTLink/7.0.1.66975
ccs $1 ccs.bam #$1 is your 3rd-sequencing data(usually end with .subreads.bam)

module load SAMtools/1.9
samtools view ccs.bam | awk '{OFS="\t"; print ">"$1"\n"$10}' - > ccs.fasta #change FORMAT from bam to fasta

module load cutadapt/1.9.1
cutadapt -g ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT -e 0 -o noadapter_ccs.fasta ccs.fasta #cut adapter

module load BLAST+/2.7.1
blastn -db /public/home/zmxiong/zmxiong/xiong/PPI-F1third/1_cell/pretreatment/zea-cdnaBase-blastdb/zea-cdnaBase -query noadapter_ccs.fasta -out noadapter_ccs.txt -evalue 1.0e-4 -num_thre
ads 4 -outfmt 7

module load BLAST/2.2.26-Linux_x86_64
python3 findadbd.py noadapter_ccs.txt noadapter_ccs.fasta PPTs.txt
#The sequence of primer, recombination site(ATTL) and their reverse complement sequence in FASTA format are necessary in current directory.
# AD.fa & AD_revcom.fa，BD.fa & BD_revcom.fa，attL.fa & attL_revcom.fa is correspond the sequence of ad-primer, bd-primer and ATTL site and their reverse complement sequence.
