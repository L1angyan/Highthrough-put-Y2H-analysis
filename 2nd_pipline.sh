module load BLAST+/2.9.0
blastn -db /public/home/zmxiong/cDNA-ref/B73-cDNA -query F-1_R1.fa -out F-1_R1.txt -evalue 1.0e-4 -num_threads 10 -outfmt 7
blastn -db /public/home/zmxiong/cDNA-ref/B73-cDNA -query F-1_R2.fa -out F-1_R2.txt -evalue 1.0e-4 -num_threads 10 -outfmt 7
python3 y2h_2nd.py F-1_R1.txt F-1_R2.txt F-1_result.txt
