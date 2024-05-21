SRR=SRR26527556
apptainer run docker://ncbi/sra-tools prefetch $SRR
apptainer run docker://ncbi/sra-tools fasterq-dump $SRR \
  --split-files --progress
rm -rf $SRR
mv ${SRR}_1.fastq R1.fq
mv ${SRR}_2.fastq R2.fq
pigz -p 12 R1.fq
pigz -p 12 R2.fq

python ../src/mutate_reads.py -r1 R1.fq.gz -f1 R2.fq.gz -r2 mutR2.fq.gz -f2 mutR2.fq.gz --num_reads_to_modify 1



