#!/bin/bash
inputfile=$1
outputfile=$2
. ~/miniconda3/bin/activate diamond
diamond blastp --very-sensitive -p $SLURM_NTASKS -d ~/RefSeq_Anno/refseq_bacteria_assembly_summary_20240702_Salmonella_cds_rmduplicated.faa.diamond -q $inputfile -o ${outputfile} -e 1e-3 -k 1 --max-hsps 1 --outfmt 6 qseqid qlen sseqid slen pident length qcovhsp scovhsp mismatch gapopen qstart qend sstart send qstrand evalue bitscore
conda deactivate
