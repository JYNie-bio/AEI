faafile=$1
outputfile=$2
outputfile2=$2

########### diamond_ARG_combine #######################
#mkdir ${outdir}
. ~/miniconda3/bin/activate diamond
diamond blastp --very-sensitive -p $SLURM_NTASKS -d ~/ARG_combine/ARG_combine_rmduplicated.faa.diamond -q $faafile -o ${outputfile} -e 1e-3 -k 1 --max-hsps 1 --outfmt 6 qseqid qlen sseqid slen pident length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore
diamond blastp --very-sensitive -p $SLURM_NTASKS -d ~/ARG_combine/ARG_combine_rmduplicated_addpoint_models.faa.diamond -q $faafile -o ${outputfile2} -e 1e-3 -k 1 --max-hsps 1 --outfmt 6 qseqid qlen sseqid slen pident length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore
conda deactivate
