faafile=$1
outdir=$2
faprefix=$3
cd ~/ARG_combine/
cut -f1 ARG_combine_info_aa_use2.txt |sed -n '2,$p'| seqkit grep -f - ARG_combine_rmduplicated_addpoint_models.faa -o ARG_combine_rmduplicated.faa
. ~/miniconda3/bin/activate diamond
diamond blastp --very-sensitive -p 2 -d ~/ARG_combine/ARG_combine_rmduplicated_addpoint_models.faa.diamond -q $faafile -o ${outdir}ARG_combine_diamond/${faprefix}.diamond.out_addpoint_models -e 1e-3 -k 1 --max-hsps 1 --outfmt 6 qseqid qlen sseqid slen pident length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore
diamond blastp --very-sensitive -p 2 -d ~/ARG_combine/ARG_combine_rmduplicated.faa.diamond -q $faafile -o ${outdir}ARG_combine_diamond/${faprefix}.diamond.out -e 1e-3 -k 1 --max-hsps 1 --outfmt 6 qseqid qlen sseqid slen pident length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore
conda deactivate
perl ARG_analyse5_add_mutation_models_processResult.pl $faprefix 1
perl ARG_analyse5_add_mutation_models_processResult.pl $faprefix 2
