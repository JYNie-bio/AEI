#!/bin/sh
date
cd /lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno
. /lustre/home/niejingyi2023/miniconda3/bin/activate seqkit
seqkit rmdup -s -D refseq_bacteria_assembly_summary_20240702_Salmonella_cds_faa_rmduplicated.info -i refseq_bacteria_assembly_summary_20240702_Salmonella_cds.faa -o refseq_bacteria_assembly_summary_20240702_Salmonella_cds_rmduplicated.faa
conda deactivate
date
