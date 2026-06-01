#!/bin/bash
date
cd /lustre/home/niejingyi2023/project/bjCDC/RefSeq_Anno
. /lustre/home/niejingyi2023/miniconda3/bin/activate diamond
diamond makedb --in refseq_bacteria_assembly_summary_20240702_Salmonella_cds_rmduplicated.faa -d refseq_bacteria_assembly_summary_20240702_Salmonella_cds_rmduplicated.faa.diamond
conda deactivate
date
