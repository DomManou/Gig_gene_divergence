#!/bin/bash

#sort or use sorted  broad/narrow peak files (output from the atac/chip-seq nf-core pipelines)
#run chip-r using the peaks (separate jobs by peak type and case)

module load Anaconda3
source activate idr

########
#h3k27ac
chipr -i polyic_k27ac_R1_peaks_sorted.narrowPeak polyic_k27ac_R2_peaks_sorted.narrowPeak polyic_k27ac_R3_peaks_sorted.narrowPeak polyic_k27ac_R4_peaks_sorted.narrowPeak polyic_k27ac_R5_peaks_sorted.narrowPeak polyic_k27ac_R6_peaks_sorted.narrowPeak \
-m 2 \
-o chipr_polyIC_k27ac.out

chipr -i pbs_k27ac_R1_peaks_sorted.narrowPeak pbs_k27ac_R2_peaks_sorted.narrowPeak pbs_k27ac_R3_peaks_sorted.narrowPeak pbs_k27ac_R4_peaks_sorted.narrowPeak pbs_k27ac_R5_peaks_sorted.narrowPeak pbs_k27ac_R6_peaks_sorted.narrowPeak \
-m 2 \
-o chipr_pbs_k27ac.out

########
#h3k27me
chipr -i polyic_k27me3_R1_peaks_sorted.broadPeak polyic_k27me3_R2_peaks_sorted.broadPeak polyic_k27me3_R3_peaks_sorted.broadPeak polyic_k27me3_R4_peaks_sorted.broadPeak polyic_k27me3_R5_peaks_sorted.broadPeak polyic_k27me3_R6_peaks_sorted.broadPeak \
-m 2 \
-o chipr_polyIC_k27me

chipr -i pbs_k27me3_R1_peaks_sorted.broadPeak pbs_k27me3_R2_peaks_sorted.broadPeak pbs_k27me3_R3_peaks_sorted.broadPeak pbs_k27me3_R4_peaks_sorted.broadPeak pbs_k27me3_R5_peaks_sorted.broadPeak pbs_k27me3_R6_peaks_sorted.broadPeak \
-m 2 \
-o chipr_pbs_k27me

#####
#atac
chipr -i polyic_atac_R1_peaks_sorted.narrowPeak polyic_atac_R2_peaks_sorted.narrowPeak polyic_atac_R3_peaks_sorted.narrowPeak polyic_atac_R4_peaks_sorted.narrowPeak polyic_atac_R5_peaks_sorted.narrowPeak polyic_atac_R6_peaks_sorted.narrowPeak \
-m 2 \
-o chipr_polyIC_atac.out

chipr -i pbs_atac_R1_peaks_sorted.narrowPeak pbs_atac_R2_peaks_sorted.narrowPeak pbs_atac_R3_peaks_sorted.narrowPeak pbs_atac_R4_peaks_sorted.narrowPeak pbs_atac_R5_peaks_sorted.narrowPeak pbs_atac_R6_peaks_sorted.narrowPeak \
-m 2 \
-o chipr_pbs_atac.out

conda deactivate

