#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --job-name=fimo
#SBATCH --mem=80G
#SBATCH --array=1-1957%10
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_out/fimo_motif_%a

## no parallelization option: read through single style motifs
# JASPAR motifs were downloaded from https://jaspar.elixir.no/downloads/

SEQUENCE="references/Salmo_salar-GCA_905237065.2-unmasked.fa"
MFILE=$(ls fimo_motif/jaspar_individ_motif/*.meme | awk 'NR=='$SLURM_ARRAY_TASK_ID )
MOTIF=$(basename $MFILE | sed 's/.meme//')
OUTDIR=fimo_motif/singl_mtf/$MOTIF

#single motif style
singularity exec https://depot.galaxyproject.org/singularity/meme%3A5.5.2--py39pl5321h290edd5_0 fimo \
  -o $OUTDIR \
  --max-stored-scores 1000000 \
  $MFILE $SEQUENCE

# make table of all motif results and filter by Q-value (FDR) to reduce the size of the dataset
cat fimo_motif/singl_mtf/MA*/fimo.tsv | grep -v "#" | awk '{ if ($8 < 0.05) { print } }' > qval_mtf_05.tsv
