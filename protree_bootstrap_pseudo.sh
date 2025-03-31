#!/bin/bash

#SBATCH --ntasks=5 #core(CPU)
#SBATCH --nodes=1 #Use 1 node
#SBATCH --mem=20G #Default memory per CPU is 3GB
#SBATCH --mail-user=domniki.manousi@nmbu.no ##Email me when job is done. read the first
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=tsrk
#SBATCH --constraint=avx2
#SBATCH --output=slurm_out/gig_amino_treesplits


out='/net/fs-2/scale/OrionStore/Projects/MSLab/Domniki/gig/2024_and_later/data/wasabi_trees/protein_trees/tree_split/iqtree_out'

#guidance trees
pg1=/mnt/project/MSLab/Domniki/gig/2024_and_later/data/wasabi_trees/prot_tree/guidance_coltrim_gig1.fa
pg2=/mnt/project/MSLab/Domniki/gig/2024_and_later/data/wasabi_trees/prot_tree/guidance_coltrim_gig2.fa

#construct a tree use 1000 bootstrap rouds and the best-fit model finder
singularity exec https://depot.galaxyproject.org/singularity/iqtree%3A2.3.3--h21ec9f0_0 iqtree -s ${pg1} -m MFP -redo -B 1000
singularity exec https://depot.galaxyproject.org/singularity/iqtree%3A2.3.3--h21ec9f0_0 iqtree -s ${pg2} -m MFP -redo -B 1000

#remove any potential pseudogenes from there
singularity exec https://depot.galaxyproject.org/singularity/treeshrink%3A1.3.9--py39r42hdfd78af_0 run_treeshrink -t ${pg1}.treefile -c -q 0.05 -o ${out}/pseudo_out -O g1s1_pseudo --force --k 22
singularity exec https://depot.galaxyproject.org/singularity/treeshrink%3A1.3.9--py39r42hdfd78af_0 run_treeshrink -t ${pg2}.treefile -c -q 0.05 -o ${out}/pseudo_out -O g1s2_pseudo --force --k 22
