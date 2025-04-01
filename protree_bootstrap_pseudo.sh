#!/bin/bash

#SBATCH --ntasks=5 #core(CPU)
#SBATCH --nodes=1 #Use 1 node
#SBATCH --mem=20G #Default memory per CPU is 3GB
#SBATCH --mail-user=domniki.manousi@nmbu.no ##Email me when job is done. read the first
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_out/gig_amino_tree

#Use the produced alignments after amino-acid residue reiability filtering (GUIDANCE2): format guidance output into fasta file
fold -w60 protein_trees/guidance_gig1.fa.txt > files/guidance_coltrim_gig1.fa
fold -w60 protein_trees/guidance_gig2.fa.txt > files/guidance_coltrim_gig2.fa

#variables
pg1=files/guidance_coltrim_gig1.fa
pg2=files/guidance_coltrim_gig2.fa
out='protein_trees/tree_split/iqtree_out'

#construct a tree use 1000 bootstrap rouds and the best-fit model finder option
singularity exec https://depot.galaxyproject.org/singularity/iqtree%3A2.3.3--h21ec9f0_0 iqtree -s ${pg1} -m MFP -redo -B 1000
singularity exec https://depot.galaxyproject.org/singularity/iqtree%3A2.3.3--h21ec9f0_0 iqtree -s ${pg2} -m MFP -redo -B 1000

#remove any potential pseudogenes from the constructed trees
singularity exec https://depot.galaxyproject.org/singularity/treeshrink%3A1.3.9--py39r42hdfd78af_0 run_treeshrink -t ${pg1}.treefile -c -q 0.05 -o ${out}/pseudo_out -O g1s1_pseudo --force --k 22
singularity exec https://depot.galaxyproject.org/singularity/treeshrink%3A1.3.9--py39r42hdfd78af_0 run_treeshrink -t ${pg2}.treefile -c -q 0.05 -o ${out}/pseudo_out -O g1s2_pseudo --force --k 22
