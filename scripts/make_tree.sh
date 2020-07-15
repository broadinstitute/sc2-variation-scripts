#!/bin/bash

set -e -o pipefail

# rename all files to remove ".merged"
#for i in $(ls -1 ./*.fasta);do mv $i $(basename $i|perl -lape 's/(.*)(?:\.merged).fasta/$1.fasta/g'); done

# rename each sequence to remove ".merged"
#perl -i.bak -lape 's/\.merged//g' *.fasta

#combine all seqs to single fasta
#cat ./G*.fasta > all_seqs_unaligned.fasta

# perform MSA
#mafft --localpair --maxiterate 1000 --reorder --thread 8 --ep 0.123 all_seqs_unaligned.fasta > aligned.fasta

# trim with trimal, convert to phy format
trimal -phylip -automated1 -in aligned.fasta -out trimmed.phy

# run IQ-tree
iqtree -s trimmed.phy -bb 1000 -nt AUTO

# get IDs in tree
#cat aligned.fasta | grep '>' | perl -lape 's/(?:.*)(G\d+\.\d)(?:.*)/$1/g' | sort > ids_in_tree.txt
#cat aligned.fasta | grep '>' | perl -lape 's/(?:>)([^\s]*).*/$1/g' | sort > ids_in_tree.txt
