#Align ORF1 sequences using MAFFT.#
mafft --auto --amino --thread 8 anello_ORF1.fasta > anello_ORF1_mafft-auto.fasta

#Trim alignment with trimAl.#
trimal -in anello_ORF1_mafft-auto.fasta -out anello_ORF1_mafft-auto_trimal.fasta -gappyout

#Run IQ-TREE 20 separate times in parallel on the trimmed alignment, with 50 bootstrap replicates each.#
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot01
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot02
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot03
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot04
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot05
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot06
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot07
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot08
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot09
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot10
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot11
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot12
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot13
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot14
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot15
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot16
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot17
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot18
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot19
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -st AA -m LG+F+G4 -nt 4 -bo 50 -pre stdBoot20

#Concatenate the output trees and generate consensus tree.#
cat stdBoot*.boottrees > anello_alltrees;
iqtree -con -t anello_alltrees;
iqtree -s anello_ORF1_mafft-auto_trimal.fasta -te anello_alltrees.contree -pre anello_alltrees.contree -m LG+F+G4
