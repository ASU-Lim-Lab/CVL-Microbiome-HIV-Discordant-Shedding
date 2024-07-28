#Align sequences with MAFFT.#
mafft --auto --nuc --thread 88 papilloma_E1E2L2L1.fasta > papilloma_E1E2L2L1_mafft-auto.fasta

#Trim alignment with trimAl.#
trimal -in papilloma_E1E2L2L1_mafft-auto.fasta -out papilloma_E1E2L2L1_mafft-auto_trimal.fasta -gappyout

#Run IQ-TREE 40 times in parallel on the trimmed alignment, with 25 bootstrap replicates each time.#
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot01
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot02
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot03
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot04
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot05
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot06
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot07
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot08
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot09
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot10
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot11
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot12
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot13
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot14
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot15
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot16
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot17
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot18
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot19
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot20
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot21
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot22
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot23
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot24
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot25
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot26
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot27
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot28
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot29
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot30
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot31
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot32
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot33
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot34
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot35
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot36
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot37
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot38
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot39
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -st DNA -m GTR+G -nt 2 -bo 25 -pre E1E2L2L1_stdBoot40

#Concatenate output trees and generate consensus tree.#
cat E1E2L2L1_stdBoot*.boottrees > E1E2L2L1_ABGMN_stdBoot_alltrees;
iqtree -con -t E1E2L2L1_ABGMN_stdBoot_alltrees;
iqtree -s papilloma_E1E2L2L1_mafft-auto_trimal.fasta -te E1E2L2L1_ABGMN_stdBoot_alltrees.contree -pre E1E2L2L1_ABGMN_stdBoot_alltrees.contree -m GTR+G

