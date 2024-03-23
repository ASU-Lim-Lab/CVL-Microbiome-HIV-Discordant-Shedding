#Usage: python parse_krakenuniqreports_mpa.py -i inputfile -o outputfile
import sys
import argparse
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True,
    dest='input_file')
parser.add_argument('-o', required=True,
    dest='output_file')
args=parser.parse_args()
input_file = open(args.input_file, 'r')
df = pd.read_table(input_file, delimiter="\t")
taxa_df = df[df['#Classification'].str.contains('s__')] #Select rows with species classification or lower.
drop_subs = taxa_df[~taxa_df['#Classification'].str.contains('s__.*x__')] #Drop rows with subspecies classification.
drop_taxa = drop_subs[~drop_subs['#Classification'].str.contains('x__Cyanobacteria/Melainabacteria_group|o__Rickettsiales|d__Eukaryota|d__Archaea|d__Viruses|x__plasmids')] #Drop unwanted taxa
no_zeros = drop_taxa.loc[(drop_taxa.sum(axis=1, numeric_only=True) != 0)]
no_zeros.to_csv(args.output_file, sep="\t", index=False) #Write output file.
input_file.close()