###!/usr/bin/env 
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 14:35:31 2018
First Version was written on Tue Aug  9 08:40:55 2016

@author: alexanderkutschera
"""
# Script edited by me chatGPT 13/05/2025, changes:
## Converted from python2 to python3
## Removed soon to be deprecated Bio.Applications dependencies

import os
import csv
from Bio import SeqIO
import subprocess
from Bio.Blast import NCBIXML

# Set paths
data_path = "/home/callum/rotation2/genome_comparison/"
proteomes_path = os.path.join(data_path, "proteome_database")
query_path = os.path.join(data_path, "query_PAO1")
blastp_path = "/home/callum/blast/ncbi-blast-2.16.0+/bin/blastp"

# Gene lists
raetz_genes = ["", "LpxA", "LpxC", "LpxD", "LpxH", "LpxB", "LpxK", "WaaA", "PSPTO_0182", "HtrB"]
core_genes = ["", "MsbA", "WaaL", "WapR", "putative_glycosyl_transferase", "DnpA", "mig-14", "WapH", "putative_carbamoyltransferase", "putative_kinase", "WapQ", "WaaP", "WaaG", "WaaC", "WaaF"]
osa_genes = ["", "Wzz2", "Wzz", "WbpA", "WbpB", "WbpC", "WbpD", "WbpE", "Wzy", "Wzx", "HisH2", "HisF2", "WbpG", "WbpH", "WbpI", "WbpJ", "WbpK", "WbpL", "WbpM"]
cpa_genes = ["", "WbpZ", "WbpY", "WbpX", "Wzt", "Wzm", "WbpW", "Gmd", "Rmd", "PA5455", "PA5456", "PA5457", "PA5458", "PA5459"]
mod_genes = ["", "ArnT", "EptA", "PagL", "LpxO", "LpxF"]

# All genes
all_genes = [gene[:-6] for gene in os.listdir(query_path) if gene.endswith(".fasta")]

# Choose gene list
name = "biosynthesis"
genelist = cpa_genes

#use if you want to blast specific proteomes
#proteomes = ["Pci_JBC1_proteome", "Pst_DC3000_proteome"] 

# Get list of proteomes
proteomes = [filename[:-4] for filename in os.listdir(proteomes_path) if filename.endswith(".faa")]
print(f"You have {len(proteomes)} proteomes in your list")

# ---- blast ----

# Initialize dictionaries
d = {name: genelist}
n = {name: genelist}

# Run BLAST
for count, proteome in enumerate(proteomes, start=1):
    # Intialise lists
    percentage_list = []
    genename_list = []

    # Iterate through list of genes
    for gene in genelist[1:]:
        print(f"BLASTing {gene} on {proteome} ...")

        # Get query sequence
        qseq = next(SeqIO.parse(os.path.join(query_path, gene + ".fasta"), "fasta"))
        print(f"Query: {len(qseq.seq)}")

        # Now depends on subprocess
        query_file = os.path.join(query_path, gene) + ".fasta"
        db_path = os.path.join(proteomes_path, proteome)
        out_file = os.path.join(data_path, "result.xml")

        cmd = [
            blastp_path,
            "-query", query_file,
            "-db", db_path,
            "-evalue", "1",
            "-outfmt", "5",  # XML format
            "-out", out_file
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Error running BLAST: {result.stderr}")
        else:
            print("BLAST run completed.")

        # Check the results
        with open(os.path.join(data_path, "result.xml")) as result_handle:
            blast_record = NCBIXML.read(result_handle)
        # If there are no alignments, there are no hits
        if len(blast_record.alignments) == 0:
            percentage = 0
            genename = "no hit"
        # Otherwise this code just records the first hit by default
        else:
            length = blast_record.alignments[0].length
            title = blast_record.alignments[0].title
            # I don't understand what the PE bit is for
            if "PE=" not in title:
                genename = title[title.find("GN=")+3:title.find("SV=")]
            else:
                genename = title[title.find("GN=")+3:title.find("PE=")]

            hsp = blast_record.alignments[0].hsps[0]
            identity = hsp.identities
            print(f"identities: {hsp.identities}")
            print(f"length: {length}")
            print(f"alignlength: {hsp.align_length}")
            percentage = (float(identity) / float(length)) * 100
        # Get save gene and percentage similarity
        print(f"Identity percentage is {percentage:.2f} %\n")
        percentage_list.append(percentage)
        genename_list.append(genename)

    percentage_list.insert(0, proteome)
    genename_list.insert(0, proteome)
    d[proteome] = percentage_list
    n[proteome] = genename_list

# Save results
with open(os.path.join(data_path, "final_results.csv"), "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(d[name])
    for proteome in proteomes:
        writer.writerow(d[proteome])

with open(os.path.join(data_path, "final_results_genenames.csv"), "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(n[name])
    for proteome in proteomes:
        writer.writerow(n[proteome])

print("\ndone!")
