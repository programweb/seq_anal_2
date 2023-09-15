# seq_anal_2
Biological Sequence Analysis Script

This Python script is designed for comprehensive biological sequence analysis and is particularly useful for researchers and bioinformaticians working with large FASTA files. Powered by the Biopython library, it automates a range of tasks related to biological sequences.

The script begins by reading sequences from a user-specified FASTA file and then computes essential sequence statistics, including length, mean length, median length, and standard deviation. It identifies whether each sequence is a protein or nucleotide sequence, allowing it to compute amino acid or nucleotide composition percentages accordingly.

Another standout feature is its ability to calculate pairwise sequence alignment scores. Using the Smith-Waterman algorithm, it determines the similarity or dissimilarity between sequences, offering insights into evolutionary relationships or structural similarities.

Additionally, the script incorporates functional domain prediction using InterProScan, enhancing its utility for researchers interested in understanding the potential functions of protein sequences.

The GFOD2.fasta file is used but can be replaced with any other FASTA file.  Scripts like this, with their versatility and automation capabilities, save valuable time and effort in biological sequence analysis, making it a valuable tool for bioinformaticians and biologists alike.
