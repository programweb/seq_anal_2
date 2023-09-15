from Bio import SeqIO
from Bio.SeqUtils import molecular_weight, gc_fraction
from Bio.Align import PairwiseAligner
from Bio import ExPASy
from Bio import SwissProt
import statistics

# Replace 'larger.fasta' with the path to your larger FASTA file
fasta_file = "GFOD2.fasta"

# Create empty lists to store various analysis results
sequence_ids = []
sequence_descriptions = []
sequence_lengths = []
sequence_means = []
sequence_medians = []
sequence_std_devs = []
at_contents = []
gc_contents = []
amino_acid_compositions = []
nucleotide_compositions = []
alignment_scores = []
functional_domain_predictions = []

# Read all sequences from the FASTA file
sequences = list(SeqIO.parse(fasta_file, "fasta"))

# Mapping of single-letter amino acid codes to full names
single_letter_to_full_name = {
    'A': 'Alanine',
    'C': 'Cysteine',
    'D': 'Aspartic acid',
    'E': 'Glutamic acid',
    'F': 'Phenylalanine',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'K': 'Lysine',
    'L': 'Leucine',
    'M': 'Methionine',
    'N': 'Asparagine',
    'P': 'Proline',
    'Q': 'Glutamine',
    'R': 'Arginine',
    'S': 'Serine',
    'T': 'Threonine',
    'V': 'Valine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine'
}

# Mapping of single-letter amino acid codes to three-letter abbreviations
single_letter_to_three_letter = {
    'A': 'Ala',
    'C': 'Cys',
    'D': 'Asp',
    'E': 'Glu',
    'F': 'Phe',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'K': 'Lys',
    'L': 'Leu',
    'M': 'Met',
    'N': 'Asn',
    'P': 'Pro',
    'Q': 'Gln',
    'R': 'Arg',
    'S': 'Ser',
    'T': 'Thr',
    'V': 'Val',
    'W': 'Trp',
    'Y': 'Tyr'
}


# Calculate sequence statistics
sequence_lengths = [len(record.seq) for record in sequences]
sequence_means = [statistics.mean(sequence_lengths)]
sequence_medians = [statistics.median(sequence_lengths)]
sequence_std_devs = [statistics.stdev(sequence_lengths)]

# Calculate molecular weight for the first sequence
mw = molecular_weight(sequences[0].seq)

# Calculate GC content for all sequences
gc_contents = [gc_fraction(record.seq) for record in sequences]

# Check if it's a protein sequence or a nucleotide sequence
for record in sequences:
    sequence_ids.append(record.id)
    sequence_descriptions.append(record.description)
    sequence_length = len(record.seq)
    at_content = (record.seq.count("A") + record.seq.count("T")) / sequence_length * 100
    at_contents.append(at_content)
    
    if all(base in "ACDEFGHIKLMNPQRSTVWY" for base in record.seq):
        # Amino acid composition (assuming it's a protein sequence)
        amino_acid_composition = {}
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            amino_acid_composition[aa] = record.seq.count(aa) / sequence_length * 100
        amino_acid_compositions.append(amino_acid_composition)
    else:
        amino_acid_compositions.append(None)

    # Nucleotide composition
    nucleotide_composition = {}
    for nt in "ACGT":
        nucleotide_composition[nt] = record.seq.count(nt) / sequence_length * 100
    nucleotide_compositions.append(nucleotide_composition)

# Create a PairwiseAligner object
aligner = PairwiseAligner()

# Calculate pairwise Smith-Waterman alignment scores for all sequences
for i in range(len(sequences)):
    for j in range(i+1, len(sequences)):
        alignment = aligner.align(sequences[i].seq, sequences[j].seq)
        alignment_score = alignment[0].score
        alignment_scores.append(alignment_score)

# Retrieve functional domain predictions using InterProScan
for record in sequences:
    sequence_id = record.id
    sequence = record.seq
    try:
        # Use ExPASy to get Swiss-Prot ID
        swissprot_id = ExPASy.get_sprot_raw(sequence_id)
        swissprot_record = SwissProt.read(swissprot_id)
        
        # Access InterProScan results for the Swiss-Prot record
        interproscan_results = swissprot_record.annotations.get('DR', [])
        functional_domains = [result for result in interproscan_results if 'InterPro' in result]
        functional_domain_predictions.append(functional_domains)
    except Exception as e:
        # Handle exceptions if the Swiss-Prot ID or InterProScan results are not available
        functional_domain_predictions.append('Prep sequence file and run InterProScan')

# Display the analysis results
for i in range(len(sequence_ids)):
    print("Sequence ID:", sequence_ids[i])
    print("Sequence Description:", sequence_descriptions[i])
    print("Sequence Length:", sequence_lengths[i])
    print("Sequence Mean Length:", sequence_means[0])
    print("Sequence Median Length:", sequence_medians[0])
    print("Sequence Standard Deviation:", sequence_std_devs[0])
    print("AT Content:", at_contents[i])
    print("GC Content:", gc_contents[i])
    print("Amino Acid Composition:", amino_acid_compositions[i])
    print("Nucleotide Composition:", nucleotide_compositions[i])
    print("Functional Domain Predictions:", functional_domain_predictions[i])

# Molecular Weight for the first sequence
print("Molecular Weight of First Sequence:", mw)

# Display pairwise alignment scores
for i in range(len(alignment_scores)):
    print(f"Alignment Score between Sequence {i//len(sequences)+1} and Sequence {i%len(sequences)+1}: {alignment_scores[i]}")

print("\n")
for k, v in single_letter_to_full_name.items():
    print('{:>20s} {:>10s} {:>14f}'.format(single_letter_to_full_name[k], single_letter_to_three_letter[k], amino_acid_compositions[0][k]))

print("\n")
