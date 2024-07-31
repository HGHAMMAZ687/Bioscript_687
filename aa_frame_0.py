import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

# Load sequences
sequences = list(SeqIO.parse("/home/uatrs-24/aligned.fasta", "fasta"))
reference = SeqIO.read("/home/uatrs-24/Documents/Covid_2026/WG.fasta", "fasta")

# Define gene annotations
gene_annotations = {
    "nsp1": (266, 805),
    "nsp2": (806, 2719),
    "nsp3": (2720, 8554),
    "nsp4": (8555, 10054),
    "nsp5": (10055, 10972),
    "nsp6": (10973, 11842),
    "nsp7": (11843, 12091),
    "nsp8": (12092, 12685),
    "nsp9": (12686, 13024),
    "nsp10": (13025, 13441),
    "nsp11": (13442, 13480),
    "nsp12": (13442, 16236),
    "nsp13": (16237, 18039),
    "nsp14": (18040, 19620),
    "nsp15": (19621, 20658),
    "nsp16": (20659, 21552),
    "S": (21563, 25384),
    "ns3a": (25393, 26220),
    "E": (26245, 26472),
    "M": (26523, 27191),
    "ns6": (27202, 27387),
    "ns7a": (27394, 27759),
    "ns7b": (27756, 27887),
    "ns8": (27894, 28259),
    "N": (28274, 29533),
    "ns10": (29558, 29674),
}

# Define a function to check for synonymous mutations
def is_synonymous(codon1, codon2):
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    aa1 = table.forward_table.get(codon1, 'Stop' if codon1 in table.stop_codons else '?')
    aa2 = table.forward_table.get(codon2, 'Stop' if codon2 in table.stop_codons else '?')
    return aa1 == aa2, aa1, aa2

# Function to determine gene based on position
def get_gene(position):
    for gene, (start, end) in gene_annotations.items():
        if start <= position < end:
            return gene
    return "Intergenic"

# Define the CSV file to write to
csv_file = "aa_frame_0.csv"

# Compare each sequence to the reference and write the results to the CSV file
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Sequence ID', 'Position', 'Ref Nucleotide', 'Seq Nucleotide', 'Ref Codon', 'Seq Codon', 'Ref AA', 'Seq AA', 'Type', 'Gene', 'AA Position'])

    for seq_record in sequences:
        current_gene = None
        aa_position = 1
        
        for i in range(0, len(reference.seq) - 3, 3):  # Adjusted frame: start at 1
            codon_start = i
            codon_ref = str(reference.seq[codon_start:codon_start + 3])
            codon_seq = str(seq_record.seq[codon_start:codon_start + 3])
            gene = get_gene(codon_start + 1)  # Determine the gene for the current position
            
            if gene != current_gene:
                current_gene = gene
                aa_position = 1  # Reset amino acid position counter for a new gene
            
            if codon_ref != codon_seq:
                for j in range(3):  # Check each nucleotide within the codon
                    ref_nucleotide = reference.seq[codon_start + j]
                    seq_nucleotide = seq_record.seq[codon_start + j]
                    
                    if ref_nucleotide != seq_nucleotide:
                        is_syn, aa_ref, aa_seq = is_synonymous(codon_ref, codon_seq)
                        writer.writerow([
                            seq_record.id, 
                            codon_start + j + 1, 
                            ref_nucleotide, 
                            seq_nucleotide, 
                            codon_ref, 
                            codon_seq, 
                            aa_ref, 
                            aa_seq, 
                            "Synonymous" if is_syn else "Non-synonymous", 
                            gene, 
                            aa_position
                        ])
            
            aa_position += 1  # Increment the amino acid position counter for the current gene

print(f"Mutations have been written to {csv_file}")
