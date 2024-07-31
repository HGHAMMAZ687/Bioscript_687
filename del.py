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
    "nsp11": (13442, 13468),
    "nsp12": (13468, 16236),
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
    aa1 = table.forward_table.get(codon1, '?')
    aa2 = table.forward_table.get(codon2, '?')
    return aa1 == aa2, aa1, aa2

# Function to determine gene based on position
def get_gene(position):
    for gene, (start, end) in gene_annotations.items():
        if start <= position < end:
            return gene
    return "Intergenic"

# Define the CSV file to write to
csv_file = "del.csv"

# Compare each sequence to the reference and write the results to the CSV file
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Sequence ID', 'Position', 'Ref Nucleotide', 'Seq Nucleotide', 'Ref Codon', 'Seq Codon', 'Ref AA', 'Seq AA', 'Type', 'Gene'])
    
    for seq_record in sequences:
        i = 0
        while i < len(reference.seq):
            ref_nucleotide = reference.seq[i]
            seq_nucleotide = seq_record.seq[i]
            
            if ref_nucleotide != seq_nucleotide:
                if seq_nucleotide == "-":
                    # Handle deletions
                    start_del = i
                    while i < len(reference.seq) and seq_record.seq[i] == "-":
                        i += 1
                    end_del = i
                    writer.writerow([seq_record.id, f"del{start_del + 1}_{end_del}", reference.seq[start_del:end_del], "-", "-", "-", "-", "-", "Deletion", get_gene(start_del)])
                else:
                    # Handle SNPs
                    codon_start = (i // 3) * 3
                    codon_ref = str(reference.seq[codon_start:codon_start + 3])
                    codon_seq = str(seq_record.seq[codon_start:codon_start + 3])
                    is_syn, aa_ref, aa_seq = is_synonymous(codon_ref, codon_seq)
                    gene = get_gene(i)
                    writer.writerow([seq_record.id, i + 1, ref_nucleotide, seq_nucleotide, codon_ref, codon_seq, aa_ref, aa_seq, "Synonymous" if is_syn else "Non-synonymous", gene])
            i += 1

print(f"Mutations have been written to {csv_file}")
