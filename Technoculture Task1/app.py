from Bio import SeqIO
from Bio.Seq import Seq
import csv


def melting_temp(dna_seq):
    dna_seq = dna_seq.upper()
    nA = dna_seq.count("A")
    nC = dna_seq.count("C")
    nG = dna_seq.count("G")
    nT = dna_seq.count("T")
    Tm = (nA + nT) * 2 + (nG + nC) * 4
    return round(Tm, 2)


def annealing_temp(aseq, bseq):
    a_temp = melting_temp(aseq)
    b_temp = melting_temp(bseq)
    return 0.3 * a_temp + 0.7 * b_temp - 14.9


def index_with_lowest_at(cdna):
    annealing_temp_diffs = []
    for i in range(2, len(cdna) - 1):
        a = cdna[:i]
        b = cdna[i:]
        a_temp = annealing_temp(a, a.transcribe())
        b_temp = annealing_temp(b, b.transcribe())
        annealing_temp_diffs.append(abs(a_temp - b_temp))
    return annealing_temp_diffs.index(min(annealing_temp_diffs))


def get_padlock_arms(miRNA):
    dna = miRNA.back_transcribe().reverse_complement()
    split = index_with_lowest_at(dna) + 2
    arm_a = dna[:split]
    arm_b = dna[split:]
    return arm_a, arm_b


def design_padlock_probe(common_csv_filename):
    padlock_probes = []

    with open(common_csv_filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header row
        for row in reader:
            miRNA_seq = row[0]
            miRNA = Seq(miRNA_seq)
            arm1, arm2 = get_padlock_arms(miRNA)
            padlock_probe = {
                'miRNA': miRNA_seq,
                'Arm1': str(arm1),
                'Arm2': str(arm2)
            }
            padlock_probes.append(padlock_probe)

    return padlock_probes


if __name__ == '__main__':
    common_csv_filename = 'final_filtered_common.csv'
    padlock_probes = design_padlock_probe(common_csv_filename)

    for probe in padlock_probes:
        print(f"miRNA: {probe['miRNA']}")
        print(f"Arm 1: {probe['Arm1']}")
        print(f"Arm 2: {probe['Arm2']}")
        print("------------")




