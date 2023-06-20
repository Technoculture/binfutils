import click
import csv
import random
from Bio import SeqIO
from Bio.Seq import Seq


def short_formula(dna_seq):
    """
    Melting temperature calculation for a sequence strictly shorter than 14 nucleotides
    """
    dna_seq = dna_seq.upper()
    nA = dna_seq.count("A")
    nC = dna_seq.count("C")
    nG = dna_seq.count("G")
    nU = dna_seq.count("U")
    nT = dna_seq.count("T")

    Tm = (nA + nU + nT) * 2 + (nG + nC) * 4
    return round(Tm, 2)


def long_formula(dna_seq):
    """
    Melting temperature calculation for a sequence strictly longer than 13 nucleotides
    """
    dna_seq = dna_seq.upper()
    nA = dna_seq.count("A")
    nC = dna_seq.count("C")
    nG = dna_seq.count("G")
    nT = dna_seq.count("T")
    nU = dna_seq.count("U")

    Tm = 64.9 + 41 * (nG + nC - 16.4) / (nA + nU + nT + nG + nC)
    return round(Tm, 2)


def melting_temp(dna_seq):
    if len(dna_seq) < 14:
        return short_formula(dna_seq)
    else:
        return long_formula(dna_seq)


def annealing_temp(aseq, bseq):
    a_temp = melting_temp(aseq)
    b_temp = melting_temp(bseq)
    return max(0.3 * a_temp + 0.7 * b_temp - 14.9, 0)


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
    split = index_with_lowest_at(dna) + 3
    arm_a = dna[:split]
    arm_b = dna[split:]
    return arm_a, arm_b


def get_reporter_sequence():
    reporter_sequences = [
        "TTCCTTTTACGACCTCAATGCTGCTGCTGTACTACTCTT",
        "TTCCTTTTACGATCGCGCTTGGTATAATCGCTACTTCTT"
    ]
    return random.choice(reporter_sequences)


@click.command()
@click.option('--output-file', default='padlock_probes.csv', help='Output CSV file')
@click.argument('input-file', type=click.Path(exists=True))
def design_padlock_probe(output_file, input_file):
    target_sequences = []
    target_data = []

    if input_file.endswith('.csv'):
        with open(input_file, 'r') as file:
            reader = csv.DictReader(file)
            target_data = list(reader)
            target_sequences = [row['Sequence'] for row in target_data]
    elif input_file.endswith('.fa'):
        with open(input_file, 'r') as file:
            sequences = SeqIO.parse(file, 'fasta')
            for sequence in sequences:
                target_sequences.append(str(sequence.seq))
                target_data.append({'Sequence': str(sequence.seq)})
    else:
        click.echo('Unsupported file format. Only CSV and FASTA files are supported.')
        return

    padlock_probes = []
    for i, target_seq in enumerate(target_sequences):
        miRNA = Seq(target_seq.upper())
        arm1, arm2 = get_padlock_arms(miRNA)
        reporter_seq = get_reporter_sequence()
        res = str(arm2) + reporter_seq.lower() + str(arm1)
        padlock_probe = target_data[i].copy()
        padlock_probe['Padlock_Probe'] = res
        padlock_probe['Arm1'] = arm1
        padlock_probe['Arm2'] = arm2
        padlock_probe['Arm1_Melting_Temp'] = melting_temp(str(arm1))
        padlock_probe['Arm2_Melting_Temp'] = melting_temp(str(arm2))
        padlock_probe['Arm1_Annealing_Temp'] = annealing_temp(str(arm1), reporter_seq)
        padlock_probe['Arm2_Annealing_Temp'] = annealing_temp(str(arm2), reporter_seq)
        padlock_probe['Melting_Temp'] = melting_temp(target_seq)
        padlock_probe['Annealing_Temp'] = annealing_temp(target_seq, reporter_seq)
        padlock_probes.append(padlock_probe)

    fieldnames = list(target_data[0].keys())
    fieldnames.insert(fieldnames.index('Sequence') + 1, 'Melting_Temp')
    fieldnames.insert(fieldnames.index('Melting_Temp') + 1, 'Annealing_Temp')
    fieldnames.insert(fieldnames.index('Annealing_Temp') + 1, 'Arm1')
    fieldnames.insert(fieldnames.index('Arm1') + 1, 'Arm1_Melting_Temp')
    fieldnames.insert(fieldnames.index('Arm1_Melting_Temp') + 1, 'Arm1_Annealing_Temp')
    fieldnames.insert(fieldnames.index('Arm1_Annealing_Temp') + 1, 'Arm2')
    fieldnames.insert(fieldnames.index('Arm2') + 1, 'Arm2_Melting_Temp')
    fieldnames.insert(fieldnames.index('Arm2_Melting_Temp') + 1, 'Arm2_Annealing_Temp')
    fieldnames.insert(fieldnames.index('Arm2_Annealing_Temp') + 1, 'Padlock_Probe')

    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(padlock_probes)


if __name__ == '__main__':
    design_padlock_probe()
