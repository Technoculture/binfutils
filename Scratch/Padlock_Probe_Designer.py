'''
1. Updated the CLI code to incorporate the desired modification in the melting temperature calculation.
2. Explored the Bio.SeqUtils.MeltingTemp module and gained an understanding of the TmNN formula.
3. Implemented the TmNN formula in the code to calculate the melting temperature accurately.
4. Conducted rigorous testing to validate the correctness and reliability of the new melting temperature calculation.
5. Ongoing work on further refinements and optimizations to improve the overall performance of the CLI.

'''

import click
import csv
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp


def calculate_melting_temp(dna_seq):
    return MeltingTemp.Tm_NN(dna_seq)


def melting_temp(dna_seq):
    return round(calculate_melting_temp(dna_seq), 2)


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
@click.argument('input-file', type=click.Path(exists=True))
@click.option('--output-file', default='padlock_probes.csv', help='Output CSV file')
def design_padlock_probe(input_file, output_file):
    target_sequences = []
    target_data = []

    if input_file.endswith('.csv'):
        with open(input_file, 'r') as file:
            reader = csv.DictReader(file)
            target_data = list(reader)
            target_sequences = [row['Sequence'] for row in target_data]
    elif input_file.endswith(('.fa', '.fasta', '.txt')):
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
