from Bio.Seq import Seq
import click
import csv


def melting_temp(dna_seq):
    dna_seq = dna_seq.upper()
    nA = dna_seq.count("A")
    nC = dna_seq.count("C")
    nG = dna_seq.count("G")
    nU = dna_seq.count("U")

    Tm = (nA + nU) * 2 + (nG + nC) * 4
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


@click.command()
@click.option('--reporter-seq', default='ACGT', help='Reporter sequence')
@click.option('--output-file', default='padlock_probes.csv', help='Output CSV file')
def design_padlock_probe(reporter_seq, output_file):
    target_file = 'final_filtered_common.csv'
    target_sequences = []
    target_data = []

    with open(target_file, 'r') as file:
        reader = csv.DictReader(file)
        target_data = list(reader)
        target_sequences = [row['Sequence'] for row in target_data]

    padlock_probes = []
    for i, target_seq in enumerate(target_sequences):
        miRNA = Seq(target_seq.upper())
        arm1, arm2 = get_padlock_arms(miRNA)
        res = str(arm2) + reporter_seq + str(arm1)
        padlock_probe = target_data[i].copy()
        padlock_probe['Padlock_Probe'] = res
        padlock_probe['Arm1'] = str(arm1)
        padlock_probe['Arm2'] = str(arm2)
        padlock_probe['Melting_Temp'] = melting_temp(target_seq)
        padlock_probe['Annealing_Temp'] = annealing_temp(arm1, arm2)
        padlock_probes.append(padlock_probe)

    fieldnames = list(target_data[0].keys())
    sequence_index = fieldnames.index('Sequence')
    fieldnames.insert(sequence_index + 1, 'Melting_Temp')
    fieldnames.insert(sequence_index + 2, 'Annealing_Temp')
    fieldnames.insert(sequence_index + 3, 'Padlock_Probe')
    fieldnames.insert(sequence_index + 4, 'Arm1')
    fieldnames.insert(sequence_index + 5, 'Arm2')

    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(padlock_probes)


if __name__ == '__main__':
    design_padlock_probe()


