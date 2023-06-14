from flask import Flask
from flask_restful import Api, Resource, reqparse
from Bio.Seq import Seq
from flask_jsonpify import jsonify
import csv

app = Flask(__name__)
api = Api(app)

def load_common_csv(file_path):
    sequences = []
    with open(common.csv, 'r') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            sequences.append(row['Sequence'])
    return sequences

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

class PadlockProbeGenerator(Resource):
    def __init__(self):
        self.sequences = load_common_csv('common.csv')

    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('Biomarker_Name', type=str, help='Target sequence is required', required=True)
        parser.add_argument('Disease_Category', type=str, help='Reporter sequence is required', required=True)
        args = parser.parse_args()

        target = args['Biomarker_Name']
        reporter_seq = args['Disease_Category']

        if target not in self.sequences:
            return jsonify({'status': 400, 'message': 'Invalid target sequence'})

        miRNA = Seq(target.upper())
        arm1, arm2 = get_padlock_arms(miRNA)
        res = str(arm2) + reporter_seq + str(arm1)

        return jsonify({'status': 200, 'output': res})

class HelloWorld(Resource):
    def get(self):
        return 'Hello World!'

api.add_resource(PadlockProbeGenerator, '/padlockGen')
api.add_resource(HelloWorld, '/')

if __name__ == '__main__':
    app.run(debug=True)
