import pandas as pd
from Bio import SeqIO

# Convert "mature.fa" to CSV
fasta_file = "mature.fa"
csv_file_mature = "mature.csv"

records = []
for record in SeqIO.parse(fasta_file, "fasta"):
    records.append([record.id, str(record.seq)])

df_mature = pd.DataFrame(records, columns=["Name", "Sequence"])
df_mature.to_csv(csv_file_mature, index=False)
print("Conversion of mature.fa to CSV complete. File saved as mature.csv")

# Convert "mirna.tsv" to CSV
tsv_file = "mirna.tsv"
csv_file_mirna = "mirna.csv"

df_mirna = pd.read_csv(tsv_file, delimiter="\t")
df_mirna.to_csv(csv_file_mirna, index=False)
print("Conversion of mirna.tsv to CSV complete. File saved as mirna.csv")
