import pandas as pd
import urllib.request
import gzip
import csv
from Bio import SeqIO

# Download "mature.fa" from miRBase
miRBase_url = "https://www.mirbase.org/ftp/CURRENT/mature.fa.gz"
miRBase_file = "mature.fa.gz"
urllib.request.urlretrieve(miRBase_url, miRBase_file)

# Extract the downloaded gzip file
extracted_file = "mature.fa"
with gzip.open(miRBase_file, 'rb') as gz_file:
    with open(extracted_file, 'wb') as out_file:
        out_file.write(gz_file.read())

print("Download and extraction of mature.fa complete. File saved as mature.fa")

# Download "mirna.tsv" from salivaDB
salivaDB_url = "https://webs.iiitd.edu.in/raghava/salivadb/mirna.tsv"
salivaDB_file = "mirna.tsv"
urllib.request.urlretrieve(salivaDB_url, salivaDB_file)

print("Download of mirna.tsv complete. File saved as mirna.tsv")

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

# Read the mirna.csv file
csv_file_mirna = "mirna.csv"
df_mirna = pd.read_csv(csv_file_mirna)

# Filter rows based on Disease_Category
df_filtered = df_mirna[df_mirna['Disease_Category'] != 'Healthy']

# Save the filtered DataFrame to a new CSV file
csv_file_filtered = "mirna_filtered.csv"
df_filtered.to_csv(csv_file_filtered, index=False)
print("Filtered data saved to mirna_filtered.csv")

# Read the mature.csv and mirna_filtered.csv files
csv_file_mature = "mature.csv"
csv_file_filtered = "mirna_filtered.csv"

df_mature = pd.read_csv(csv_file_mature)
df_filtered = pd.read_csv(csv_file_filtered)

# Find the common data based on all columns
common_df = pd.merge(df_mature, df_filtered)

# Apply the first filter to select desired columns
filtered_df = common_df[['Name', 'Biomarker_Name', 'Sequence', 'Sampling_Method', 'Collection_site', 'Disease_name', 'Regulation', 'Exosomal/Non-exosomal']]

# Apply the second filter to keep rows with the same 'Name' and 'Biomarker_Name'
filtered_df = filtered_df[filtered_df['Name'] == filtered_df['Biomarker_Name']]

# Apply the third filter to remove the 'Name' column
filtered_df = filtered_df.drop('Name', axis=1)

# Keep only one row for each unique 'Biomarker_Name'
filtered_df = filtered_df.drop_duplicates(subset='Biomarker_Name')

# Save the final filtered DataFrame to a new CSV file
csv_file_final_filtered = "filtered_data.csv"
filtered_df.to_csv(csv_file_final_filtered, index=False)
print("Final filtered data saved to filtered_data.csv")
