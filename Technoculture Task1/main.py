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
