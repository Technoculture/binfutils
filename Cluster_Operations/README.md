# Fasta Sequence Clustering Analysis

The Jupyter Notebook '*final_cluster.ipynb*' analyzes a fasta file containing sequences for clustering purposes. It utilizes the Biopython module to calculate the total number of sequences present in the fasta file and the sum of cluster sequence counts. It also verifies if the total number of lines in the fasta file matches the sum of the cluster sequence count. The code saves the cluster number and cluster sequence count in a CSV file and saves the cluster in a clusters.clstr file. Finally, it finds a subsequence of length 14 from the cluster with the highest number of sequences.

# Prerequisites

  1. Python 3.9

  2. Biopython module (install using pip install biopython)

# Instructions
  
  1. Ensure you have the necessary prerequisites installed.
  
  2. Download the fasta file that you want to analyze and place it in the same directory as this Jupyter Notebook.
 
  3. Open the Jupyter Notebook in JupyterLab or Jupyter Notebook.
  
  4. Modify the fasta_file_path variable to specify the filename and path of your fasta file.
  
  5. Run the code cells in the Jupyter Notebook.

# Note

  1. Replace "your_fasta_file.fasta" in the fasta_file_path variable with the actual filename and path of your fasta file.
  
  2. The CSV file will be saved as "clusters.csv" in the same directory as the Jupyter Notebook.
  
  3. The clusters.clstr file will be saved as "clusters.clstr" in the same directory as the Jupyter Notebook.
  
  4. The subsequence of length 14 from the cluster with the highest number of sequences will be printed on the console.





