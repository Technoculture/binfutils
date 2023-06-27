# Technoculture_Task1
# Padlock Probe Designer CLI
The Padlock Probe Designer CLI is a command-line tool developed in Python using the Click library and Biopython module. It enables users to design Padlock Probes by automating the process of downloading data from miRBase and salivaDB, converting them into CSV format, performing analysis, filtering, and generating Padlock Probes. This README provides a detailed overview of the project and guides users on how to install, set up, and utilize the CLI effectively.

# Installation
1. Clone the repository to your local machine:
```
git clone https://github.com/palak178/Technoculture_Task1.git 
```

2. Navigate to the project directory: 
```
cd Technoculture_Task1
```

3. Install the required dependencies using pip:
```
pip install -r requirements.txt
```

# Usage
1. Run the following commands to execute the Padlock Probe Designer CLI:
```
python input_data_preparation.py
``` 
```
python Padlock_Probe_Designer.py --output-file padlock_probes.csv filtered_data.csv
```

   The CLI will utilize the provided filtered_data.csv file as input and generate Padlock Probes.

2. Run the following commands to execute the Padlock Probe Designer CLI if having fasta file as input:
```
python Padlock_Probe_Designer.py --output-file padlock_probes.csv input.fa
```

# Input Data Preparation
   Before designing padlock probes, the input data needs to be prepared using the input_data_preparation.py script, which performs the following tasks:

   1. Retrieves the mature.fa file from miRBase and mirna.tsv file from salivaDB.
  
   2. Converts the downloaded data to CSV format.
  
   3. Analyze and filters the miRNA data based on Disease Categories.
   
   4. Refines the data and generates the 'filtered_data.csv' file as input for the Padlock Probe Designer CLI.

 To prepare the input data, run the following command:
 ```
  python input_data_preparation.py
 ```
This will download the required files and generate a filtered_data.csv file, which serves as the input file for the padlock probe designer.

# Padlock Probe Design
The Padlock_Probe_Designer.py script designs padlock probes based on the data in the filtered_data.csv file. It calculates melting temperatures, annealing temperatures, and generates padlock probe sequences for each target miRNA.

To design padlock probes, run the following command:
```
python Padlock_Probe_Designer.py --output-file <output_filename.csv> <input_filename.csv>
```
Replace <output_filename.csv> with the desired name for the output file that will contain the designed padlock probes.

Replace <input_filename.csv> with the desired name for the input file that can be in csv or fasta format.

The script will generate the output file containing the designed padlock probes, along with relevant information such as melting temperatures, annealing temperatures, and arm sequences.

# Contributing
Contributions are welcome! If you'd like to contribute to this project, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature/bug fix.
3. Make the necessary changes and commit them.
4. Push the changes to your branch.
5. Open a pull request, describing the changes you made.

# Support
If you have any questions, issues, or need assistance, please feel free to reach out at palak.namdev1701@gmail.com

**Thank you for using the Padlock Probe Designer CLI!**


