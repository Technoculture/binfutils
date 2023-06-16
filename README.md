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
1. Run the following command to execute the Padlock Probe Designer CLI:
```
python Padlock_Probe_Designer.py --input-file final_filtered_common.csv 
```
2. The CLI will utilize the provided final_filtered_common.csv file as input and generate Padlock Probes.

# Input Data Preparation
1. Ensure you have the following files available in the project directory: 

     a. mature.fa from miRBase: This file contains mature microRNA sequences.
     
     b. mirna.tsv from salivaDB: This file contains miRNA data related to salivary biomarkers.

2. Execute the following steps to prepare the input data:

    a. Run the Python script to download the required files and convert them into CSV format using Biopython and Pandas:
    ```
    python main.py
    python convert.py
    ```
   This script retrieves the mature.fa file from miRBase and mirna.tsv file from salivaDB, converting both files into mature.csv and mirna.csv, respectively.

    b. Run the script filter.py to analyze and filter the mature.csv and mirna.csv files, extracting the common data into a new CSV file named common.csv:
    ```
    python filter.py
    python compare.py
    ``` 
    This step compares the sequences in mature.csv and mirna.csv, identifying the common data and saving it as common.csv.

    c. Finally, filter the common.csv file using Pandas to obtain the final filtered data, and save it as final_filtered_common.csv:
    ```
    python filter_common.py
    ```
    This step applies additional filters to the common.csv file, refining the data and generating the final_filtered_common.csv file as input for the Padlock
    Probe Designer CLI.

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


