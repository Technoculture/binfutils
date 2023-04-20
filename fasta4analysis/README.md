# DNA Sequences Data Analysis tool using fasta, txt, and csv files data to create histogram, barh and pie chart plots

# Package for fasta matplotlib Data Analysis

# Technoculture Matplotlib project task by Swati Mishra

# Requires Python 3.10 and above

# How to use the package and the script:
1. Change your current directory to the name of the package folder i.e. fasta4analysis
2. Place all your data.fasta and data files processed from it inside the above directory
3. pip install fasta4analysis or pip install .
4. Types of file formats to be used: fasta, txt, csv
5. fastacli is the cli command to invoke other commands
6. Commands to be used:

    1. '--help' command shows all command options and choices and purpose of the cli script
        fastacli --help
            Usage: fastacli [OPTIONS]

            DNA Sequences Data Analysis tool using fasta, txt, and csv files data to
            create histogram, barh and pie chart plots

            Options:
                -n, --files TEXT                Files
                -t, --type [fasta|csv|histo|barh|pie]
                                                Type of operation
                --help                          Show this message and exit.

        
    2. '-n' command is to receive file input types and '-t' command follows file input and is to receive the command choices for type of operation to be done on those input files. -n requires you to give three file inputs types .fasta, .txt and .csv followed by the command -t and command choice for operation. Using only -n command followed by three inputs invokes 'fasta' command choice by default.
    
    fastacli -n data.fasta -n data.txt -n sequences.csv

    3. 'fasta' command choice:
    It is used by writing fastacli command followed by '-n' command followed by three inputs file.fasta file.txt file1.csv followed by '-t' command followed by 'fasta'
    or is invoked by default without using '-t' command. It is used to extract large fasta files into a quickly accessible txt file.
    Once extracted, this command isn't needed again for other operations.
    
    e.g. fastacli -n data.fasta -n data.txt -n sequences.csv
    OR
    fastacli -n data.fasta -n data.txt -n sequences.csv -t fasta

    4. 'csv' command choice:
    It is used by writing fastacli command followed by '-n' command followed by inputs file.fasta file.txt file1.csv followed by '-t' command followed by 'csv'.
    It is used to generate two types of csv files, sequences_lengthwise and sequences_occurrence wise depending on certain prompted keyword choices.
    
    e.g. fastacli -n data.fasta -n data.txt -n filename.csv -t csv

    5. 'histo' command choice:
    It is used by writing fastacli command followed by '-n' command followed by inputs file1.csv followed by '-t' command followed by 'histo'.
    It is used to generate histogram plot using certain input prompts and the csv file containing lengthwise occurence data.
    Note: here you don't need to input three files. Just filename.csv is enough.
    
    e.g. fastacli -n lengthwise.csv -t histo

    6. 'barh' command choice:
    It is used by writing fastacli command followed by '-n' command followed by inputs file1.csv followed by '-t' command followed by 'barh'.
    It is used to generate barh plot using certain input prompts and the csv file containing sequence wise occurence data.
    Note: here you don't need to input three files. Just filename.csv is enough.
    
    e.g. fastacli -n sequencewise.csv -t barh

    7. 'pie' command choice:
    It is used by writing fastacli command followed by '-n' command followed by inputs file1.csv followed by '-t' command followed by 'pie'.
    It is used to generate pie-chart using certain input prompts and the csv file containing sequence wise occurence data.
    Note: here you don't need to input three files. Just filename.csv is enough.
    
    e.g. fastacli -n sequencewise.csv -t pie
