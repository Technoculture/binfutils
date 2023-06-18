# How to Use the Package and the Script

1. Change your current directory to the name of the package folder, i.e. `fasta4analysis`.
2. Place all your `data.fasta` and processed `data` files inside the package directory.
3. Install the package using pip: `pip install fasta4analysis` or `pip install .`
4. Supported file formats: `fasta`, `txt`, `csv`
5. Use `fastacli` to invoke other commands.
6. Available commands:

    1. `--help`: displays all command options, choices, and the purpose of the CLI script.

            Usage: fastacli [OPTIONS]
            
            DNA Sequences Data Analysis tool using fasta, txt, and csv files data to create histogram, barh and pie chart plots
            
            Options:
                -n, --files TEXT                Files
                -t, --type [fasta|csv|histo|barh|pie]
                                                Type of operation
                --help                          Show this message and exit.

        
    2. `-n` command: receive file input types.
    3. `-t` command: follows file input and receives the command choices for the type of operation to be done on those input files. Use `-n` to give three file inputs types `.fasta`, `.txt`, and `.csv` followed by `-t` and command choice for operation. Using only `-n` invokes the `fasta` command choice by default.
    
            fastacli -n data.fasta -n data.txt -n sequences.csv

    4. `fasta` command choice: used to extract large fasta files into a quickly accessible txt file. Once extracted, this command isn't needed again for other operations. Invoke it by writing `fastacli -n data.fasta -n data.txt -n sequences.csv -t fasta` or use the default by only using `-n` command.
    
    e.g. `fastacli -n data.fasta -n data.txt -n sequences.csv`
    

    5. `csv` command choice: used to generate two types of csv files, `sequences_lengthwise` and `sequences_occurrence` wise, depending on certain prompted keyword choices. Invoke it by writing `fastacli -n data.fasta -n data.txt -n filename.csv -t csv`.
    
    e.g. `fastacli -n data.fasta -n data.txt -n filename.csv -t csv`

    6. `histo` command choice: used to generate a histogram plot using certain input prompts and the csv file containing lengthwise occurrence data. Note: you don't need to input three files. Just `filename.csv` is enough. Invoke it by writing `fastacli -n lengthwise.csv -t histo`.
    
    e.g. `fastacli -n lengthwise.csv -t histo`

    7. `barh` command choice: used to generate a barh plot using certain input prompts and the csv file containing sequence-wise occurrence data. Note: you don't need to input three files. Just `filename.csv` is enough. Invoke it by writing `fastacli -n sequencewise.csv -t barh`.
    
    e.g. `fastacli -n sequencewise.csv -t barh`

    8. `pie` command choice: used to generate a pie-chart using certain input prompts and the csv file containing sequence-wise occurrence data. Note: you don't need to input three files. Just `filename.csv` is enough. Invoke it by writing `fastacli -n sequencewise.csv -t pie`.
    
    e.g. `fastacli -n sequencewise.csv -t pie`
