## Overview
A FASTA file parse written in python

## REGEX Use
- Used r"^>" to find the start of a sequence.
- Used r"T" to check if T is present in a sequence, and r'U' to check if  U is present in a sequence. If only T is present, it's a valid DNA, and if only U is present, it's a valid RNA. If both or none are present, it's an invalid sequence.
- Used r"ATG" to find the beginning of a DNA ORF, and r"TAG|TAA|TGA" to find the end. Now, if for all pairs of beginning and end of an ORF where beginning + 2 < end, it is a valid ORF. (We add 2 to make sure there is no overlap between beginning and end ORFs)
- Used r"AUG" to find the beginning of an RNA ORF, and r"UAG|UAA|UGA" to find the end. Now, if for all pairs of beginning and end of an ORF where beginning + 2 < end, it is a valid ORF. (We add 2 to make sure there is no overlap between beginning and end ORFs)

## Assumptions
- Nucleotides are always depicted with Upper case letters.
- Any necluotide not in (A, U, G, T, C) is ambiguous.
- A valid DNA or RNA sequence may or may not have any valid ORFs.

## Output Example
\>SEQ_DNA_6\
ATGAAACCCGGGTTTTAA\
Type: DNA\
Length: 18\
Counts: A=6, C=3, G=4, T=5\
Found 1 ORF(s):\
\- Start at 0, Stop at 15, ORF: ATGAAACCCGGGTTTTAA\

--SUMMARY--\
Valid DNA sequences: 6\
Mean DNA length: 10.17\
Valid RNA sequences: 4\
Mean RNA length: 12.25\
Valid Invalid sequences: 2\

## Requirements
- Python 3.8+

## Usage
The program accepts FASTA file path as an argument. To run the code, navigate to the root directory and execute the main script with the FAST file path as shown below:
```bash
python3 main.py fasta_file_path
```

## Results
The results of the experiments can be found in the `results.txt` file.
For all three files, the output is divided by ----file_name---- lines.