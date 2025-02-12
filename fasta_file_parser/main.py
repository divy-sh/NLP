import sys
import re
from collections import Counter
from typing import List, Tuple, Dict

# All required regex patterns, compiled for efficiency
SEQUENCE_START_REGEX = re.compile(r"^>")
DNA_REGEX = re.compile(r"T")
RNA_REGEX = re.compile(r"U")
DNA_START_ORF_REGEX = re.compile(r"ATG")
DNA_END_ORF_REGEX = re.compile(r"TAG|TAA|TGA")
RNA_START_ORF_REGEX = re.compile(r"AUG")
RNA_END_ORF_REGEX = re.compile(r"UAA|UAG|UGA")

# main function which reads the input file, processes the sequences, and prints the results
def main(file_path: str):
    # Dictionary to store sequence metadata - sequnece type count, and total length sum so that we can calculate average length
    sequenceMetadata = {"DNA": {"lenSum": 0, "count": 0}, "RNA": {"lenSum": 0, "count": 0}, "Invalid": {"lenSum": 0, "count": 0}}

    for sequence_id, sequence_data in parse_fasta(file_path):
        printSequenceDetails(sequence_id, sequence_data, sequenceMetadata)
    
    printSummary(sequenceMetadata)

def parse_fasta(filepath: str) -> Tuple[str, str]:
    # Read a FASTA file and process each sequence
    with open(filepath, "r") as fasta_file:
        sequence_id = ""
        sequence_data = []
        for line in fasta_file:
            line = line.strip()
            # If new sequnece starts, yield the previous sequence
            if SEQUENCE_START_REGEX.match(line): # Header line
                if sequence_data:
                    yield sequence_id, ''.join(sequence_data)
                sequence_id = line
                sequence_data = []
            else:
                sequence_data.append(line)
        if sequence_data:  # Process the last sequence
            yield sequence_id, ''.join(sequence_data)

def getSequenceType(sequence: str) -> str:
    # Returns sequence type - DNA, RNA, or Invalid
    if DNA_REGEX.search(sequence) and RNA_REGEX.search(sequence):
        return "Invalid"
    elif DNA_REGEX.search(sequence):
        return "DNA"
    elif RNA_REGEX.search(sequence):
        return "RNA"
    return "Invalid"

def find_ORFs(sequence: str, start_regex, end_regex) -> List[Tuple[int, int, str]]:
    # Finds ORFs and calculates the start and end positions
    orfs = []
    start_positions = list(start_regex.finditer(sequence))
    end_positions = list(end_regex.finditer(sequence))

    # If the start position is before the end position, and there is enough space between 
    # the two for a valid start codon, then it is a valid ORF
    for start in start_positions:
        for end in end_positions:
            if start.start() + 2 < end.start():
                orfs.append((start.start(), end.start(), sequence[start.start():end.start() + 3]))
    return orfs

def getORFs(sequence: str, sequence_type: str) -> List[Tuple[int, int, str]]:
    # Returns ORFs for a given sequence
    if sequence_type == "DNA":
        return find_ORFs(sequence, DNA_START_ORF_REGEX, DNA_END_ORF_REGEX)
    elif sequence_type == "RNA":
        return find_ORFs(sequence, RNA_START_ORF_REGEX, RNA_END_ORF_REGEX)
    return []

def computeCounts(sequence: str) -> str:
    # Returns counts of each nucleotide in the sequence
    count = Counter(sequence)
    # Get core nucleotides with zero count from the counter, everything else is ambiguous
    counts = [f"{base}={count[base]}" for base in "ACGTU" if count[base] > 0]
    # calculate ambiguous count
    ambiguous_count = len(sequence) - sum(count[base] for base in "ACGTU")
    if ambiguous_count > 0:
        counts.append(f"Ambiguous={ambiguous_count}")
    return ", ".join(counts)

def printSequenceDetails(sequence_id: str, sequence_data: str, sequenceMetadata: Dict[str, Dict[str, int]]):
    # Process and print sequence details
    print(sequence_id)
    print(sequence_data)
    
    sequence_type = getSequenceType(sequence_data)
    
    # Update sequence metadata
    sequenceMetadata[sequence_type]["count"] += 1
    sequenceMetadata[sequence_type]["lenSum"] += len(sequence_data)
    
    print(f"Type: {sequence_type}")
    print(f"Length: {len(sequence_data)}")
    print(f"Counts: {computeCounts(sequence_data)}")
    
    # We only find ORFs for valid DNA and RNA sequences
    if sequence_type != "Invalid":
        ORF_data = getORFs(sequence_data, sequence_type)
        if ORF_data:
            print(f"Found {len(ORF_data)} ORF(s):")
            for start, end, orf in ORF_data:
                print(f"- Start at {start}, Stop at {end}, ORF: {orf}")
    print()

def printSummary(sequenceMetadata: Dict[str, Dict[str, int]]):
    """Print a summary of sequence counts and lengths."""
    print("--SUMMARY--")
    # Print counts and mean lengths for DNA, RNA, and Invalid sequences
    for seq_type in ["DNA", "RNA"]:
        count = sequenceMetadata[seq_type]["count"]
        total_length = sequenceMetadata[seq_type]["lenSum"]
        mean_length = 0 if count == 0 else total_length / count
        print(f"Valid {seq_type} sequences: {count}")
        print(f"Mean {seq_type} length: {mean_length:.2f}")

    print(f"Invalid sequences: {sequenceMetadata['Invalid']['count']}")

if __name__ == "__main__":
    # Take the FASTA file path as an argument
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
        main(file_path)
    else:
        print("No file path was passed.")