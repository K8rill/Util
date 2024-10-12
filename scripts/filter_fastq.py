from typing import Dict, Tuple
import os


def calculate_gc_content(sequence: str) -> float:
    """Calculates the GC content percentage in the sequence."""
    gc_count = sequence.upper().count("G") + sequence.upper().count("C")
    return (gc_count / len(sequence)) * 100


def calculate_average_quality(quality_str: str) -> float:
    """Calculates the average quality score."""
    return sum(ord(ch) - 33 for ch in quality_str) / len(quality_str)


def read_fastq(input_fastq: str) -> Dict[str, Tuple[str, str]]:
    """
    Reads a FASTQ file and converts it to a dictionary.
    Args:
        input_fastq: Path to the FASTQ file.
    Returns:
        A dictionary where keys are sequence names and values are
        tuples of (sequence, quality).
    """
    seqs = {}
    with open(input_fastq, "r") as file:
        while True:
            name = file.readline().strip()
            if name == "":
                print(name)
                break
            sequence = file.readline().strip()
            file.readline()  # plus line, ignore
            quality = file.readline().strip()
            seqs[name] = (sequence, quality)
    return seqs


def write_fastq(sequences: Dict[str, Tuple[str, str]], output_fastq: str):
    """
    Writes the filtered sequences to a FASTQ file.
    Args:
        sequences: Filtered sequences in dictionary format.
        output_fastq: Output FASTQ file path.
    """
    output_dir = os.path.dirname(output_fastq)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_fastq, "a+") as file:
        for name, (sequence, quality) in sequences.items():
            file.write(f"{name}\n")
            file.write(f"{sequence}\n")
            file.write("+\n")
            file.write(f"{quality}\n")
