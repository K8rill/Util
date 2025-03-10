from abc import ABC, abstractmethod
from typing import Any
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """
    Abstract class for nucleic acid.
    """

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index: Any):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

    @abstractmethod
    def validate_alphabet(self):
        """Checks whether the sequence contains only valid characters."""
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Base class for nucleic acids.
    """

    def __init__(self, sequence: str):
        self.sequence = sequence.upper()

    def complement(self):
        complement_map = str.maketrans("ATGCU", "TACGA")
        return self.sequence.translate(complement_map)

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.reverse().complement()

    def validate_alphabet(self):
        """Checks whether the sequence contains only valid nucleotides."""
        if not set(self.sequence).issubset({"A", "T", "G", "C", "U"}):
            raise ValueError("Invalid sequence: contains invalid characters.")


class DNASequence(NucleicAcidSequence):
    """
    Class for DNA sequences.
    """

    def complement(self):
        complement_map = str.maketrans("AGCU", "TCGA")
        return self.sequence.translate(complement_map)

    def validate_alphabet(self):
        if not set(self.sequence).issubset({"A", "T", "G", "C"}):
            raise ValueError(
                "Invalid DNA sequence: contains non-DNA characters."
                )

    def transcribe(self):
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    """
    Class for RNA sequences.
    """

    def validate_alphabet(self):
        if not set(self.sequence).issubset({"A", "U", "G", "C"}):
            raise ValueError(
                "Invalid RNA sequence: contains non-RNA characters."
                )


class AminoAcidSequence(BiologicalSequence):
    """
    Class for amino acid sequences.
    """

    def __init__(self, sequence: str):
        self.sequence = sequence.upper()

    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

    def validate_alphabet(self):
        if not set(self.sequence).issubset(self.valid_amino_acids):
            raise ValueError(
                "Invalid protein sequence: contains non-amino acid characters."
            )

    def molecular_weight(self):
        """Calculates the approximate molecular weight of the sequence."""
        weights = {
            "A": 89.09,
            "C": 121.15,
            "D": 133.10,
            "E": 147.13,
            "F": 165.19,
            "G": 75.07,
            "H": 155.16,
            "I": 131.17,
            "K": 146.19,
            "L": 131.17,
            "M": 149.21,
            "N": 132.12,
            "P": 115.13,
            "Q": 146.15,
            "R": 174.20,
            "S": 105.09,
            "T": 119.12,
            "V": 117.15,
            "W": 204.23,
            "Y": 181.19,
        }
        return sum(weights[aa] for aa in self.sequence)


def filter_fastq(
    input_file, output_file, min_length=50, min_quality=20, gc_bounds=(0, 100)
):
    """
    Filters the FASTQ file by GC composition, length, and average quality.
    Args:
        input_fastq: The path to the input FASTQ file.
        output_fastq: The path to the output FASTQ file.
        gc_bounds: Range of GC composition (percentages).
        length_bounds: The range of allowed read lengths.
        quality_threshold: Minimum average Phred score (Phred33).
    """
    with open(output_file, "w") as out_handle:
        filtered_records = []

        for record in SeqIO.parse(input_file, "fastq"):
            avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record)
            gc_content = gc_fraction(record.seq) * 100

            if (
                len(record) >= min_length
                and avg_quality >= min_quality
                and gc_bounds[0] <= gc_content <= gc_bounds[1]
            ):
                filtered_records.append(record)

        SeqIO.write(filtered_records, out_handle, "fastq")
