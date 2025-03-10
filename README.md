# DNA/RNA and FASTQ Filtering — A Simple Tool for Working with Sequences

Welcome to **DNA/RNA Tools** — a set of tools for working with nucleic acid sequences and filtering data from fastq files. This project is designed for easy handling, analysis, and filtering of data.

## Installation and Usage

**Clone the repository**:

``` bash
git clone git@github.com:K8rill/Util.git
```

**Run**: To use the functions, write in your script:

``` python
import main_script
```

or

``` python
from main_script import DNASequence, RNASequence, AminoAcidSequence, filter_fastq
```

# Usage

## Sequence Handling

Import the necessary classes and create instances:

``` python
# DNA sequence operations
dna_seq = DNASequence("ATGCGT")
print(dna_seq.complement())  
# Output: TACGCA
print(dna_seq.reverse_complement())  
# Output: ACGCAT
rna_seq = dna_seq.transcribe()
print(rna_seq)  
# Output: AUGCGU

# Amino acid sequence operations
protein_seq = AminoAcidSequence("MKAT")
print(protein_seq.molecular_weight())  
# Output: 503.61
```

## FASTQ Filtering

Filter FASTQ sequences based on length, GC content, and quality:

``` python
filter_fastq(
    input_file="input.fastq",
    output_file="filtered.fastq",
    min_length=50,
    min_quality=20,
    gc_bounds=(40, 60)
)
```

## Dependencies

biopython

## How to Contribute

If you have suggestions for improvement or new features to add, feel free to create a pull request or open an issue. We welcome any contributions!

------------------------------------------------------------------------
