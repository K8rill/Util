# DNA/RNA and FASTQ Filtering — A Simple Tool for Working with Sequences

Welcome to **DNA/RNA Tools** — a set of tools for working with nucleic acid sequences and filtering data from fastq files. This project is designed for easy handling, analysis, and filtering of data using Python.

## Installation and Usage

1.  **Clone the repository**:

    ``` bash
    git clone https://github.com/your-repo/dna-rna-tools.git
    ```

2.  **Run**: To use the functions, write in your script:

    ``` python
    import main_script
    ```

    or

    ``` python
    from main_script import run_dna_rna_tools, filter_fastq
    ```

3.  **Example usage of the `run_dna_rna_tools` function**: You can use the `run_dna_rna_tools` function for DNA/RNA sequence operations:

    ``` python
    result = run_dna_rna_tools("ATGCGT", "transcribe")
    print(result)  # Output: AUGCGU
    ```

4.  **FASTQ data filtering**: Example usage of filtering with the `filter_fastq` function:

    ``` python
    seqs_example = {
        "seq1": ("ATGCGTAC", "IIIIIIII"),
        "seq2": ("GCGCGCGC", "IIIIIIHH"),
        "seq3": ("ATATATAT", "IIIIIIII"),
    }
    filtered_seqs = filter_fastq(seqs_example, gc_bounds=(40, 60), length_bounds=(6, 8), quality_threshold=30)
    print(filtered_seqs)  # Output: {'seq1': ('ATGCGTAC', 'IIIIIIII')}
    ```

## Functions and Features

### 1. **DNA and RNA Sequence Operations**

The project implements basic functions for sequence processing:

-   **transcribe**: Transcribes the sense DNA strand into RNA.
-   **reverse**: Returns the reversed sequence.
-   **complement**: Returns the complementary sequence.
-   **reverse_complement**: Returns the complementary and reversed sequence.
-   **gc_cont**: Calculates the GC content percentage.
-   **DNA_or_RNA**: Determines whether the sequence is DNA or RNA.

### 2. **FASTQ Data Filtering**

The `filter_fastq` function allows you to filter fastq sequences by the following criteria: 
-   **GC content**: The range of guanine and cytosine content in the sequence.
-   **Length**: The length of the sequence.
-   **Quality**: The minimum average read quality on the phred33 scale.

## Dependencies

The project uses only standard Python libraries.

## How to Contribute

If you have suggestions for improvement or new features to add, feel free to create a pull request or open an issue. We welcome any contributions!

------------------------------------------------------------------------