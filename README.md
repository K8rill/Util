# DNA/RNA and FASTQ Filtering — A Simple Tool for Working with Sequences

Welcome to **DNA/RNA Tools** — a set of tools for working with nucleic acid sequences and filtering data from fastq files. This project is designed for easy handling, analysis, and filtering of data.

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

    Input file:

    ```         
    @SRX079804:1:SRR292678:1:1101:21885:21885 1:N:0:1 BH:ok 
    ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA 
    +SRX079804:1:SRR292678:1:1101:21885:21885 1:N:0:1 BH:ok 
    FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD
    @SRX079804:1:SRR292678:1:1101:24563:24563 1:N:0:1 BH:failed
    ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG
    +SRX079804:1:SRR292678:1:1101:24563:24563 1:N:0:1 BH:failed
    BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D
    @SRX079804:1:SRR292678:1:1101:1266246:1266246 1:N:0:1 BH:failed
    AA
    +SRX079804:1:SRR292678:1:1101:1266246:1266246 1:N:0:1 BH:failed
    C@
    @SRX079804:1:SRR292678:1:1101:1269735:1269735 1:N:0:1 BH:failed
    C
    +SRX079804:1:SRR292678:1:1101:1269735:1269735 1:N:0:1 BH:failed
    G
    ```

    ``` python
    filter_fastq_file(input_fastq="input.fastq", output_fastq="filtered/filtered_output.fastq", gc_bounds=(40, 60), length_bounds=(10), quality_threshold=30)
     # Output: {'seq1': ('ATGCGTAC', 'IIIIIIII')}
    ```

    Output file

    ```         
    @SRX079804:1:SRR292678:1:1101:1266246:1266246 1:N:0:1 BH:failed
    AA
    +
    C@
    @SRX079804:1:SRR292678:1:1101:1269735:1269735 1:N:0:1 BH:failed
    C
    +
    G
    ```

### 1. **DNA and RNA Sequence Operations**

The project implements basic functions for sequence processing:

-   **transcribe**: Transcribes the sense DNA strand into RNA.
-   **reverse**: Returns the reversed sequence.
-   **complement**: Returns the complementary sequence.
-   **reverse_complement**: Returns the complementary and reversed sequence.
-   **gc_cont**: Calculates the GC content percentage.
-   **dna_or_rna**: Determines whether the sequence is DNA or RNA.

### 2. **FASTQ Data Filtering**

The `filter_fastq` function allows you to filter fastq sequences by the following criteria: - **GC content**: The range of guanine and cytosine content in the sequence. - **Length**: The length of the sequence. - **Quality**: The minimum average read quality on the phred33 scale.

## Bioinformatics File Handling

### 1. **Converting Multiline FASTA to Oneline**:

To convert a multiline FASTA file into a one-line FASTA format:

``` python
from bio_files_processor import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")
```

If no output file is provided, it will automatically generate one with "\_oneline" added to the input filename.

### 2. **Parsing BLAST Output**:

To extract top matches from BLAST output:

``` python
from bio_files_processor import parse_blast_output

parse_blast_output("example_blast_results.txt", "top_blast_hits.txt")
```

This will generate a sorted list of top matches from the BLAST output.

## Dependencies

The project uses only standard Python libraries.

## How to Contribute

If you have suggestions for improvement or new features to add, feel free to create a pull request or open an issue. We welcome any contributions!

------------------------------------------------------------------------
