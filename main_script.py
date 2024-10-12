import scripts.dna_rna_tools as nk_tools
import scripts.filter_fastq as fi_fa
import os
from typing import List, Union, Dict, Tuple


def run_dna_rna_tools(*args: str) -> Union[str, List[str]]:
    """
    Executes the specified operation on the given sequences.
    Args: Arbitrary number of strings with sequences,
    the last argument should be the operation name
    (transcribe, reverse, complement, reverse_complement, gc_cont, DNA_or_RNA).
    Returns:
    If one sequence is passed — str,
    If multiple sequences are passed — list[str, str].
    """
    operation = args[- 1]
    res = []
    for i in range(len(args) - 1):
        if operation == 'transcribe':
            res.append(nk_tools.transcribe(args[i]))
        elif operation == 'reverse':
            res.append(nk_tools.reverse(args[i]))
        elif operation == 'complement':
            res.append(nk_tools.complement(args[i]))
        elif operation == 'reverse_complement':
            res.append(nk_tools.reverse_complement(args[i]))
        elif operation == 'gc_cont':
            res.append(nk_tools.gc_cont(args[i]))
        elif operation == 'dna_or_rna':
            res.append(nk_tools.dna_or_rna(args[i]))
    if len(res) == 1:
        return res[0]
    else:
        return res


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0
) -> Dict[str, Tuple[str, str]]:
    """
    Filters sequences by GC content, length, and average quality.
    Args:
    seqs: Dictionary of sequences, where the key is the sequence name,
          and the value is a tuple (sequence, quality)
    gc_bounds: Range of GC content in percentage
    length_bounds: Range of sequence lengths
    quality_threshold: Minimum average read quality (Phred33 scale)
    Returns:
        Filtered dictionary of sequences that meet all criteria.
    """
    seqs = fi_fa.read_fastq(input_fastq)
    filtered_seqs = {}
    # Adjust boundaries if a single number is provided
    if isinstance(gc_bounds, (float, int)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    for name, (sequence, quality) in seqs.items():
        # Filter by length
        seq_length = len(sequence)
        if not (length_bounds[0] <= seq_length <= length_bounds[1]):
            continue
        # Filter by GC content
        gc_content = fi_fa.calculate_gc_content(sequence)
        if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
            continue
        # Filter by average quality
        avg_quality = fi_fa.calculate_average_quality(quality)
        if avg_quality < quality_threshold:
            continue
        # If all conditions are met, add the sequence
        filtered_seqs[name] = (sequence, quality)
    output_fastq = os.path.join("filtered", output_fastq)
    fi_fa.write_fastq(filtered_seqs, output_fastq)
    return filtered_seqs
