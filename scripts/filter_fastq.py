def calculate_gc_content(sequence: str) -> float:
    """Calculates the GC content percentage in the sequence."""
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100


def calculate_average_quality(quality_str: str) -> float:
    """Calculates the average quality score."""
    return sum(ord(ch) - 33 for ch in quality_str) / len(quality_str)
