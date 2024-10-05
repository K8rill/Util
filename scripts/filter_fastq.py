def calculate_gc_content(sequence: str) -> float:
    """Вычисление процента содержания G и C в последовательности."""
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

def calculate_average_quality(quality_str: str) -> float:
    """Вычисление среднего значения качества."""
    return sum(ord(ch) - 33 for ch in quality_str) / len(quality_str)
