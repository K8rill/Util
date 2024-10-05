def transcribe(sequence: str) -> str:
    """
    Transcribes the coding DNA sequence into RNA.
    """
    return sequence.replace('T', 'U').replace('t', 'u')


def reverse(sequence: str) -> str:
    """
    Returns the reversed sequence.
    """
    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Returns the complementary sequence.
    """
    complement_map = str.maketrans('ATGCUatgcu', 'TACGAtacga')
    return sequence.translate(complement_map)


def reverse_complement(sequence: str) -> str:
    """
    Returns the reverse complement of the sequence.
    """
    return reverse(complement(sequence))


def gc_cont(sequence: str) -> float:
    """
    Calculates the GC content percentage in the sequence.
    """
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100


def DNA_or_RNA(sequence: str) -> str:
    """
    Determines whether the sequence is DNA or RNA.
    """
    if ('U' in sequence.upper()) == ('T' in sequence.upper()):
        return 'Unknown'
    elif 'U' in sequence.upper():
        return 'RNA'
    elif 'T' in sequence.upper():
        return 'DNA'
