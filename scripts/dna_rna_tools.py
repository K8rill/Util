
def transcribe(sequence: str) -> str:
    """
    Транискрипция со смысловой ДНК в РНК.
    """
    return sequence.replace('T', 'U').replace('t', 'u')


def reverse(sequence: str) -> str:
    """
    Возвращает перевёрнутую последовательность
    """
    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Возвращает комплементарную последовательность.
    """
    complement_map = str.maketrans('ATGCUatgcu', 'TACGAtacga')
    return sequence.translate(complement_map)


def reverse_complement(sequence: str) -> str:
    """
    Возвращает комплементарную, перевёрнутую последовательность.
    """
    return reverse(complement(sequence))


def gc_cont(sequence: str) -> float:
    """
    Вычисляет процент содержания G и C в последовательности.
    """
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100


def DNA_or_RNA(sequence: str) -> str:
    """
    Определяет, является ли последовательность ДНК или РНК.
    """
    if ('U' in sequence.upper()) == ('T' in sequence.upper()):
        return 'Не знаю'
    elif 'U' in sequence.upper():
        return 'РНК'
    elif 'T' in sequence.upper():
        return 'ДНК'
