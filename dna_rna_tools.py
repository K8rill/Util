def run_dna_rna_tools(*args):
    dl = len(args)
    Fun = args[dl - 1]
    res = []
    for i in range(dl - 1):
        if Fun == 'transcribe':
            res.append(transcribe(args[i]))
        if Fun == 'reverse':
            res.append(reverse(args[i]))
        if Fun == 'complement':
            res.append(complement(args[i]))
        if Fun == 'reverse_complement':
            res.append(reverse_complement(args[i]))
        if Fun == 'gc_cont':
            res.append(gc_cont(args[i]))
        if Fun == 'DNA_or_RNA':
            res.append(DNA_or_RNA(args[i]))
    if len(res) == 1:
        return res[0]
    else:
        return res


def transcribe(args):
    return args.replace('T', 'U').replace('t', 'u')


def reverse(args):
    return args[::-1]


def complement(args):
    complement_map = str.maketrans('ATGCUatgcu', 'TACGAtacga')
    return args.translate(complement_map)


def reverse_complement(args):
    return reverse(complement(args))


def gc_cont(args):
    """Вычисление процента содержания G и C в последовательности."""
    gc_count = args.upper().count('G') + args.upper().count('C')
    return (gc_count / len(args)) * 100


def DNA_or_RNA(args):
    """Проверка типа последовательности."""
    if ('U' in args.upper()) == ('T' in args.upper()):
        return 'Не знаю'
    elif 'U' in args.upper():
        return 'РНК'
    elif 'T' in args.upper():
        return 'ДНК'


if __name__ == "__main__":
    run_dna_rna_tools()
