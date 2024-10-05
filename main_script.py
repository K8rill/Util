from scripts.dna_rna_tools import transcribe, reverse, complement, reverse_complement
from scripts.filter_fastq import calculate_gc_content, calculate_average_quality
from typing import List, Union, Dict, Tuple

def run_dna_rna_tools(*args: str) -> Union[str, List[str]]:
    """
    Выполняет указанную операцию над переданными последовательностями ДНК или РНК.
    
    Args:произвольное количество строк с последовательностями, последним аргументом передается название процедуры 
    (transcribe, reverse, complement, reverse_complement, gc_cont, DNA_or_RNA)
    
    Return: если передана одна последовательность — str,
            если несколько последовательностей — list[str, str].
    """
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


def filter_fastq(
    seqs: Dict[str, Tuple[str, str]], 
    gc_bounds: Union[Tuple[float, float], float] = (0, 100), 
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32), 
    quality_threshold: float = 0
) -> Dict[str, Tuple[str, str]]:
    """
    Фильтрация последовательностей по GC-составу, длине и среднему качеству.
    
    Args:
    seqs: словарь последовательностей, где ключ — имя, а значение — кортеж (последовательность, качество)
    gc_bounds: интервал GC-состава в процентах (нижняя и верхняя границы)
    length_bounds: интервал длины последовательностей (нижняя и верхняя границы)
    quality_threshold: пороговое значение среднего качества рида
    
    Return: отфильтрованный словарь сиквенсов
    """
    filtered_seqs = {}

    # Приведение границ, если передано одно число
    if isinstance(gc_bounds, (float, int)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    for name, (sequence, quality) in seqs.items():
        # Фильтр по длине
        seq_length = len(sequence)
        if not (length_bounds[0] <= seq_length <= length_bounds[1]):
            continue
        
        # Фильтр по GC-составу
        gc_content = calculate_gc_content(sequence)
        if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
            continue
        
        # Фильтр по среднему качеству
        avg_quality = calculate_average_quality(quality)
        if avg_quality < quality_threshold:
            continue
        
        # Если все условия соблюдены, добавляем последовательность в словарь
        filtered_seqs[name] = (sequence, quality)
    
    return filtered_seqs


