from collections import Counter, defaultdict
import itertools as it

import baseconv
import regex as re


DNA_BASE4 = baseconv.BaseConverter('ACGT')
COMPLEMENTS = str.maketrans('ACGT', 'TGCA')


def pattern_count(genome, pattern):
    '''
    >>> genome = 'GCGCG'
    >>> pattern = 'GCG'
    >>> pattern_count(genome, pattern)
    2

    '''
    return len(re.findall(pattern, genome, overlapped=True))


def kmer_counts(genome, kmer_length):
    pattern = fr'(.{{{kmer_length}}}).*\1'
    repeaters = re.findall(pattern, genome, overlapped=True)
    kmers = Counter(repeaters)

    counts = defaultdict(list)
    for k, v in kmers.items():
        counts[v+1].append(k)
    return counts


def frequent_kmers(genome, kmer_length, min_freq=None):
    '''
    >>> genome = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    >>> kmer_length = 4
    >>> sorted(frequent_kmers(genome, kmer_length))
    ['CATG', 'GCAT']
    >>> sorted(frequent_kmers(genome, kmer_length, min_freq=2))
    ['ATGA', 'CATG', 'GCAT', 'TGCA']
    >>> frequent_kmers('GCGAT', 3)
    []

    '''
    counts = kmer_counts(genome, kmer_length)

    if not counts:
        return []
    elif min_freq is None:
        return counts[max(counts)]
    else:
        return list(it.chain.from_iterable(
            [counts[i]
             for i in it.takewhile(lambda i: i >= min_freq,
                                   sorted(counts.keys(), reverse=True))]))


def pattern_to_number(pattern):
    '''
    >>> pattern_to_number('GT')
    11
    >>> pattern_to_number('AGT')
    11

    '''
    return int(DNA_BASE4.decode(pattern))


def number_to_pattern(number, n_bases=0):
    '''
    >>> number_to_pattern(11)
    'GT'
    >>> number_to_pattern(11, n_bases=3)
    'AGT'

    '''
    return f'{DNA_BASE4.encode(number):A>{n_bases}}'


def window(genome, n_bases):
    '''
    >>> genome = 'GCGCG'
    >>> list(window(genome, 2))
    ['GC', 'CG', 'GC', 'CG']
    >>> list(window(genome, 3))
    ['GCG', 'CGC', 'GCG']

    '''
    for i in range(len(genome) - n_bases + 1):
        yield genome[i:i + n_bases]


def computing_frequencies(genome, n_bases):
    '''
    >>> genome = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    >>> computing_frequencies(genome, 2)
    [0, 1, 2, 4, 3, 0, 2, 1, 3, 4, 0, 2, 0, 1, 5, 1]

    '''
    frequencies = [0] * (4**n_bases)
    for substring in window(genome, n_bases):
        frequencies[pattern_to_number(substring)] += 1
    return frequencies


def dict_frequencies(genome, n_bases):
    '''
    >>> genome = 'GATTACA'
    >>> sorted(dict(dict_frequencies(genome, 2)).items())
    [('AC', 1), ('AT', 1), ('CA', 1), ('GA', 1), ('TA', 1), ('TT', 1)]

    '''
    frequencies = defaultdict(int)
    for substring in window(genome, n_bases):
        frequencies[substring] += 1
    return frequencies


def reverse_complement(genome):
    '''
    >>> reverse_complement('AAAACCCGGT')
    'ACCGGGTTTT'

    '''
    return genome.translate(COMPLEMENTS)[::-1]


def start_positions(genome, pattern):
    '''
    >>> start_positions('GATATATGCATATACTT', 'ATAT')
    [1, 3, 9]

    '''
    return [m.start() for m in re.finditer(pattern, genome, overlapped=True)]


def clumping_naive(genome, kmer_length, window_size, min_freq):
    '''
    >>> genome = 'GATCAGCATAAGGGTCCCTGCAATGCATGACAAGCCTGCAGTTGTTTTAC'
    >>> clumping_naive(genome, 4, 25, 3)
    {'TGCA'}

    '''
    patterns = set()
    for substring in window(genome, window_size):
        for kmer in frequent_kmers(substring, kmer_length, min_freq):
            patterns.add(kmer)
    return patterns


def clumping(genome, kmer_length, window_size, min_freq):
    '''
    >>> genome = 'GATCAGCATAAGGGTCCCTGCAATGCATGACAAGCCTGCAGTTGTTTTAC'
    >>> clumping(genome, 4, 25, 3)
    {'TGCA'}

    '''
    freqs = dict_frequencies(genome[:window_size], kmer_length)
    patterns = {k for k, v in freqs.items() if v >= min_freq}
    for substring in window(genome, window_size+1):
        first_pattern = substring[:kmer_length]
        last_pattern = substring[-kmer_length:]

        freqs[first_pattern] -= 1
        freqs[last_pattern] += 1

        if freqs[last_pattern] >= min_freq:
            patterns.add(last_pattern)
    return patterns


def skew(genome):
    '''
    >>> list(skew('CATGGGCATCGGCCATACGCC'))
    [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]

    '''
    diff = 0
    yield diff
    for base in genome:
        if base == 'C':
            diff -= 1
        elif base == 'G':
            diff += 1
        yield diff


def min_skew(genome):
    '''
    >>> genome = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
    >>> list(min_skew(genome))
    [11, 24]

    '''
    skew_list = list(skew(genome))
    minimum = min(skew_list)
    yield from (i for i, v in enumerate(skew_list) if v == minimum)


def hamming_distance(s1, s2):
    '''
    >>> s1 = 'GGGCCGTTGGT'
    >>> s2 = 'GGACCGTTGAC'
    >>> hamming_distance(s1, s2)
    3

    '''
    return sum(l != r for l, r in zip(s1, s2))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
