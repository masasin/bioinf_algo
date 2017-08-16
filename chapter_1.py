from collections import Counter, defaultdict
from functools import lru_cache
import itertools as it

import baseconv
import regex as re


BASES = 'ACGT'
DNA_BASE4 = baseconv.BaseConverter(BASES)
COMPLEMENTS = str.maketrans(BASES, 'TGCA')


def pattern_count(genome, pattern):
    '''
    >>> genome = 'GCGCG'
    >>> pattern = 'GCG'
    >>> pattern_count(genome, pattern)
    2

    '''
    return len(re.findall(pattern, genome, overlapped=True))


def kmer_counts(genome, kmer_length):
    '''
    >>> genome = 'GCGCG'
    >>> sorted(dict(kmer_counts(genome, 2)).items())
    [(2, ['GC', 'CG'])]

    '''
    freqs = dict_frequencies(genome, kmer_length)
    counts = defaultdict(list)
    for k, v in freqs.items():
        counts[v].append(k)
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

    if max(counts) <= 1:
        return []
    elif min_freq is None:
        return counts[max(counts)]
    else:
        return list(it.chain.from_iterable([counts[i]
            for i in it.takewhile(lambda i: i >= min_freq,
                                  sorted(counts, reverse=True))]))


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
    >>> list(start_positions('GATATATGCATATACTT', 'ATAT'))
    [1, 3, 9]

    '''
    yield from (m.start()
                for m in re.finditer(pattern, genome, overlapped=True))


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


@lru_cache()
def hamming_distance(s1, s2):
    '''
    >>> s1 = 'GGGCCGTTGGT'
    >>> s2 = 'GGACCGTTGAC'
    >>> hamming_distance(s1, s2)
    3

    '''
    return sum(l != r for l, r in zip(s1, s2))


def start_positions_approx(genome, pattern, dist_max, method=hamming_distance):
    '''
    >>> genome = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
    >>> pattern = 'ATTCTGGA'
    >>> distance = 3
    >>> list(start_positions_approx(genome, pattern, distance))
    [6, 7, 26, 27]

    '''
    for i, substring in enumerate(window(genome, len(pattern))):
        if method(pattern, substring) <= dist_max:
            yield i


def pattern_count_approx(genome, pattern, dist_max, method=hamming_distance):
    '''
    >>> pattern_count_approx('AACAAGCTGATAAACATTTAAAGAG', 'AAAAA', 1)
    4
    >>> pattern_count_approx('AACAAGCTGATAAACATTTAAAGAG', 'AAAAA', 2)
    11

    '''
    return sum(method(pattern, substring) <= dist_max
               for substring in window(genome, len(pattern)))


def neighbors(pattern, dist_max, method=hamming_distance):
    '''
    >>> sorted(neighbors('ACG', 1))
    ['AAG', 'ACA', 'ACC', 'ACG', 'ACT', 'AGG', 'ATG', 'CCG', 'GCG', 'TCG']

    '''
    if dist_max == 0:
        return set([pattern])
    if len(pattern) == 1:
        return set(BASES)

    neighborhood = set()
    head, tail = pattern[0], pattern[1:]
    neighbors_tail = neighbors(tail, dist_max, method)

    for neighbor in neighbors_tail:
        if method(tail, neighbor) < dist_max:
            neighborhood.update({base + neighbor for base in BASES})
        else:
            neighborhood.add(head + neighbor)
    return neighborhood


def neighbors_iterative(pattern, dist_max, method=hamming_distance):
    '''
    >>> sorted(neighbors_iterative('ACG', 1))
    ['AAG', 'ACA', 'ACC', 'ACG', 'ACT', 'AGG', 'ATG', 'CCG', 'GCG', 'TCG']

    '''
    def immediate_neighbor(pattern):
        return {pattern[:i] + base + pattern[i+1:]
                for i in range(len(pattern))
                for base in BASES}

    neighborhood = set([pattern])
    for dist in range(dist_max):
        for neighbor in neighborhood.copy():
            neighborhood.update(immediate_neighbor(neighbor))
    return neighborhood


if __name__ == '__main__':
    import doctest
    doctest.testmod()
