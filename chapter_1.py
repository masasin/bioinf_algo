from collections import Counter, defaultdict
from functools import lru_cache
import itertools as it

import baseconv


BASES = 'ACGT'
DNA_BASE4 = baseconv.BaseConverter(BASES)
COMPLEMENTS = str.maketrans(BASES, 'TGCA')


@lru_cache()
def hamming_distance(s1, s2):
    '''
    >>> s1 = 'GGGCCGTTGGT'
    >>> s2 = 'GGACCGTTGAC'
    >>> hamming_distance(s1, s2)
    3

    '''
    return sum(l != r for l, r in zip(s1, s2))


def neighbors(pattern, dist_max, *, method=hamming_distance):
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
    neighbors_tail = neighbors(tail, dist_max, method=method)

    for neighbor in neighbors_tail:
        if method(tail, neighbor) < dist_max:
            neighborhood.update({base + neighbor for base in BASES})
        else:
            neighborhood.add(head + neighbor)
    return neighborhood


def pattern_count(genome, pattern, *, dist_max=0, method=hamming_distance):
    '''
    >>> pattern_count('GCGCG', 'GCG')
    2
    >>> pattern_count('AACAAGCTGATAAACATTTAAAGAG', 'AAAAA', dist_max=1)
    4
    >>> pattern_count('AACAAGCTGATAAACATTTAAAGAG', 'AAAAA', dist_max=2)
    11

    '''
    return sum(method(pattern, substring) <= dist_max
               for substring in window(genome, len(pattern)))


def kmer_counts(genome, kmer_length, *, dist_max=0, method=hamming_distance, reverse=False):
    '''
    >>> genome = 'GCGCG'
    >>> sorted(dict(kmer_counts(genome, 2)).items())
    [(2, ['GC', 'CG'])]
    >>> counts = kmer_counts(genome, 2, dist_max=1)
    >>> sorted(counts) == [2, 4]
    True
    >>> sorted(sorted(v) for v in counts.values())
    [['AC', 'AG', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TC', 'TG'], ['CC', 'GG']]
    >>> counts = kmer_counts(genome, 2, dist_max=1, reverse=True)
    >>> sorted(counts) == [4, 8]
    True

    '''
    freqs = dict_frequencies(genome, kmer_length, dist_max=dist_max, method=method, reverse=reverse)
    counts = defaultdict(list)
    for k, v in freqs.items():
        counts[v].append(k)
    return counts


def frequent_kmers(genome, kmer_length, *, dist_max=0, min_freq=None, method=hamming_distance, reverse=False):
    '''
    >>> genome = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    >>> kmer_length = 4
    >>> sorted(frequent_kmers(genome, kmer_length))
    ['CATG', 'GCAT']
    >>> sorted(frequent_kmers(genome, kmer_length, min_freq=2))
    ['ATGA', 'CATG', 'GCAT', 'TGCA']
    >>> sorted(frequent_kmers('GCGAT', 3))
    []
    >>> sorted(frequent_kmers(genome, kmer_length, dist_max=1))
    ['ATGC', 'ATGT', 'GATG']
    >>> sorted(frequent_kmers(genome, kmer_length, dist_max=1, min_freq=4))
    ['AATG', 'ACAT', 'AGCA', 'ATGA', 'ATGC', 'ATGG', 'ATGT', 'CATG', 'CATT', 'CGTG', 'CTTG', 'GATG', 'GCGT', 'GGAT', 'TATG', 'TCAT', 'TGAA', 'TGCT']
    >>> sorted(frequent_kmers('GCGAT', 4, dist_max=1))
    []
    >>> sorted(frequent_kmers(genome, kmer_length, dist_max=1, reverse=True))
    ['ACAT', 'ATGT']

    '''
    counts = kmer_counts(genome, kmer_length, dist_max=dist_max, method=method, reverse=reverse)

    if max(counts) <= 1:
        return
    elif min_freq is None:
        yield from counts[max(counts)]
    else:
        yield from it.chain.from_iterable([counts[i]
            for i in it.takewhile(lambda i: i >= min_freq,
                                  sorted(counts, reverse=True))])


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


def list_frequencies(genome, n_bases, *, dist_max=0, method=hamming_distance):
    '''
    >>> genome = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    >>> list_frequencies(genome, 2)
    [0, 1, 2, 4, 3, 0, 2, 1, 3, 4, 0, 2, 0, 1, 5, 1]
    >>> genome = 'AAGCAAAGGTGGG'
    >>> list_frequencies(genome, 2, dist_max=1)
    [6, 6, 9, 6, 4, 2, 7, 2, 9, 5, 8, 5, 5, 2, 6, 2]

    '''
    frequencies = [0] * (4**n_bases)
    for substring in window(genome, n_bases):
        for neighbor in neighbors(substring, dist_max, method=method):
            frequencies[pattern_to_number(neighbor)] += 1
    return frequencies


def dict_frequencies(genome, n_bases, *, dist_max=0, method=hamming_distance, reverse=False):
    '''
    >>> genome = 'GATTACA'
    >>> sorted(dict(dict_frequencies(genome, 2)).items())
    [('AC', 1), ('AT', 1), ('CA', 1), ('GA', 1), ('TA', 1), ('TT', 1)]
    >>> sorted(dict(dict_frequencies(genome, 2, dist_max=1)).items())
    [('AA', 5), ('AC', 2), ('AG', 2), ('AT', 3), ('CA', 3), ('CC', 2), ('CG', 1), ('CT', 3), ('GA', 3), ('GC', 2), ('GG', 1), ('GT', 3), ('TA', 4), ('TC', 3), ('TG', 2), ('TT', 3)]
    >>> sorted(dict(dict_frequencies(genome, 2, dist_max=1, reverse=True)).items())
    [('AA', 8), ('AC', 5), ('AG', 5), ('AT', 6), ('CA', 5), ('CC', 3), ('CG', 2), ('CT', 5), ('GA', 6), ('GC', 4), ('GG', 3), ('GT', 5), ('TA', 8), ('TC', 6), ('TG', 5), ('TT', 8)]

    '''
    frequencies = defaultdict(int)
    for substring in window(genome, n_bases):
        for neighbor in neighbors(substring, dist_max, method=method):
            frequencies[neighbor] += 1
            if reverse:
                frequencies[reverse_complement(neighbor)] += 1
    return frequencies


def reverse_complement(genome):
    '''
    >>> reverse_complement('AAAACCCGGT')
    'ACCGGGTTTT'

    '''
    return genome.translate(COMPLEMENTS)[::-1]


def start_positions(genome, pattern, *, dist_max=0, method=hamming_distance):
    '''
    >>> list(start_positions('GATATATGCATATACTT', 'ATAT'))
    [1, 3, 9]
    >>> genome = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
    >>> list(start_positions(genome, 'ATTCTGGA', dist_max=3))
    [6, 7, 26, 27]

    '''
    for i, substring in enumerate(window(genome, len(pattern))):
        if method(pattern, substring) <= dist_max:
            yield i


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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
