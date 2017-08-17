from collections import Counter, defaultdict
from functools import lru_cache
import itertools as it

import baseconv


BASES = 'ACGT'
DNA_BASE4 = baseconv.BaseConverter(BASES)
COMPLEMENTS = str.maketrans(BASES, 'TGCA')


@lru_cache()
def hamming_distance(s1, s2):
    return sum(l != r for l, r in zip(s1, s2))


def neighbors(pattern, dist_max, *, method=hamming_distance):
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
    return sum(method(pattern, substring) <= dist_max
               for substring in window(genome, len(pattern)))


def kmer_counts(genome, kmer_length, *, dist_max=0, method=hamming_distance,
                reverse=False):
    freqs = frequencies(genome, kmer_length,
                        dist_max=dist_max, method=method, reverse=reverse)
    counts = defaultdict(set)
    for k, v in freqs.items():
        counts[v].add(k)
    return counts


def frequent_kmers(genome, kmer_length, *, dist_max=0, method=hamming_distance,
                   reverse=False, min_freq=None):
    counts = kmer_counts(genome, kmer_length, dist_max=dist_max, method=method,
                         reverse=reverse)

    if max(counts) <= 1:
        return
    elif min_freq is None:
        yield from counts[max(counts)]
    else:
        yield from it.chain.from_iterable([counts[i]
            for i in it.takewhile(lambda i: i >= min_freq,
                                  sorted(counts, reverse=True))])


def pattern_to_number(pattern):
    return int(DNA_BASE4.decode(pattern))


def number_to_pattern(number, *, n_bases=0):
    return f'{DNA_BASE4.encode(number):A>{n_bases}}'


def window(genome, n_bases):
    for i in range(len(genome) - n_bases + 1):
        yield genome[i:i + n_bases]


def vectorize(frequency_dict):
    n_bases = len(list(frequency_dict.keys())[0])
    return [frequency_dict[number_to_pattern(i, n_bases=n_bases)]
            for i in range(len(BASES)**n_bases)]


def frequencies(genome, n_bases, *, dist_max=0, method=hamming_distance,
                reverse=False):
    freqs = defaultdict(int)
    for substring in window(genome, n_bases):
        for neighbor in neighbors(substring, dist_max, method=method):
            freqs[neighbor] += 1
            if reverse:
                freqs[reverse_complement(neighbor)] += 1
    return freqs


def reverse_complement(genome):
    return genome.translate(COMPLEMENTS)[::-1]


def start_positions(genome, pattern, *, dist_max=0, method=hamming_distance):
    for i, substring in enumerate(window(genome, len(pattern))):
        if method(pattern, substring) <= dist_max:
            yield i


def clumpers(genome, kmer_length, *, window_size=None, min_freq=0):
    if window_size is None:
        window_size = len(genome)

    freqs = frequencies(genome[:window_size], kmer_length)

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
    diff = 0
    yield diff
    for base in genome:
        diff += -1 if base == 'C' else 1 if base == 'G' else 0
        yield diff


def min_skew(genome):
    skew_list = list(skew(genome))
    minimum = min(skew_list)
    yield from (i for i, v in enumerate(skew_list) if v == minimum)
