from collections import Counter, defaultdict

import baseconv
import regex as re


DNA_BASE4 = baseconv.BaseConverter('ACGT')
COMPLEMENTS = str.maketrans('ACGT', 'TGCA')


def pattern_count(string, pattern):
    '''
    >>> string = 'GCGCG'
    >>> pattern = 'GCG'
    >>> pattern_count(string, pattern)
    2

    '''
    return len(re.findall(pattern, string, overlapped=True))


def frequent_words(string, kmer_length):
    '''
    >>> string = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    >>> kmer_length = 4
    >>> sorted(frequent_words(string, kmer_length))
    ['CATG', 'GCAT']

    '''
    pattern = fr'(.{{{kmer_length}}}).*\1'
    repeaters = re.findall(pattern, string, overlapped=True)
    kmers = Counter(repeaters)

    counts = defaultdict(list)
    for k, v in kmers.items():
        counts[v].append(k)

    return counts[max(counts)]


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


def window(string, n_bases):
    '''
    >>> string = 'GCGCG'
    >>> list(window(string, 2))
    ['GC', 'CG', 'GC', 'CG']
    >>> list(window(string, 3))
    ['GCG', 'CGC', 'GCG']

    '''
    for i in range(len(string) - n_bases + 1):
        yield string[i:i + n_bases]


def computing_frequencies(string, n_bases):
    '''
    >>> string = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    >>> computing_frequencies(string, 2)
    [0, 1, 2, 4, 3, 0, 2, 1, 3, 4, 0, 2, 0, 1, 5, 1]

    '''
    frequencies = [0] * (4**n_bases)
    for substring in window(string, n_bases):
        frequencies[pattern_to_number(substring)] += 1
    return frequencies


def dict_frequencies(string, n_bases):
    '''
    >>> string = 'GATTACA'
    >>> sorted(dict_frequencies(string, 2).items())
    [('AC', 1), ('AT', 1), ('CA', 1), ('GA', 1), ('TA', 1), ('TT', 1)]

    '''
    frequencies = defaultdict(int)
    for substring in window(string, n_bases):
        frequencies[substring] += 1
    return dict(frequencies)


def reverse_complement(string):
    '''
    >>> reverse_complement('AAAACCCGGT')
    'ACCGGGTTTT'

    '''
    return string.translate(COMPLEMENTS)[::-1]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
