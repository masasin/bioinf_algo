from functools import partial
import math
import random

import chapter_1 as c


class TestDistance:
    def test_hamming_distance(self):
        s1 = 'GGGCCGTTGGT'
        s2 = 'GGACCGTTGAC'
        assert c.hamming_distance(s1, s2) == 3

    def test_neighbors(self):
        assert c.neighbors('ACG', 1) == {
            'ACG', 'CCG', 'GCG', 'TCG',
            'AAG', 'ACG', 'AGG', 'ATG',
            'ACA', 'ACC', 'ACG', 'ACT',
        }


class TestPatternCount:
    def setup(self):
        genome = 'AACAAGCTGATAAACATTTAAAGAG'
        pattern = 'AAAAA'
        self.f = partial(c.pattern_count, genome, pattern)

    def test_simple(self):
        assert c.pattern_count('GCGCG', 'GCG') == 2
        assert self.f() == 0

    def test_mutable(self):
        assert self.f(dist_max=1) == 4
        assert self.f(dist_max=2) == 11


class TestKmerCounts:
    def setup(self):
        self.genome = 'GCGCG'
        self.kmer_length = 2
        self.f = partial(c.kmer_counts, self.genome, self.kmer_length)

    def test_simple(self):
        assert self.f() == {
            2: {'GC', 'CG'},
        }

    def test_mutable(self):
        assert self.f(dist_max=1) == {
            2: {'AC', 'GC', 'TC', 'GT', 'GA', 'TG', 'CG', 'CT', 'CA', 'AG'},
            4: {'GG', 'CC'},
        }

    def test_mutable_reverse(self):
        assert self.f(dist_max=1, reverse=True) == {
            4: {'AC', 'GC', 'TC', 'GT', 'GA', 'TG', 'CG', 'CT', 'CA', 'AG'},
            8: {'GG', 'CC'},
        }


class TestFrequentKmers:
    def setup(self):
        self.genome = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        self.kmer_length = 4
        self.f = partial(c.frequent_kmers, self.genome, self.kmer_length)

    def test_simple(self):
        assert set(self.f()) == {'CATG', 'GCAT'}
        assert set(self.f(min_freq=2)) == {'ATGA', 'CATG', 'GCAT', 'TGCA'}
        assert set(self.f(min_freq=5)) == set()
    
    def test_mutable(self):
        counts = c.kmer_counts(self.genome, self.kmer_length, dist_max = 1)
        maximum = max(counts)
        top_two = counts[maximum] | counts[maximum - 1]

        assert set(self.f(dist_max=1)) == {'ATGC', 'ATGT', 'GATG'}
        assert set(self.f(dist_max=1, min_freq=maximum - 1)) == top_two
        assert set(self.f(dist_max=1, min_freq=maximum + 1)) == set()

    def test_mutable_reverse(self):
        assert set(self.f(dist_max=1, reverse=True)) == {'ACAT', 'ATGT'}


class TestTransformations:
    def test_pattern_to_number(self):
        assert c.pattern_to_number('GT') == 11
        assert c.pattern_to_number('AGT') == 11

    def test_number_to_pattern(self):
        assert c.number_to_pattern(11) == 'GT'
        assert c.number_to_pattern(11, n_bases=3) == 'AGT'

    def test_reverse_complement(self):
        assert c.reverse_complement('AAAACCCGGT') == 'ACCGGGTTTT'

        genome = ''.join(random.choices(c.BASES, k=random.randint(0, 100)))
        assert c.reverse_complement(c.reverse_complement(genome)) == genome


def test_window():
    genome = 'GCGCG'
    assert list(c.window(genome, 2)) == ['GC', 'CG', 'GC', 'CG']
    assert list(c.window(genome, 3)) == ['GCG', 'CGC', 'GCG']


class TestFrequencies:
    def setup(self):
        self.genome = 'AAGCAAAGGTGGG'
        self.n_bases = 2
        self.f = partial(c.frequencies, self.genome, self.n_bases)

    def test_simple(self):
        frequency_dict = self.f()
        assert frequency_dict == {
            'AA': 3, 'AG': 2,
            'CA': 1,
            'GC': 1, 'GG': 3, 'GT': 1,
            'TG': 1,
        }
        assert c.vectorize(frequency_dict) == [
            3, 0, 2, 0,
            1, 0, 0, 0,
            0, 1, 3, 1,
            0, 0, 1, 0,
        ]

    def test_mutable(self):
        frequency_dict = self.f(dist_max=1)
        assert frequency_dict == {
            'AA': 6, 'AC': 6, 'AG': 9, 'AT': 6,
            'CA': 4, 'CC': 2, 'CG': 7, 'CT': 2,
            'GA': 9, 'GC': 5, 'GG': 8, 'GT': 5,
            'TA': 5, 'TC': 2, 'TG': 6, 'TT': 2,
        }
        assert c.vectorize(frequency_dict) == [
            6, 6, 9, 6,
            4, 2, 7, 2,
            9, 5, 8, 5,
            5, 2, 6, 2,
        ]

    def test_mutable_reverse(self):
        frequency_dict = self.f(dist_max=1, reverse=True)
        assert frequency_dict == {
            'AA': 8, 'AC': 11, 'AG': 11, 'AT': 12,
            'CA': 10, 'CC': 10, 'CG': 14, 'CT': 11,
            'GA': 11, 'GC': 10, 'GG': 10, 'GT': 11,
            'TA': 10, 'TC': 11, 'TG': 10, 'TT': 8,
        }
        assert c.vectorize(frequency_dict) == [
            8, 11, 11, 12,
            10, 10, 14, 11,
            11, 10, 10, 11,
            10, 11, 10, 8,
        ]


class TestSkew:
    def test_skew(self):
        genome = 'CATGGGCATCGGCCATACGCC'
        assert list(c.skew(genome)) == [
                   0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0,
                   1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2
                ]

    def test_min_skew(self):
        genome = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
        assert list(c.min_skew(genome)) == [11, 24]


class TestStartPositions:
    def test_simple(self):
        assert list(c.start_positions('GATATATGCATATACTT', 'ATAT')) == [1, 3, 9]

    def mutable(self):
        genome = ('CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGT'
                  'CAATCAAAT')
        assert (list(c.start_positions(genome, 'ATTCTGGA', dist_max=3))
                == [6, 7, 26, 27])


class TestClumpers:
    def setup(self):
        self.genome = 'GATCAGCATAAGGGTCCCTGCAATGCATGACAAGCCTGCAGTTGTTTTAC'
        self.kmer_length = 4
        self.f = partial(c.clumpers, self.genome, self.kmer_length)

    def test_full_window(self):
        assert self.f(min_freq=3) == {'TGCA'}

    def test_windowed(self):
        assert self.f(window_size=25, min_freq=3) == {'TGCA'}


class TestProbabilities:
    def test_prob_pattern(self):
        f = partial(c.prob_pattern, length=4, bases='01')
        assert f('01', exact=True) == 11/16
        assert f('11', exact=True) == 8/16
        assert f('01', min_freq=2, exact=True) == 1/16
        assert f('11', min_freq=2, exact=True) == 3/16
        assert f('11', min_freq=2) == 1/16
        assert math.isclose(c.prob_pattern('ACTAT', 30, min_freq=3), 7.599e-7,
                            rel_tol=1e-5, abs_tol=1e-8)

    def test_prob_kmer(self):
        assert c.prob_kmer(9, 500, min_freq=3) == 17861900 / 68719476736


class TestOverlapCorrelation:
    def setup(self):
        self.s1 = '011011'
        self.s2 = '110110'

    def test_overlap_corr(self):
        assert (list(c.overlap_corr(self.s1, self.s2))
                == [False, True, False, False, True, True])


    def test_overlap_corr_freq(self):
        assert c.overlap_corr_freq(self.s1, self.s2) == 19/32
