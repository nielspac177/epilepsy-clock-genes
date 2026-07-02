import math

from epicirc.mr.region_check import _pearson


def test_pearson_perfect_positive():
    assert _pearson([1, 2, 3, 4], [2, 4, 6, 8]) == 1.0


def test_pearson_perfect_negative():
    assert _pearson([1, 2, 3, 4], [4, 3, 2, 1]) == -1.0


def test_pearson_too_few():
    assert math.isnan(_pearson([1, 2], [1, 2]))
