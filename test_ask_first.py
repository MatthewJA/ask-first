"""Tests for ask_first.py.

Matthew Alger <matthew.alger@anu.edu.au>
Research School of Astronomy and Astrophysics
The Australian National University
2018
"""

from __future__ import print_function, division

import unittest

import numpy

import ask_first


# Path for the FIRST data. Since this library is all about querying this
# data, it's a dependency. Change this path before running tests.
FIRST_PATH = '/Volumes/Alger/archive.stsci.edu/pub/vla_first/data/'


def fun(x):
    return x + 1

class TestGetImage(unittest.TestCase):
    def test(self):
        self.assertEqual(fun(3), 4)
        # These coordinates are the defaults for the FIRST cutout server.
        coords = 162.5302917, 30.6770889
        width = 3 / 60
        im = ask_first.get_image(coords, width, FIRST_PATH)
        reference_im = numpy.load('test_data_162.5302917_30.6770889.npy')
        numpy.testing.assert_allclose(reference_im, im)


if __name__ == '__main__':
    unittest.main()
