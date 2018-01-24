"""Tests for ask_first.py.

Matthew Alger <matthew.alger@anu.edu.au>
Research School of Astronomy and Astrophysics
The Australian National University
2018
"""

from __future__ import print_function, division

import os.path
import unittest

import numpy

import ask_first

# Path to directory containing this code.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# Path to directory containing the test data.
DATA_PATH = os.path.join(THIS_DIR, 'data')


class TestGetImage(unittest.TestCase):

    def test(self):
        # These coordinates are the defaults for the FIRST cutout server.
        coords = 162.5302917, 30.6770889
        width = 3 / 60

        im = ask_first.get_image(coords, width, DATA_PATH)

        reference_im = numpy.load(os.path.join(
            DATA_PATH,
            'test_data_162.5302917_30.6770889.npy'))
        numpy.testing.assert_allclose(reference_im, im)


class TestReadCatalogue(unittest.TestCase):

    def test(self):
        cat = ask_first.read_catalogue(os.path.join(
            DATA_PATH,
            'test_catalog_14dec17.bin'))
        self.assertEqual(len(cat), 27)
        self.assertTrue(all(cat['Field'] == '10510+30456E'))


if __name__ == '__main__':
    unittest.main()
