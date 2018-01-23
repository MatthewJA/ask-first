"""Query images from Faint Images of the Radio Sky at Twenty-Centimeters.

Matthew Alger <matthew.alger@anu.edu.au>
Research School of Astronomy and Astrophysics
The Australian National University
2018
"""

from __future__ import print_function, division

import argparse
import collections
import json
import os
import warnings

import astropy.io.fits
import astropy.wcs
import numpy
import scipy.spatial.distance


def read_paths(first_path):
    """Read the paths of the FIRST images.

    Notes
    -----
    This function is necessary so we can pre-generate a list of centres
    based on the filenames. The fact we are aggregating paths instead of
    centres is just a convenience.

    Parameters
    ----------
    first_path : str
        Path to the top-level FIRST directory.

    Returns
    -------
    dict
        Map from (float, float) centre coordinate tuples (in degrees)
        to str paths to FITS image files centred on those coordinates.
    """
    centre_to_path = collections.defaultdict(list)
    for dirpath, dirnames, filenames in os.walk(first_path):
        for filename in filenames:
            if filename.startswith('.'):
                continue

            if not filename.endswith('.fits'):
                continue

            # Filenames are of the form 09210-07233R.fits.
            # Get the centre by slicing the filename.
            ra, dec = filename[:5], filename[5:11]
            # Convert the centre into a more useful form, i.e. degrees.
            # Into hours...
            ra = int(ra[:2]) + int(ra[2:4]) / 60 + int(ra[4]) / 60 / 60
            # ...Then into degrees.
            ra *= 360 / 24
            # dec goes straight into degrees.
            dec = int(dec[:3]) + int(dec[3:5]) / 60 + int(dec[5]) / 60 / 60
            centre = [ra, dec]

            path = os.path.join(dirpath, filename)
            centre_to_path[tuple(centre)].append(path)

    return centre_to_path


def get_image(coord, width, paths):
    """Get an image from FIRST at a coordinate.

    Notes
    -----
    This function currently does not support stitching multiple images
    together, so the image must be within a single image patch.

    Parameters
    ----------
    coord : (float, float)
        Centre of image (RA, dec).

    width : float
        Width in degrees.

    paths : dict | str
        Map from centre coordinates (in degrees) to FITS filenames OR
        path to the FIRST data.

    Returns
    -------
    numpy.ndarray
    """
    # Handle the dict/str options for the paths argument.
    try:
        centres = list(paths.keys())
    except AttributeError:
        paths = read_paths(paths)
        centres = list(paths.keys())
    dists = scipy.spatial.distance.cdist([coord], centres)
    closest = centres[dists.argmin()]
    assert isinstance(closest, tuple)
    # There are, for some reason, four of these images. I can't find
    # any documentation on what the difference is between them, so I'll just
    # pick one ~arbitrarily~.
    # TODO(MatthewJA): Figure out what to do here properly.
    path = paths[closest][0]
    with astropy.io.fits.open(path) as fits:
        header = fits[0].header
        image = fits[0].data
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            wcs = astropy.wcs.WCS(header).dropaxis(3).dropaxis(2)
        ra, dec = coord
        # RA increases to the left.
        min_ra, max_ra = ra + width / 2, ra - width / 2
        height = numpy.cos(numpy.deg2rad(dec)) * width
        min_dec, max_dec = dec - height / 2, dec + height / 2
        ((min_x, min_y), (max_x, max_y)) = wcs.all_world2pix(
            ((min_ra, min_dec), (max_ra, max_dec)), 1)
        assert min_x < max_x
        assert min_y < max_y
        patch = image[0, 0, int(min_y):int(max_y), int(min_x):int(max_x)]
        return patch


def main(first_path, centres_path=None):
    """Loads paths and caches if a path is provided.

    Parameters
    ----------
    first_path : str
        Path to the top-level FIRST directory.
    """


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('first_path', help='Path to FIRST data')
    args = parser.parse_args()

    paths = read_paths(args.first_path)
    im = get_image((162.5302917, 30.6770889), 3 / 60, paths)
