"""Query images from Faint Images of the Radio Sky at Twenty-Centimeters.

Matthew Alger <matthew.alger@anu.edu.au>
Research School of Astronomy and Astrophysics
The Australian National University
2018
"""

from __future__ import print_function, division

import argparse
import json
import os


def read_paths(centres_path, first_path):
    """Read the paths of the FIRST images into a file.

    Notes
    -----
    This function is necessary so we can pre-generate a list of centres
    based on the filenames. The fact we are aggregating paths instead of
    centres is just a convenience.

    Parameters
    ----------
    centres_path : str
        Path to the output centres file.

    first_path : str
        Path to the top-level FIRST directory.
    """
    centre_to_path = {}
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
            centre = (ra, dec)

            path = os.path.join(dirpath, filename)
            # Filter out the FIRST path (since this is common - saves space).
            path = os.path.relpath(path, first_path)
            # Convert key to str for JSON dump.
            centre_to_path[str(centre)] = path

    with open(centres_path, 'w') as output_file:
        json.dump({'first_path': first_path, 'paths': centre_to_path},
                  output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('first_path', help='Path to FIRST data')
    parser.add_argument('centres_path', help='Path to store centres')
    args = parser.parse_args()
    read_paths(args.centres_path, args.first_path)
