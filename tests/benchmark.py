"""Benchmark ask_first.

Matthew Alger <matthew.alger@anu.edu.au>
Research School of Astronomy and Astrophysics
The Australian National University
2018
"""

from __future__ import print_function, division

import logging
import os.path
import random
import statistics
import string
import time

import ask_first

# Path to directory containing this code.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# Path to directory containing the benchmark data.
DATA_PATH = os.path.join(THIS_DIR, 'data')

random.seed(0)
logger = logging.getLogger(__name__)

def main(step=2):
    assert step in {0, 1, 2}
    ra_step = [80, 100, 200][step]  # Must be a factor of 6400.

    # We want to pretend to have much more data than we really have for
    # testing and benchmarking because we expect the query time to be
    # non-constant. So we will generate a paths dict mapping centres to paths,
    # including paths that don't exist, and ensure that the test data we _do_
    # have is not overlapped by any fake paths.
    # We will use a grid like the real data. Our real data for testing
    # contains just one centre, 10510+30456, so we will use require that this
    # is in the grid.
    # All of the multiplying by 10 keeps everything as an integer.
    # All values are in arcmin.
    ras = [ra for ra in range(110, 13800, ra_step)]
    decs = [dec for dec in range(-53344, 53400, 200)]
    assert 10 * 60 * 10 + 510 in ras
    assert 30 * 60 * 10 + 456 in decs
    centres = [((ra / 10 / 60 / 24 * 360, dec / 10 / 60),
                '{}{}{}{}{}'.format(str(ra // 600).zfill(2),
                                    str(ra % 600).zfill(3),
                                    '-' if dec < 0 else '+',
                                    str(abs(dec) // 600).zfill(2),
                                    str(abs(dec) % 600).zfill(3)))
               for ra in ras for dec in decs]
    assert ((162.75, 30.759999999999998), '10510+30456') in centres
    logger.debug('Generated %d centres.', len(centres))
    paths = {}
    for centre, centre_str in centres:
        # This fakes the epoch multiplicity.
        epochs = [random.choice(string.ascii_uppercase) for _ in range(4)]
        path = [os.path.join(DATA_PATH, centre_str[:5],
                             '{}{}.fits'.format(centre_str, epoch))
                for epoch in epochs]
        paths[centre] = path
    logger.debug('Generated paths.')
    # Now ensure that the real path is included.
    paths[162.75, 30.759999999999998] = [
        os.path.join(DATA_PATH, '10510', p) for p in [
            '10510+30456E.fits',
            '10510+30456F.fits',
            '10510+30456S.fits',
            '10510+30456T.fits']]
    # Benchmarking: Query the test image 1000 times.
    logger.info('Beginning benchmarking with %d centres.', len(centres))
    coords = 162.5302917, 30.6770889
    width = 3 / 60
    times = []
    for i in range(1000):
        t = time.time()
        im = ask_first.get_image(coords, width, paths)
        times.append(time.time() - t)
    logger.info('Done benchmarking.')
    print('Average time: {} +- {} seconds'.format(
        statistics.mean(times), statistics.stdev(times)))


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s-%(name)s-%(levelname)s: %(message)s')
    main()
