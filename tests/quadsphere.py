from itertools import izip
import math
import pdb
import sys
import unittest

import lsst.utils.tests as utilsTests
import lsst.afw.image as afwImage
import lsst.geom.geometry as g
import lsst.skypix.quadsphere as q


class QuadSpherePixelizationTestCase(unittest.TestCase):
    """Tests for quad-sphere based sky pixelization.
    """
    def testIdAndCoord(self):
        for R in (3,4,16,17):
            qs = q.QuadSpherePixelization(R, 0.0)
            for i in xrange(6 * R * R):
                root, x, y = qs.coords(i)
                self.assertEqual(qs.id(root, x, y), i)

    def _checkNeighbors(self, qs, root, x, y, tolerance):
        i = qs.id(root, x, y)
        neighbors = qs.getNeighbors(i)
        centers = [qs.getCenter(n) for n in neighbors]
        c = qs.getCenter(i)
        for v in centers:
            self.assertTrue(g.cartesianAngularSep(c, v) < tolerance)

    def _checkNeighborCommutativity(self, qs, root, x, y):
        i = qs.id(root, x, y)
        neighbors = qs.getNeighbors(i)
        for n in neighbors:
            self.assertTrue(i in qs.getNeighbors(n))

    def testNeighborProximity(self):
        """Tests that neighbors of a pixel P are close to P.
        """
        tolerance = 1.0
        for R in (180, 181):
            qs = q.QuadSpherePixelization(R, 0.0)
            for root in xrange(6):
                # test root pixel edges and corners
                for i in xrange(R - 1):
                    self._checkNeighbors(qs, root, i, 0, tolerance)
                    self._checkNeighbors(qs, root, 0, i, tolerance)
                    self._checkNeighbors(qs, root, i, R - 1, tolerance)
                    self._checkNeighbors(qs, root, R - 1, i, tolerance)
                # test root pixel centers
                self._checkNeighbors(qs, root, R / 2, R / 2, tolerance)

    def testNeighborCommutativity(self):
        """Tests that if pixel P1 is a neighbor of pixel P2,
        then P2 is a neighbor of P1.
        """
        for R in (16, 17):
            qs = q.QuadSpherePixelization(R, 0.0)
            for root in xrange(6):
                # test root pixel edges and corners
                for i in xrange(R - 1):
                    self._checkNeighborCommutativity(qs, root, i, 0)
                    self._checkNeighborCommutativity(qs, root, 0, i)
                    self._checkNeighborCommutativity(qs, root, i, R - 1)
                    self._checkNeighborCommutativity(qs, root, R - 1, i)
                # test root pixel centers
                self._checkNeighborCommutativity(qs, root, R / 2, R / 2)

    def testSpiral(self):
        """Tests that pixels are ordered in an approximate spiral.
        """
        for R in (17, 32):
            qs = q.QuadSpherePixelization(R, 0.0)
            lag = 8 * R
            for i in xrange(6 * R * R):
                neighbors = qs.getNeighbors(i)
                maxN = max(neighbors)
                minN = min(neighbors)
                self.assertTrue(maxN - i < lag)
                self.assertTrue(i - minN < lag)

    def testIntersect(self):
        """Tests polygon sky-pixel intersection.
        """
        qs1 = q.QuadSpherePixelization(3, 0.0)
        polygons = [qs1.getGeometry(qs1.id(r, 1, 1), True) for r in xrange(6)]
        qs2 = q.QuadSpherePixelization(4, 0.0)
        results = [set([qs2.id(r, 1, 1), qs2.id(r, 2, 1),
                        qs2.id(r, 1, 2), qs2.id(r, 2, 2)]) for r in xrange(6)]
        for poly, res in izip(polygons, results):
            self.assertEqual(set(qs2.intersect(poly)), res)

    def testGeometry(self):
        """Verifies that pixels are contained in their padded selves.
        """
        for R in (4, 5):
            qs = q.QuadSpherePixelization(R, 1.0)
            for i in xrange(6 * R ** 2):
                pixel = qs.getGeometry(i, True)
                paddedPixel = qs.getGeometry(i, False)
                self.assertTrue(paddedPixel.contains(pixel))

def suite():
    utilsTests.init()
    return unittest.TestSuite(unittest.makeSuite(QuadSpherePixelizationTestCase))

def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run()
