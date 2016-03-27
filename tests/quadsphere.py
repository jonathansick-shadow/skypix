#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

from itertools import izip
import math
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.skypix as skypix


class QuadSpherePixelizationTestCase(unittest.TestCase):
    """Tests for quad-sphere based sky pixelization.
    """

    def testIdAndCoord(self):
        for R in (3, 4, 16, 17):
            qs = skypix.QuadSpherePixelization(R, 0.0)
            self.assertEqual(len(qs), 6 * R**2)
            for i in xrange(6 * R**2):
                root, x, y = qs.coords(i)
                self.assertEqual(qs.id(root, x, y), i)

    def _checkNeighbors(self, qs, root, x, y, tolerance):
        i = qs.id(root, x, y)
        neighbors = qs.getNeighbors(i)
        centers = [qs.getCenter(n) for n in neighbors]
        c = qs.getCenter(i)
        for v in centers:
            self.assertTrue(geom.cartesianAngularSep(c, v) < tolerance)

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
            qs = skypix.QuadSpherePixelization(R, 0.0)
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
            qs = skypix.QuadSpherePixelization(R, 0.0)
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
            qs = skypix.QuadSpherePixelization(R, 0.0)
            lag = 8 * R
            for i in qs:
                neighbors = qs.getNeighbors(i)
                maxN = max(neighbors)
                minN = min(neighbors)
                self.assertTrue(maxN - i < lag)
                self.assertTrue(i - minN < lag)

    def testIntersect(self):
        """Tests polygon sky-pixel intersection.
        """
        qs1 = skypix.QuadSpherePixelization(3, 0.0)
        polygons = [qs1.getGeometry(qs1.id(r, 1, 1), True) for r in xrange(6)]
        qs2 = skypix.QuadSpherePixelization(4, 0.0)
        results = [set([qs2.id(r, 1, 1), qs2.id(r, 2, 1),
                        qs2.id(r, 1, 2), qs2.id(r, 2, 2)]) for r in xrange(6)]
        for poly, res in izip(polygons, results):
            self.assertEqual(set(qs2.intersect(poly)), res)

    def testGeometry(self):
        """Verifies that pixels are contained in their padded selves.
        """
        for R in (4, 5):
            qs = skypix.QuadSpherePixelization(R, math.radians(1.0))
            for i in xrange(6 * R ** 2):
                pixel = qs.getGeometry(i, True)
                paddedPixel = qs.getGeometry(i, False)
                self.assertTrue(paddedPixel.contains(pixel))

    def testPixel(self):
        """Verifies that the center of pixel P is mapped back to P.
        """
        for R in (64, 65):
            qs = skypix.QuadSpherePixelization(R, 0.0)
            for i in qs:
                c = geom.sphericalCoords(qs.getCenter(i))
                pixelId = qs.pixel(math.radians(c[0]), math.radians(c[1]))
                self.assertEqual(i, pixelId)

    def testImageToPixels(self):
        """Tests intersection of an image (WCS and dimensions) with a
        quad-sphere pixelization.
        """
        # metadata taken from CFHT data v695856-e0/v695856-e0-c000-a00.sci_img.fits
        metadata = dafBase.PropertySet()
        metadata.set("SIMPLE", "T")
        metadata.set("BITPIX", -32)
        metadata.set("NAXIS", 2)
        metadata.set("NAXIS1", 1024)
        metadata.set("NAXIS2", 1153)
        metadata.set("RADECSYS", "FK5")
        metadata.set("EQUINOX", 2000.0)
        metadata.set("CRVAL1", 215.604025685476)
        metadata.set("CRVAL2", 53.1595451514076)
        metadata.set("CRPIX1", 1109.99981456774)
        metadata.set("CRPIX2", 560.018167811613)
        metadata.set("CTYPE1", "RA---TAN")
        metadata.set("CTYPE2", "DEC--TAN")
        metadata.set("CD1_1", 5.10808596133527E-05)
        metadata.set("CD1_2", 1.85579539217196E-07)
        metadata.set("CD2_2", -5.10281493481982E-05)
        metadata.set("CD2_1", -8.27440751733828E-07)
        wcs = afwImage.makeWcs(metadata)
        qs = skypix.createQuadSpherePixelization()
        poly = skypix.imageToPolygon(wcs, 1024, 1153)
        pixels = qs.intersect(poly)
        self.assertEqual(len(pixels), 1)
        self.assertEqual(pixels[0], 182720)
        self.assertTrue(qs.getGeometry(182720).contains(poly))


def suite():
    utilsTests.init()
    return unittest.TestSuite(unittest.makeSuite(QuadSpherePixelizationTestCase))


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)

