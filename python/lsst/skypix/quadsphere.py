import math
import pdb

import lsst.afw.coord.coordLib as coord
import lsst.geom.geometry as g


# Methods for rotating vectors around the x, y, z and -x, -y, -z axes
def _rotX(v, sin, cos):
    return (v[0], v[1] * cos - v[2] * sin, v[1] * sin + v[2] * cos)
def _rotY(v, sin, cos):
    return (v[0] * cos + v[2] * sin, v[1], -v[0] * sin + v[2] * cos)
def _rotZ(v, sin, cos):
    return (v[0] * cos - v[1] * sin, v[0] * sin + v[1] * cos, v[2])
def _rotNX(v, sin, cos):
    return _rotX(v, -sin, cos)
def _rotNY(v, sin, cos):
    return _rotY(v, -sin, cos)
def _rotNZ(v, sin, cos):
    return _rotZ(v, -sin, cos)


class QuadSpherePixelization(object):
    """A quad-sphere based pixelization of the unit sphere. Each of the 6
    root pixels are divided into an R by R grid of sky-pixels, where R
    is the resolution (a positive integer). Sky-pixel identifiers are
    contiguous integers in range [0, 6 * R**2), and are ordered such that
    enumerating the pixels 0, 1, ..., 6 * R**2 results in a roughly spiral
    traversal (from south to north pole) of the unit sphere. In particular,
    iterating over pixels in this order guarantees that when arriving at
    pixel P:

      - all pixels with id P' < P have been visited
      - all neighbors of pixels with id P' < P - 8R have been visited

    Each pixel may optionally be expanded by a padding radius A, such that
    points on the unit sphere outside of a pixel but within angular separation
    A of the fiducial pixel boundaries are also considered to be part of that
    pixel. Pixel (and padded pixel) edges are great circles on the unit
    sphere.

    This class provides methods for computing fiducial or padded geometric
    representations of a sky-pixel - spherical convex polygons in both
    cases. It also supports finding the list of sky-pixels intersecting
    an input spherical convex polygon, as well as determining the neighbors
    and center of a sky-pixel.

    TODO: Right now, pixels are defined by a simple tangent plane projection
    of each cube face onto the sphere, but pixel boundaries could be adjusted
    to produce equal area pixels.
    """
    def __init__(self, resolution, padding):
        """Creates a new quad-sphere sky pixelisation.

        resolution: the number of pixels along the x and y axes of each root
                    pixel.

        padding:    the angular separation (deg) by which fiducial sky-pixels
                    are to be padded.
        """
        if not isinstance(resolution, (int, long)):
            raise TypeError(
                'Quad-sphere resolution parameter must be an int or long')
        if resolution < 3 or resolution > 16384:
            raise RuntimeError(
                'Quad-sphere resolution must be in range [3, 16384]')
        if not isinstance(padding, float):
            raise TypeError(
                'Quad-sphere pixel padding radius must be a float')
        if padding < 0.0 or padding >= 45.0:
            raise RuntimeError(
                'Quad-sphere pixel padding radius must be in range [0, 45) deg')
        R = resolution
        self.resolution = R
        self.padding = padding
        x, y, z = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
        nx, ny, nz = (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0)
        # center of each root pixel
        self.center = [z, x, y, nx, ny, nz]
        # x basis vector of each root pixel
        self.x = [ny, y, nx, ny, x, ny]
        # y basis vector of each root pixel
        self.y = [x, z, z, z, z, nx]
        # functions for rotating vectors in direction of increasing x (per root pixel)
        self.xrot = [_rotX, _rotZ, _rotZ, _rotZ, _rotZ, _rotNX]
        # functions for rotating vectors in direction of increasing y (per root pixel)
        self.yrot = [_rotY, _rotNY, _rotX, _rotY, _rotNX, _rotY]
        # Compute fiducial and splitting x/y planes for each root pixel
        self.xplane = []
        self.yplane = []
        sp = math.sin(math.radians(self.padding))
        for root in xrange(6):
            xplanes = []
            yplanes = []
            c, x, y = self.center[root], self.x[root], self.y[root]
            for i in xrange(R + 1):
                xfp = self._fiducialXPlane(root, i)
                yfp = self._fiducialYPlane(root, i)
                f = 2.0 * float(i) / float(R) - 1.0
                thetaX = math.radians(0.5 * g.cartesianAngularSep(
                    (c[0] + f * x[0] + y[0],
                     c[1] + f * x[1] + y[1],
                     c[2] + f * x[2] + y[2]),
                    (c[0] + f * x[0] - y[0],
                     c[1] + f * x[1] - y[1],
                     c[2] + f * x[2] - y[2])))
                thetaY = math.radians(0.5 * g.cartesianAngularSep(
                    (c[0] + x[0] + f * y[0],
                     c[1] + x[1] + f * y[1],
                     c[2] + x[2] + f * y[2]),
                    (c[0] - x[0] + f * y[0],
                     c[1] - x[1] + f * y[1],
                     c[2] - x[2] + f * y[2])))
                sinrx = sp / math.cos(thetaX)
                sinry = sp / math.cos(thetaY)
                cosrx = math.sqrt(1.0 - sinrx * sinrx)
                cosry = math.sqrt(1.0 - sinry * sinry)
                if i == 0:
                    xlp = None
                    ybp = None
                else:
                    xlp = self.xrot[root](xfp, sinrx, cosrx)
                    ybp = self.yrot[root](yfp, sinry, cosry)
                    xlp = (-xlp[0], -xlp[1], -xlp[2])
                    ybp = (-ybp[0], -ybp[1], -ybp[2])
                if i == R:
                    xrp = None
                    ytp = None
                else:
                    xrp = self.xrot[root](xfp, -sinrx, cosrx)
                    ytp = self.yrot[root](yfp, -sinry, cosry)
                xplanes.append((xfp, xlp, xrp))
                yplanes.append((yfp, ybp, ytp))
            self.xplane.append(xplanes)
            self.yplane.append(yplanes)
        # Corner pixel neighbors.
        self.cornerNeighbors = (
            # root pixel 0
            (
                # x, y = 0, 0
                (self.id(0, 1, 0),
                 self.id(0, 1, 1),
                 self.id(0, 0, 1),
                 self.id(2, R - 1, R - 1),
                 self.id(2, R - 2, R - 1),
                 self.id(3, 0, R - 1),
                 self.id(3, 1, R - 1)),
                # x, y = 0, R-1
                (self.id(0, 1, R - 1),
                 self.id(0, 1, R - 2),
                 self.id(0, 0, R - 2),
                 self.id(1, R - 1, R - 1),
                 self.id(1, R - 2, R - 1),
                 self.id(2, 0, R - 1),
                 self.id(2, 1, R - 1)),
                # x, y = R-1, 0
                (self.id(0, R - 2, 0),
                 self.id(0, R - 2, 1),
                 self.id(0, R - 1, 1),
                 self.id(3, R - 2, R - 1),
                 self.id(3, R - 1, R - 1),
                 self.id(4, 0, R - 1),
                 self.id(4, 1, R - 1)),
                # x, y = R-1, R-1
                (self.id(0, R - 2, R - 1),
                 self.id(0, R - 2, R - 2),
                 self.id(0, R - 1, R - 2),
                 self.id(1, 0, R - 1),
                 self.id(1, 1, R - 1),
                 self.id(4, R - 2, R - 1),
                 self.id(4, R - 1, R - 1))
            ),
            # root pixel 1
            (
                # x, y = 0, 0
                (self.id(1, 1, 0),
                 self.id(1, 1, 1),
                 self.id(1, 0, 1),
                 self.id(4, R - 1, 0),
                 self.id(4, R - 1, 1),
                 self.id(5, R - 2, 0),
                 self.id(5, R - 1, 0)),
                # x, y = 0, R-1
                (self.id(1, 1, R - 1),
                 self.id(1, 1, R - 2),
                 self.id(1, 0, R - 2),
                 self.id(4, R - 1, R - 2),
                 self.id(4, R - 1, R - 1),
                 self.id(0, R - 2, R - 1),
                 self.id(0, R - 1, R - 1)),
                # x, y = R-1, 0
                (self.id(1, R - 2, 0),
                 self.id(1, R - 2, 1),
                 self.id(1, R - 1, 1),
                 self.id(2, 0, 0),
                 self.id(2, 0, 1),
                 self.id(5, 0, 0),
                 self.id(5, 1, 0)),
                # x, y = R-1, R-1
                (self.id(1, R - 2, R - 1),
                 self.id(1, R - 2, R - 2),
                 self.id(1, R - 1, R - 2),
                 self.id(2, 0, R - 2),
                 self.id(2, 0, R - 1),
                 self.id(0, 0, R - 1),
                 self.id(0, 1, R - 1))
            ),
            # root pixel 2
            (
                # x, y = 0, 0
                (self.id(2, 1, 0),
                 self.id(2, 1, 1),
                 self.id(2, 0, 1),
                 self.id(1, R - 1, 0),
                 self.id(1, R - 1, 1),
                 self.id(5, 0, 0),
                 self.id(5, 0, 1)),
                # x, y = 0, R-1
                (self.id(2, 1, R - 1),
                 self.id(2, 1, R - 2),
                 self.id(2, 0, R - 2),
                 self.id(1, R - 1, R - 2),
                 self.id(1, R - 1, R - 1),
                 self.id(0, 0, R - 2),
                 self.id(0, 0, R - 1)),
                # x, y = R-1, 0
                (self.id(2, R - 2, 0),
                 self.id(2, R - 2, 1),
                 self.id(2, R - 1, 1),
                 self.id(3, 0, 0),
                 self.id(3, 0, 1),
                 self.id(5, 0, R - 2),
                 self.id(5, 0, R - 1)),
                # x, y = R-1, R-1
                (self.id(2, R - 2, R - 1),
                 self.id(2, R - 2, R - 2),
                 self.id(2, R - 1, R - 2),
                 self.id(3, 0, R - 2),
                 self.id(3, 0, R - 1),
                 self.id(0, 0, 0),
                 self.id(0, 0, 1))
            ),
            # root pixel 3
            (
                # x, y = 0, 0
                (self.id(3, 1, 0),
                 self.id(3, 1, 1),
                 self.id(3, 0, 1),
                 self.id(2, R - 1, 0),
                 self.id(2, R - 1, 1),
                 self.id(5, 0, R - 1),
                 self.id(5, 1, R - 1)),
                # x, y = 0, R-1
                (self.id(3, 1, R - 1),
                 self.id(3, 1, R - 2),
                 self.id(3, 0, R - 2),
                 self.id(2, R - 1, R - 2),
                 self.id(2, R - 1, R - 1),
                 self.id(0, 0, 0),
                 self.id(0, 1, 0)),
                # x, y = R-1, 0
                (self.id(3, R - 2, 0),
                 self.id(3, R - 2, 1),
                 self.id(3, R - 1, 1),
                 self.id(4, 0, 0),
                 self.id(4, 0, 1),
                 self.id(5, R - 2, R - 1),
                 self.id(5, R - 1, R - 1)),
                # x, y = R-1, R-1
                (self.id(3, R - 2, R - 1),
                 self.id(3, R - 2, R - 2),
                 self.id(3, R - 1, R - 2),
                 self.id(4, 0, R - 2),
                 self.id(4, 0, R - 1),
                 self.id(0, R - 2, 0),
                 self.id(0, R - 1, 0))
            ),
            # root pixel 4
            (
                # x, y = 0, 0
                (self.id(4, 1, 0),
                 self.id(4, 1, 1),
                 self.id(4, 0, 1),
                 self.id(3, R - 1, 0),
                 self.id(3, R - 1, 1),
                 self.id(5, R - 1, R - 2),
                 self.id(5, R - 1, R - 1)),
                # x, y = 0, R-1
                (self.id(4, 1, R - 1),
                 self.id(4, 1, R - 2),
                 self.id(4, 0, R - 2),
                 self.id(3, R - 1, R - 2),
                 self.id(3, R - 1, R - 1),
                 self.id(0, R - 1, 0),
                 self.id(0, R - 1, 1)),
                # x, y = R-1, 0
                (self.id(4, R - 2, 0),
                 self.id(4, R - 2, 1),
                 self.id(4, R - 1, 1),
                 self.id(1, 0, 0),
                 self.id(1, 0, 1),
                 self.id(5, R - 1, 0),
                 self.id(5, R - 1, 1)),
                # x, y = R-1, R-1
                (self.id(4, R - 2, R - 1),
                 self.id(4, R - 2, R - 2),
                 self.id(4, R - 1, R - 2),
                 self.id(1, 0, R - 2),
                 self.id(1, 0, R - 1),
                 self.id(0, R - 1, R - 2),
                 self.id(0, R - 1, R - 1))
            ),
            # root pixel 5
            (
                # x, y = 0, 0
                (self.id(5, 1, 0),
                 self.id(5, 1, 1),
                 self.id(5, 0, 1),
                 self.id(1, R - 2, 0),
                 self.id(1, R - 1, 0),
                 self.id(2, 0, 0),
                 self.id(2, 1, 0)),
                # x, y = 0, R-1
                (self.id(5, 1, R - 1),
                 self.id(5, 1, R - 2),
                 self.id(5, 0, R - 2),
                 self.id(2, R - 2, 0),
                 self.id(2, R - 1, 0),
                 self.id(3, 0, 0),
                 self.id(3, 1, 0)),
                # x, y = R-1, 0
                (self.id(5, R - 2, 0),
                 self.id(5, R - 2, 1),
                 self.id(5, R - 1, 1),
                 self.id(1, 0, 0),
                 self.id(1, 1, 0),
                 self.id(4, R - 2, 0),
                 self.id(4, R - 1, 0)),
                # x, y = R-1, R-1
                (self.id(5, R - 2, R - 1),
                 self.id(5, R - 2, R - 2),
                 self.id(5, R - 1, R - 2),
                 self.id(3, R - 2, 0),
                 self.id(3, R - 1, 0),
                 self.id(4, 0, 0),
                 self.id(4, 1, 0))
            ),
        )

    def getGeometry(self, pixelId, fiducial=False):
        """Returns a spherical convex polygon corresponding to the fiducial
        or padded boundaries of the sky-pixel with the specified id.
        """
        root, ix, iy = self.coords(pixelId)
        xl = 2.0 * float(ix) / float(self.resolution) - 1.0
        xr = 2.0 * float(ix + 1) / float(self.resolution) - 1.0
        yb = 2.0 * float(iy) / float(self.resolution) - 1.0
        yt = 2.0 * float(iy + 1) / float(self.resolution) - 1.0
        c, x, y = self.center[root], self.x[root], self.y[root]
        v = map(g.normalize, [(c[0] + xl * x[0] + yb * y[0],
                               c[1] + xl * x[1] + yb * y[1],
                               c[2] + xl * x[2] + yb * y[2]),
                              (c[0] + xr * x[0] + yb * y[0],
                               c[1] + xr * x[1] + yb * y[1],
                               c[2] + xr * x[2] + yb * y[2]),
                              (c[0] + xr * x[0] + yt * y[0],
                               c[1] + xr * x[1] + yt * y[1],
                               c[2] + xr * x[2] + yt * y[2]),
                              (c[0] + xl * x[0] + yt * y[0],
                               c[1] + xl * x[1] + yt * y[1],
                               c[2] + xl * x[2] + yt * y[2])])
        if not fiducial and self.padding > 0.0:
            # Determine angles by which edge planes must be rotated outwards
            sp = math.sin(math.radians(self.padding))
            theta = map(lambda x: 0.5 * g.cartesianAngularSep(x[0], x[1]),
                        ((v[0],v[3]), (v[1],v[0]), (v[2],v[1]), (v[3],v[2])))
            sina = map(lambda x: sp / math.cos(math.radians(x)), theta)
            cosa = map(lambda x: math.sqrt(1.0 - x * x), sina)
            # find plane equations of fiducial pixel boundaries
            xlp = self.xplane[root][ix][0]
            ybp = self.yplane[root][iy][0]
            xrp = self.xplane[root][ix + 1][0]
            ytp = self.yplane[root][iy + 1][0]
            # rotate edge planes outwards
            xlp = self.xrot[root](xlp, -sina[0], cosa[0])
            ybp = self.yrot[root](ybp, -sina[1], cosa[1])
            xrp = self.xrot[root](xrp, sina[2], cosa[2])
            ytp = self.yrot[root](ytp, sina[3], cosa[3])
            # intersect rotated planes to find vertices of padded sky-pixel
            v = map(g.normalize, [g.cross(xlp, ybp),
                                  g.cross(xrp, ybp),
                                  g.cross(xrp, ytp),
                                  g.cross(xlp, ytp)])
        return g.SphericalConvexPolygon(v)

    def getCenter(self, pixelId):
        """Returns the center of a sky-pixel as a unit cartesian 3-vector.
        """
        if not isinstance(pixelId, (int, long)):
            raise TypeError(
                'Sky-pixel id must be an int or long')
        root, ix, iy = self.coords(pixelId)
        xc = 2.0 * (float(ix) + 0.5) / float(self.resolution) - 1.0
        yc = 2.0 * (float(iy) + 0.5) / float(self.resolution) - 1.0
        c, x, y = self.center[root], self.x[root], self.y[root]
        return g.normalize((c[0] + xc * x[0] + yc * y[0],
                            c[1] + xc * x[1] + yc * y[1],
                            c[2] + xc * x[2] + yc * y[2]))

    def getNeighbors(self, pixelId):
        """Returns a list of ids for sky-pixels adjacent to specified pixel.
        """
        root, ix, iy = self.coords(pixelId)
        R = self.resolution
        if ix == 0:
            # left edge
            if iy == 0:
                return self.cornerNeighbors[root][0]
            elif iy == R - 1:
                return self.cornerNeighbors[root][1]
            elif root == 0:
                neighbors = (self.id(2, R - iy - 2, R - 1),
                             self.id(2, R - iy - 1, R - 1),
                             self.id(2, R - iy,     R - 1))
            elif root == 5:
                neighbors = (self.id(2, iy - 1, 0),
                             self.id(2, iy,     0),
                             self.id(2, iy + 1, 0))
            else:
                r = root - 1
                if r == 0: r = 4
                neighbors = (self.id(r, R - 1, iy - 1),
                             self.id(r, R - 1, iy    ),
                             self.id(r, R - 1, iy + 1))
            return neighbors + (self.id(root, 1, iy - 1),
                                self.id(root, 1, iy    ),
                                self.id(root, 1, iy + 1),
                                self.id(root, 0, iy - 1),
                                self.id(root, 0, iy + 1))
        elif ix == R - 1:
            # right edge
            if iy == 0:
                return self.cornerNeighbors[root][2]
            elif iy == R - 1:
                return self.cornerNeighbors[root][3]
            elif root == 0:
                neighbors = (self.id(4, iy - 1, R - 1),
                             self.id(4, iy,     R - 1),
                             self.id(4, iy + 1, R - 1))
            elif root == 5:
                neighbors = (self.id(4, R - iy - 2, 0),
                             self.id(4, R - iy - 1, 0),
                             self.id(4, R - iy,     0))
            else:
                r = root + 1
                if r == 5: r = 1
                neighbors = (self.id(r, 0, iy - 1),
                             self.id(r, 0, iy    ),
                             self.id(r, 0, iy + 1))
            return neighbors + (self.id(root, R - 2, iy - 1),
                                self.id(root, R - 2, iy),
                                self.id(root, R - 2, iy + 1),
                                self.id(root, R - 1, iy - 1),
                                self.id(root, R - 1, iy + 1))
        elif iy == 0:
            # bottom edge
            if root == 0:
                neighbors = (self.id(3, ix - 1, R - 1),
                             self.id(3, ix    , R - 1),
                             self.id(3, ix + 1, R - 1))
            elif root == 1:
                neighbors = (self.id(5, R - ix - 2, 0),
                             self.id(5, R - ix - 1, 0),
                             self.id(5, R - ix,     0))
            elif root == 2:
                neighbors = (self.id(5, 0, ix - 1),
                             self.id(5, 0, ix    ),
                             self.id(5, 0, ix + 1))
            elif root == 3:
                neighbors = (self.id(5, ix - 1, R - 1),
                             self.id(5, ix    , R - 1),
                             self.id(5, ix + 1, R - 1))
            elif root == 4:
                neighbors = (self.id(5, R - 1, R - ix - 2),
                             self.id(5, R - 1, R - ix - 1),
                             self.id(5, R - 1, R - ix    ))
            else:
                neighbors = (self.id(1, R - ix - 2, 0),
                             self.id(1, R - ix - 1, 0),
                             self.id(1, R - ix,     0))
            return neighbors + (self.id(root, ix - 1, 1),
                                self.id(root, ix,     1),
                                self.id(root, ix + 1, 1),
                                self.id(root, ix - 1, 0),
                                self.id(root, ix + 1, 0))
        elif iy == R - 1:
            # top edge
            if root == 0:
                neighbors = (self.id(1, R - ix - 2, R - 1),
                             self.id(1, R - ix - 1, R - 1),
                             self.id(1, R - ix,     R - 1))
            elif root == 1:
                neighbors = (self.id(0, R - ix - 2, R - 1),
                             self.id(0, R - ix - 1, R - 1),
                             self.id(0, R - ix,     R - 1))
            elif root == 2:
                neighbors = (self.id(0, 0, R - ix - 2),
                             self.id(0, 0, R - ix - 1),
                             self.id(0, 0, R - ix    ))
            elif root == 3:
                neighbors = (self.id(0, ix - 1, 0),
                             self.id(0, ix,     0),
                             self.id(0, ix + 1, 0))
            elif root == 4:
                neighbors = (self.id(0, R - 1, ix - 1),
                             self.id(0, R - 1, ix    ),
                             self.id(0, R - 1, ix + 1))
            else:
                neighbors = (self.id(3, ix - 1, 0),
                             self.id(3, ix,     0),
                             self.id(3, ix + 1, 0))
            return neighbors + (self.id(root, ix - 1, R - 2),
                                self.id(root, ix,     R - 2),
                                self.id(root, ix + 1, R - 2),
                                self.id(root, ix - 1, R - 1),
                                self.id(root, ix + 1, R - 1))
        # interior pixel
        return (self.id(root, ix - 1, iy - 1),
                self.id(root, ix    , iy - 1),
                self.id(root, ix + 1, iy - 1),
                self.id(root, ix - 1, iy),
                self.id(root, ix + 1, iy),
                self.id(root, ix - 1, iy + 1),
                self.id(root, ix    , iy + 1),
                self.id(root, ix + 1, iy + 1))

    def intersect(self, polygon):
        """Returns a list of ids for sky-pixels intersecting the given
        polygon. Intersection is computed against padded sky-pixels.
        """
        if not isinstance(polygon, g.SphericalConvexPolygon):
            raise TypeError('Expecting a SphericalConvexPolygon as input ' +
                            'to the sky-pixel intersection computation')
        pixels = []
        R = self.resolution
        for root in xrange(6):
            # clip against padded root pixel edge planes
            poly = polygon
            for p in (self.xplane[root][0][2], self.xplane[root][R][1],
                      self.yplane[root][0][2], self.yplane[root][R][1]):
                poly = poly.clip(p)
                if poly == None:
                    break
            if poly == None:
                # polygon doesn't intersect this root triangle
                continue
            self._intersect(pixels, poly, root, (0, R, 0, R))
        return pixels

    def __len__(self):
        """Returns the total number of sky-pixels in this pixelization.
        """
        return 6 * self.resolution * self.resolution

    def __iter__(self):
        """Returns an iterator over the sky-pixel ids of all the pixels
        in this pixelization.
        """
        return xrange(len(self))

    def __repr__(self):
        return ''.join([self.__class__.__name__, '(',
                        repr(self.resolution), ', ', repr(self.padding), ')'])

    def id(self, root, x, y):
        """Maps from a root pixel number and x, y pixel coordinates
        to a sky-pixel id.
        """
        if not all(isinstance(p, (int, long)) for p in (root, x, y)):
            raise TypeError('Pixel coordinates must all be of type int or long')
        if root < 0 or root > 5:
            raise RuntimeError('root pixel number must be in range [0,6)')
        R = self.resolution
        if x < 0 or x >= R or y < 0 or y >= R:
            raise RuntimeError('x and y pixel coordinates must ' +
                               'be in range [0,%d)' % R)
        if root > 0 and root < 5:
            # easy equatorial pixel case
            return R ** 2 + 4 * R * y + (root - 1) * R + x
        # assign ids to pixels in expanding concentric squares
        dx = x - 0.5 * (R - 1)
        dy = y - 0.5 * (R - 1)
        ring = max(abs(dx), abs(dy))
        r = 2.0 * ring - 1.0
        if r < 0:
            if root == 5:
                return 0
            else:
                return 6 * R ** 2 - 1
        if dy == -ring:
            i = int(r * r + ring - dx)
        elif dx == -ring:
            i = int(r * r + 3.0 * ring + dy)
        elif dy == ring:
            i = int(r * r + 5.0 * ring + dx)
        else: # dx == ring
            i = int(r * r + 7.0 * ring - dy)
        if root == 5:
            # south pole - done
            return i
        # north pole: use shrinking concentric squares and wind
        # in the opposite direction
        return 6 * R ** 2 - i - 1

    def coords(self, pixelId):
        """Maps from a sky-pixel id to a root pixel number and x, y 
        pixel coordinates.
        """
        if not isinstance(pixelId, (int, long)):
            raise TypeError(
                'Sky-pixel id must be an int or long')        
        R = self.resolution
        R2 = R ** 2
        if pixelId < 0 or pixelId >= 6 * R2:
            raise RuntimeError(
                'Invalid sky-pixel id %d - value must be in range [0, %d)' %
                (pixelId, 6 * R2))
        r = pixelId / R2
        if r > 0 and r < 5:
            # equatorial pixel case
            pixelId -= R2
            y = pixelId / (4 * R)
            pixelId -= 4 * R * y
            root = 1 + pixelId / R
            x = pixelId - (root - 1) * R
            return (root, x, y)
        # polar root pixels
        if r == 5:
            root = 0
            i = 6 * R2 - pixelId - 1
        else:
            root = 5
            i = pixelId
        s = int(math.floor(math.sqrt(i)))
        if R & 1 == 0:
            ring = 0.5 * (s + ((s & 1) ^ 1))
        else:
            ring = 0.5 * (s + (s & 1))
        r = 2.0 * ring - 1.0
        if r < 0.0:
            return (root, R / 2, R / 2)
        i -= r * r
        if i <= 2.0 * ring:
            dx = ring - i
            dy = -ring
        elif i <= 4.0 * ring:
            dx = -ring
            dy = i - 3.0 * ring
        elif i <= 6.0 * ring:
            dx = i - 5.0 * ring
            dy = ring
        else:
            dx = ring
            dy = 7.0 * ring - i
        return (root, int(dx + 0.5 * (R - 1)), int(dy + 0.5 * (R - 1)))

    def _fiducialXPlane(self, root, ix):
        assert isinstance(ix, (int,long)) and ix >= 0 and ix <= self.resolution
        x = 2.0 * float(ix) / float(self.resolution) - 1.0
        c, b = self.center[root], self.x[root]
        v = (c[0] + x * b[0], c[1] + x * b[1], c[2] + x * b[2])
        return g.normalize(self.xrot[root](v, 1.0, 0.0))

    def _fiducialYPlane(self, root, iy):
        assert isinstance(iy, (int,long)) and iy >= 0 and iy <= self.resolution
        y = 2.0 * float(iy) / float(self.resolution) - 1.0
        c, b = self.center[root], self.y[root]
        v = (c[0] + y * b[0], c[1] + y * b[1], c[2] + y * b[2])
        return g.normalize(self.yrot[root](v, 1.0, 0.0))

    def _intersect(self, pixels, poly, root, box):
        dx = box[1] - box[0]
        dy = box[3] - box[2]
        if dx == 1 and dy == 1:
            i = self.id(root, box[0], box[2])
            if poly.intersects(self.getGeometry(i)):
                pixels.append(i)
            return
        if dx >= dy:
            # Split xrange and recurse
            xsplit = box[0] + dx / 2
            p = poly.clip(self.xplane[root][xsplit][1])
            if p != None:
                self._overlap(pixels, p, root, (box[0], xsplit, box[2], box[3]))
            p = poly.clip(self.xplane[root][xsplit][2])
            if p != None:
                self._overlap(pixels, p, root, (xsplit, box[1], box[2], box[3]))
        else:
            # Split yrange and recurse
            ysplit = box[2] + dy / 2
            p = poly.clip(self.yplane[root][ysplit][1])
            if p != None:
                self._overlap(pixels, p, root, (box[0], box[1], box[2], ysplit))
            p = poly.clip(self.yplane[root][ysplit][2])
            if p != None:
                self._overlap(pixels, p, root, (box[0], box[1], ysplit, box[3]))


def imageToSkyPixels(pixelization, wcs, naxis1, naxis2, pad=0.0):
    """Given a sky pixelization, a WCS, a pair of dimensions, and an
    optional padding amount (in units of pixels), compute the list of ids of
    sky-pixel intersecting the image.
    """
    # Compute image corners
    xmin = -0.5 - pad
    ymin = xmin
    xmax = naxis1 - 0.5 + pad
    ymax = naxis2 - 0.5 + pad
    # Produce a lsst.afw.coord.coordLib.Coord object for each vertex
    verts = [wcs.pixelToSky(xmin, ymin), wcs.pixelToSky(xmax, ymin),
             wcs.pixelToSky(xmin, ymax), wcs.pixelToSky(xmax, ymax)]
    # Map these to cartesian unit vectors
    verts = map(g.cartesianUnitVector,
                ((v.getLongitude(coord.DEGREES), v.getLatitude(coord.DEGREES))
                 for v in verts))
    # then turn them into a spherical convex polygon
    convex, cc = g.convex(verts)
    if not convex:
        raise RuntimeError('Image corners do not form a convex polygon: ' + cc)
    elif not cc:
        verts.reverse()
    return pixelization.intersect(g.SphericalConvexPolygon(verts))

