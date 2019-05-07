import numpy
import pytest
from copy import deepcopy
from mesh import *

def tight_approx(b):
    return pytest.approx(b, rel=1.e-12, abs=1.e-13)

class TestPoints:
    p0 = Point2D(1., 2.)
    p1 = Point2D(3., 4.)

    def test_add(self):
        p = self.p0 + self.p1
        assert p.R == tight_approx(4.)
        assert p.Z == tight_approx(6.)

    def test_sub(self):
        p = self.p0 - self.p1
        assert p.R == tight_approx(-2.)
        assert p.Z == tight_approx(-2.)

    def test_mul(self):
        p = self.p0 * 1.5
        assert p.R == tight_approx(1.5)
        assert p.Z == tight_approx(3.)

    def test_rmul(self):
        p = 1.5 * self.p0
        assert p.R == tight_approx(1.5)
        assert p.Z == tight_approx(3.)

    def test_div(self):
        p = self.p0 / 3.
        assert p.R == tight_approx(1./3.)
        assert p.Z == tight_approx(2./3.)

    def test_iter(self):
        assert [x for x in self.p0] == tight_approx([1.,2.])

    def test_repr(self):
        assert str(self.p0) == "Point2D(1.0,2.0)"

    def test_distance(self):
        assert calc_distance(self.p0, self.p1) == tight_approx(2.*numpy.sqrt(2.))

class TestContour:
    @pytest.fixture
    def testcontour(self):
        psifunc = lambda R,Z: R**2 + Z**4
        psi_xpoint = 0.

        # make a circle, not centred on origin
        class returnObject:
            npoints = 23
            r = 1.
            R0 = .2
            Z0 = .3
            theta = numpy.linspace(0.,2.*numpy.pi,npoints)
            R = R0 + r*numpy.cos(theta)
            Z = Z0 + r*numpy.sin(theta)
            c = MeshContour([Point2D(R,Z) for R,Z in zip(R,Z)], psifunc, psi_xpoint)

        return returnObject()

    def test_distance(self, testcontour):
        segment_length = 2.*testcontour.r*numpy.sin(2.*numpy.pi/(testcontour.npoints-1)/2.)
        assert testcontour.c.distance == tight_approx(segment_length*numpy.arange(23))


    def test_iter(self, testcontour):
        clist = list(testcontour.c)

        for i,item in enumerate(clist):
            assert item.R == tight_approx(testcontour.R[i])
            assert item.Z == tight_approx(testcontour.Z[i])

    def test_getitem(self, testcontour):
        p = testcontour.c[5]
        assert p.R == tight_approx(testcontour.R[5])
        assert p.Z == tight_approx(testcontour.Z[5])

    def test_append(self, testcontour):
        c = testcontour.c
        expected_distance = c.distance[-1] + numpy.sqrt((1.-c[-1].R)**2 + (1.-c[-1].Z)**2)
        c.append(Point2D(1.,1.))
        assert c.distance[-1] == tight_approx(expected_distance)

    def test_reverse(self, testcontour):
        c = testcontour.c
        orig = deepcopy(c)

        c.reverse()

        n = len(orig)
        total_d = orig.distance[-1]
        for i in range(n):
            assert orig[n-1-i].R == tight_approx(c[i].R)
            assert orig[n-1-i].Z == tight_approx(c[i].Z)
            assert total_d - orig.distance[n-1-i] == tight_approx(c.distance[i])

    def test_refine(self, testcontour):
        # MeshContour.refine just calls MeshContour.getRefined, so this tests both

        c = testcontour.c
        c.psival = .7
        # c does not start close to a contour of psi_func, so need to use a large width
        c.refine(width=2., atol=1.e-13)
        for p in c:
            assert c.psi(p.R, p.Z) == tight_approx(.7)
