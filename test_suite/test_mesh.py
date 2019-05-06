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
    Afunc = lambda R,Z: R*Z
    A_xpoint = 0.
    c = MeshContour([Point2D(R,Z) for R,Z in [Point2D(0.,0.), Point2D(1.,0.),
        Point2D(1.,1.), Point2D(0.,1.)]], Afunc, A_xpoint)

    def test_append(self):
        copy = deepcopy(self.c)
        copy.append(Point2D(0.,0.))
        assert copy.distance[-1] == tight_approx(4.)

    def test_distance(self):
        assert self.c.distance == tight_approx(numpy.arange(4))
