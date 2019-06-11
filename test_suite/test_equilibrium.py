import numpy
import pytest
from copy import deepcopy
from ..equilibrium import *
from .utils_for_tests import *

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

    def test_as_ndarray(self):
        p = self.p0.as_ndarray()
        assert issubclass(type(p), numpy.ndarray)
        assert p == tight_approx(numpy.array([1., 2.]))

def test_find_intersectionRR1():
    l1start = Point2D(-1., -.1)
    l1end = Point2D(1., .1)
    l2start = Point2D(-1., .1)
    l2end = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionRR2():
    l1end = Point2D(-1., -.1)
    l1start = Point2D(1., .1)
    l2start = Point2D(-1., .1)
    l2end = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionRR3():
    l1start = Point2D(-1., -.1)
    l1end = Point2D(1., .1)
    l2end = Point2D(-1., .1)
    l2start = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionRR4():
    l1end = Point2D(-1., -.1)
    l1start = Point2D(1., .1)
    l2end = Point2D(-1., .1)
    l2start = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionRRNone1():
    l1start = Point2D(2., -.1)
    l1end = Point2D(4., .1)
    l2start = Point2D(-1., .1)
    l2end = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionRRNone2():
    l1end = Point2D(2., -.1)
    l1start = Point2D(4., .1)
    l2start = Point2D(-1., .1)
    l2end = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionRRNone3():
    l1start = Point2D(2., -.1)
    l1end = Point2D(4., .1)
    l2end = Point2D(-1., .1)
    l2start = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionRRNone4():
    l1end = Point2D(2., -.1)
    l1start = Point2D(4., .1)
    l2end = Point2D(-1., .1)
    l2start = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionRZ1():
    l1start = Point2D(-1., -.1)
    l1end = Point2D(1., .1)
    l2start = Point2D(-.1, 1.)
    l2end = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionRZ2():
    l1end = Point2D(-1., -.1)
    l1start = Point2D(1., .1)
    l2start = Point2D(-.1, 1.)
    l2end = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionRZ3():
    l1start = Point2D(-1., -.1)
    l1end = Point2D(1., .1)
    l2end = Point2D(-.1, 1.)
    l2start = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionRZ4():
    l1end = Point2D(-1., -.1)
    l1start = Point2D(1., .1)
    l2end = Point2D(-.1, 1.)
    l2start = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionRZNone1():
    l1start = Point2D(2., -.1)
    l1end = Point2D(4., .1)
    l2start = Point2D(-.1, 1.)
    l2end = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionRZNone2():
    l1end = Point2D(2., -.1)
    l1start = Point2D(4., .1)
    l2start = Point2D(-.1, 1.)
    l2end = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionRZNone3():
    l1start = Point2D(2., -.1)
    l1end = Point2D(4., .1)
    l2end = Point2D(-.1, 1.)
    l2start = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionRZNone4():
    l1end = Point2D(2., -.1)
    l1start = Point2D(4., .1)
    l2end = Point2D(-.1, 1.)
    l2start = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionZR1():
    l1start = Point2D(-.1, -1.)
    l1end = Point2D(.1, 1.)
    l2start = Point2D(-1., .1)
    l2end = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionZR2():
    l1end = Point2D(-.1, -1.)
    l1start = Point2D(.1, 1.)
    l2start = Point2D(-1., .1)
    l2end = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionZR3():
    l1start = Point2D(-.1, -1.)
    l1end = Point2D(.1, 1.)
    l2end = Point2D(-1., .1)
    l2start = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionZR4():
    l1end = Point2D(-.1, -1.)
    l1start = Point2D(.1, 1.)
    l2end = Point2D(-1., .1)
    l2start = Point2D(1., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionZRNone1():
    l1start = Point2D(-.1, -1.)
    l1end = Point2D(.1, 1.)
    l2start = Point2D(2., .1)
    l2end = Point2D(4., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionZRNone2():
    l1end = Point2D(-.1, -1.)
    l1start = Point2D(.1, 1.)
    l2start = Point2D(2., .1)
    l2end = Point2D(4., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionZRNone3():
    l1start = Point2D(-.1, -1.)
    l1end = Point2D(.1, 1.)
    l2end = Point2D(2., .1)
    l2start = Point2D(4., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionZRNone4():
    l1end = Point2D(-.1, -1.)
    l1start = Point2D(.1, 1.)
    l2end = Point2D(2., .1)
    l2start = Point2D(4., -.1)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionZZ1():
    l1start = Point2D(-.1, -1.)
    l1end = Point2D(.1, 1.)
    l2start = Point2D(-.1, 1.)
    l2end = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionZZ2():
    l1end = Point2D(-.1, -1.)
    l1start = Point2D(.1, 1.)
    l2start = Point2D(-.1, 1.)
    l2end = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionZZ3():
    l1start = Point2D(-.1, -1.)
    l1end = Point2D(.1, 1.)
    l2end = Point2D(-.1, 1.)
    l2start = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionZZ4():
    l1end = Point2D(-.1, -1.)
    l1start = Point2D(.1, 1.)
    l2end = Point2D(-.1, 1.)
    l2start = Point2D(.1, -1.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect.R == tight_approx(0.)
    assert intersect.Z == tight_approx(0.)

def test_find_intersectionZZNone1():
    l1start = Point2D(-.1, -1.)
    l1end = Point2D(.1, 1.)
    l2start = Point2D(-.1, 4.)
    l2end = Point2D(.1, 2.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionZZNone2():
    l1end = Point2D(-.1, -1.)
    l1start = Point2D(.1, 1.)
    l2start = Point2D(-.1, 4.)
    l2end = Point2D(.1, 2.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionZZNone3():
    l1start = Point2D(-.1, -1.)
    l1end = Point2D(.1, 1.)
    l2end = Point2D(-.1, 4.)
    l2start = Point2D(.1, 2.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

def test_find_intersectionZZNone4():
    l1end = Point2D(-.1, -1.)
    l1start = Point2D(.1, 1.)
    l2end = Point2D(-.1, 4.)
    l2start = Point2D(.1, 2.)
    intersect = find_intersection(l1start, l1end, l2start, l2end)
    assert intersect == None

class TestContour:

    @pytest.fixture
    def testcontour(self):
        # make a circle, not centred on origin
        class returnObject:
            npoints = 23
            r = 1.
            R0 = .2
            Z0 = .3
            theta = numpy.linspace(0.,numpy.pi,npoints)

            def __init__(self):
                psifunc = lambda R,Z: (R - self.R0)**2 + (Z - self.Z0)**2

                self.R = self.R0 + self.r*numpy.cos(self.theta)
                self.Z = self.Z0 + self.r*numpy.sin(self.theta)

                psi_xpoint = psifunc(self.R[0], self.Z[0])

                self.c = PsiContour([Point2D(R,Z) for R,Z in zip(self.R,self.Z)], psifunc,
                        psi_xpoint)
                self.c.refine_width = 1.e-3

        return returnObject()

    def test_distance(self, testcontour):
        segment_length = testcontour.r*numpy.pi/(testcontour.npoints-1)
        c = testcontour.c
        assert testcontour.c.distance == pytest.approx(
                segment_length*numpy.arange(testcontour.npoints), abs=1.e-4)

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
        point_to_add = c[-2]
        expected_distance = c.distance[-2]
        del c.points[-1]
        del c.points[-1]
        c.endInd = -1
        c.append(point_to_add)
        assert c.distance[-1] == pytest.approx(expected_distance, abs=4.e-5)

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
        # PsiContour.refine just calls PsiContour.getRefined, so this tests both

        c = testcontour.c
        c.psival = .7
        # c does not start close to a contour of psi_func, so need to use a large width
        c.refine(width=2., atol=1.e-13)
        for p in c:
            assert c.psi(p.R, p.Z) == tight_approx(.7)

    def test_interpFunction(self, testcontour):
        f = testcontour.c.interpFunction()
        p = f(numpy.pi*testcontour.r)
        assert p.R == pytest.approx(testcontour.R0 - testcontour.r, abs=1.e-9)
        assert p.Z == pytest.approx(testcontour.Z0, abs=1.e-5)

    def test_getRegridded(self, testcontour):
        orig = testcontour.c
        r = testcontour.r

        newNpoints = 97
        sfunc_true = lambda i: numpy.sqrt(i / (newNpoints - 1)) * numpy.pi * r
        sfunc = lambda i: numpy.sqrt(i / (newNpoints - 1)) * orig.distance[-1]
        newTheta = sfunc_true(numpy.arange(newNpoints)) / r
        newR = testcontour.R0 + r*numpy.cos(newTheta)
        newZ = testcontour.Z0 + r*numpy.sin(newTheta)

        new = orig.getRegridded(newNpoints, sfunc=sfunc, width=1.e-3)

        assert [p.R for p in new] == pytest.approx(newR, abs=4.e-4)
        assert [p.Z for p in new] == pytest.approx(newZ, abs=4.e-4)

    def test_getRegridded_extend(self, testcontour):
        orig = testcontour.c

        new = orig.getRegridded(testcontour.npoints, width=.1, extend_lower=1, extend_upper=2)

        assert numpy.array([[*p] for p in new[1:-2]]) == pytest.approx(numpy.array([[*p] for p in orig]), abs=1.e-8)

        # test the extend_lower
        assert [*new[0]] == pytest.approx(
                [orig[1].R, 2.*testcontour.Z0-orig[1].Z], abs=1.e-10)

        # test the extend_upper
        assert [*new[-2]] == pytest.approx(
                [orig[-2].R, 2.*testcontour.Z0-orig[-2].Z], abs=1.e-10)
        assert [*new[-1]] == pytest.approx(
                [orig[-3].R, 2.*testcontour.Z0-orig[-3].Z], abs=1.e-10)

    def test_contourSfunc(self, testcontour):
        c = testcontour.c
        c.startInd = 2
        c.endInd = len(c) - 2
        n = c.endInd - c.startInd + 1

        f = c.contourSfunc()

        indices = numpy.arange(n, dtype=float)
        assert f(indices) == tight_approx([d - c.distance[c.startInd] for d in
                                           c.distance[c.startInd:c.endInd+1]])
        assert f(-1.) == tight_approx(0.)
        assert f(n + 1.) == tight_approx(c.distance[c.endInd] - c.distance[c.startInd])

    def test_contourSfunc_list(self, testcontour):
        c1 = testcontour.c
        c1.startInd = 2
        c1.endInd = len(c1) - 2
        n = c1.endInd - c1.startInd

        npoints = len(c1)
        r = testcontour.r
        R0 = testcontour.R0
        Z0 = testcontour.Z0
        theta = numpy.linspace(0.,1.5*numpy.pi,npoints)

        R = R0 + r*numpy.cos(theta)
        Z = Z0 + r*numpy.sin(theta)

        c2 = PsiContour([Point2D(R,Z) for R,Z in zip(R,Z)], c1.psi,
                c1.psival)
        c2.refine_width = c1.refine_width
        c2.startInd = 2
        c2.endInd = len(c2) - 1

        c_list = [c1, c2]
        sfunc_list = []

        # This does NOT work as desired - when the function added to 'sfunc_list' is
        # evaluated, it takes 'sfunc_orig' from the enclosing scope, which means it uses
        # the last value from the loop instead of the value when 'sfunc' was added to
        # 'sfunc_list'
        for c in c_list:
            sfunc_orig = c.contourSfunc()
            sfunc = lambda i: sfunc_orig(i) - 3.
            sfunc_list.append(sfunc)

        # notice we check that the first test *fails*
        assert not sfunc_list[0](float(n)) == pytest.approx(n/(npoints - 1.) * numpy.pi*r - 3., abs=1.e-6)
        assert sfunc_list[1](float(n)) == pytest.approx(n/(npoints - 1.) * 1.5*numpy.pi*r - 3., abs=4.e-6)

        sfunc_list2 = []
        # This version does work, because when the lambda is evaluated it uses
        # 'sfunc_orig' from the scope of 'shift_sfunc' in which it was created.
        def shift_sfunc(c):
            sfunc_orig = c.contourSfunc()
            return lambda i: sfunc_orig(i) - 3.
        for c in c_list:
            sfunc_list2.append(shift_sfunc(c))

        assert sfunc_list2[0](float(n)) == pytest.approx(n/(npoints - 1.) * numpy.pi*r - 3., abs=1.e-6)
        assert sfunc_list2[1](float(n)) == pytest.approx(n/(npoints - 1.) * 1.5*numpy.pi*r - 3., abs=4.e-6)

    def test_interpSSperp(self, testcontour):
        c = testcontour.c

        # Make c.startInd > 0
        c.insert(0, Point2D(c[1].R, 2.*c[0].Z - c[1].Z))

        # 'vec' argument is in Z-direction, so 's_perp' is displacement in R-direction
        sfunc, s_perp_total = c.interpSSperp([0., 1.])
        assert sfunc(0.) == tight_approx(0.)
        assert sfunc(1.) == pytest.approx(numpy.pi/2., abs=1.e-6)
        assert sfunc(2.) == pytest.approx(numpy.pi, abs=2.e-6)
        assert s_perp_total == tight_approx(2.)

    def test_FineContour(self, testcontour):
        testcontour.c.refine_width = 1.e-2

        fc = FineContour(testcontour.c)

        assert fc.totalDistance() == pytest.approx(numpy.pi, abs=1.e-5)

        interpFunc = fc.interpFunction()
        r = testcontour.r
        for theta in numpy.linspace(0., numpy.pi, 17):
            p = interpFunc(r*theta)
            assert p.R == pytest.approx(testcontour.R0 + r*numpy.cos(theta), abs=1.e-4)
            assert p.Z == pytest.approx(testcontour.Z0 + r*numpy.sin(theta), abs=1.e-4)

class TestEquilibrium:

    @pytest.fixture
    def eq(self):
        eq = Equilibrium()
        eq.wall = [Point2D(-1., -1.), Point2D(1., -1.), Point2D(1., 1.), Point2D(-1., 1.)]
        return eq

    def test_make1dGrid(self, eq):
        n = 4
        f = lambda i: i**2
        r = eq.make1dGrid(n, f)
        assert r == tight_approx([0., 0.5, 1., 2.5, 4., 6.5, 9., 12.5, 16.])

    def test_wallIntersection(self, eq):
        intersect = eq.wallIntersection(Point2D(0., 0.), Point2D(2., 0.))
        assert intersect.R == tight_approx(1.)
        assert intersect.Z == tight_approx(0.)

        # Line goes exactly through a corner of the wall, so intersects with two wall
        # segments
        intersect = eq.wallIntersection(Point2D(0., 0.), Point2D(2., 2.))
        assert intersect.R == tight_approx(1.)
        assert intersect.Z == tight_approx(1.)

class TestEquilibriumRegion:

    @pytest.fixture
    def eqReg(self):
        options = {'nx':1, 'ny':5, 'y_boundary_guards':1}
        equilib = Equilibrium()
        equilib.psi = lambda R,Z: R - Z
        n = 11.
        points = [Point2D(i * 3. / (n - 1.), i * 3. / (n - 1.)) for i in numpy.arange(n)]
        eqReg = EquilibriumRegion(equilib, '', 1, options, points, equilib.psi, 0.)
        return eqReg

    def test_getPolynomialPoloidalDistanceFuncLinear(self, eqReg):
        L = 2.
        N = 10.
        N_norm = 1
        f = eqReg.getPolynomialPoloidalDistanceFunc(L, N, N_norm)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(10.) == tight_approx(2.)
        # f(i) = i/N*L
        assert f(3.) == tight_approx(0.6)

    def test_getPolynomialPoloidalDistanceFuncDLower(self, eqReg):
        d_lower = 0.01
        L = 2.
        N = 10.
        N_norm = 40.
        f = eqReg.getPolynomialPoloidalDistanceFunc(L, N, N_norm, d_lower = d_lower)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f = d_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(d_lower*itest/N_norm, abs=1.e-5)
        # for i<<1, d2f/di2 = 0
        assert f(itest) - f(0) == pytest.approx(f(2.*itest) - f(itest), abs=1.e-5)

    def test_getPolynomialPoloidalDistanceFuncDUpper(self, eqReg):
        d_upper = 0.01
        L = 2.
        N = 10.
        N_norm = 40.
        f = eqReg.getPolynomialPoloidalDistanceFunc(L, N, N_norm, d_upper = d_upper)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for N-i<<1, f = L - d_upper*i/N_norm
        itest = 0.01
        assert f(N - itest) == pytest.approx(L - d_upper*itest/N_norm, abs=1.e-5)
        # for N-i<<1, d2f/di2 = 0
        assert f(N - itest) - f(N) == pytest.approx(f(N - 2.*itest) - f(N - itest), abs=1.e-5)

    def test_getPolynomialPoloidalDistanceFuncDBoth(self, eqReg):
        d_lower = 0.02
        d_upper = 0.01
        L = 2.
        N = 10.
        N_norm = 40.
        f = eqReg.getPolynomialPoloidalDistanceFunc(L, N, N_norm, d_lower = d_lower, d_upper = d_upper)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f = d_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(d_lower*itest/N_norm, abs=1.e-5)
        # for i<<1, d2f/di2 = 0
        assert f(itest) - f(0) == pytest.approx(f(2.*itest) - f(itest), abs=1.e-5)
        # for N-i<<1, f = L - d_upper*i/N_norm
        assert f(N - itest) == pytest.approx(L - d_upper*itest/N_norm, abs=1.e-5)
        # for N-i<<1, d2f/di2 = 0
        assert f(N - itest) - f(N) == pytest.approx(f(N - 2.*itest) - f(N - itest), abs=1.e-5)

    def test_getSqrtPoloidalDistanceFuncLinear(self, eqReg):
        L = 2.
        N = 10.
        N_norm = 1
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(10.) == tight_approx(2.)
        # f(i) = i/N*L
        assert f(3.) == tight_approx(0.6)

    def test_getSqrtPoloidalDistanceFuncBLower(self, eqReg):
        b_lower = 0.01
        L = 2.
        N = 10.
        N_norm = 40.
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm, b_lower = b_lower)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ 2*b_lower*sqrt(i/N_norm) + b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(2.*b_lower*numpy.sqrt(itest/N_norm)
                                         + b_lower*itest/N_norm, abs=1.e-5)

    def test_getSqrtPoloidalDistanceFuncBUpper(self, eqReg):
        b_upper = 0.01
        L = 2.
        N = 10.
        N_norm = 40.
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm, b_upper = b_upper)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for (N-i)<<1, f ~ L - 2*b_upper*sqrt((N-i)/N_norm) - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(L - 2.*b_upper*numpy.sqrt((N - itest)/N_norm)
                                         - b_upper*(N - itest)/N_norm, abs=1.e-5)

    def test_getSqrtPoloidalDistanceFuncBBoth(self, eqReg):
        b_lower = 0.1
        b_upper = 0.2
        L = 2.
        N = 10.
        N_norm = 40.
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm, b_lower = b_lower, b_upper = b_upper)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ 2*b_lower*sqrt(i/N_norm) + b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(2.*b_lower*numpy.sqrt(itest/N_norm)
                                         + b_lower*itest/N_norm, abs=1.e-5)
        # for (N-i)<<1, f ~ L - 2*b_upper*sqrt((N-i)/N_norm) - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(L - 2.*b_upper*numpy.sqrt((N - itest)/N_norm)
                                         - b_upper*(N - itest)/N_norm, abs=1.e-5)

    def test_getSqrtPoloidalDistanceFuncBothLower(self, eqReg):
        b_lower = 0.01
        a_lower = 0.05
        L = 2.
        N = 10.
        N_norm = 2.
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm, b_lower = b_lower, a_lower =
                a_lower)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ 2*a_lower*sqrt(i/N_norm) + b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(2.*a_lower*numpy.sqrt(itest/N_norm)
                                         + b_lower*itest/N_norm, abs=1.e-5)

    def test_getSqrtPoloidalDistanceFuncBothUpper(self, eqReg):
        b_upper = 0.01
        a_upper = 0.05
        L = 2.
        N = 10.
        N_norm = 1
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm, b_upper = b_upper, a_upper = a_upper)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for (N-i)<<1, f ~ L - 2*a_upper*sqrt((N-i)/N_norm) - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(L - 2.*a_upper*numpy.sqrt((N - itest)/N_norm)
                                         - b_upper*(N - itest)/N_norm, abs=1.e-5)

    def test_getSqrtPoloidalDistanceFuncBothBoth(self, eqReg):
        b_lower = 0.01
        a_lower = 0.05
        b_upper = 0.2
        a_upper = 0.07
        L = 2.
        L = 2.
        N = 10.
        N_norm = 1
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm, b_lower = b_lower,
                a_lower = a_lower, b_upper = b_upper,
                a_upper = a_upper)
        # f(0) = 0
        assert f(0.) == tight_approx(0.)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ 2*a_lower*sqrt(i/N_norm) + b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(2.*a_lower*numpy.sqrt(itest/N_norm)
                                         + b_lower*itest/N_norm, abs=1.e-5)
        # for (N-i)<<1, f ~ L - 2*a_upper*sqrt((N-i)/N_norm) - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(L - 2.*a_upper*numpy.sqrt((N - itest)/N_norm)
                                         - b_upper*(N - itest)/N_norm, abs=1.e-5)

    def test_combineSfuncsPoloidalSpacingNoD(self, eqReg):
        n = len(eqReg)
        L = eqReg.totalDistance()

        sfunc_orthogonal = lambda i: i/(n - 1.) * L

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(0.) == tight_approx(0.)
        assert sfunc((n - 1.)/3.) == tight_approx(L/3.)
        assert sfunc(n - 1.) == tight_approx(L)

    def test_combineSfuncsPoloidalSpacingRangeLower(self, eqReg):
        n = len(eqReg)
        L = eqReg.totalDistance()

        eqReg.poloidalSpacingParameters.nonorthogonal_range_lower = .1
        eqReg.poloidalSpacingParameters.N_norm = 40.

        sfunc_orthogonal = lambda i: numpy.piecewise(i, [i < 0.],
                [100., lambda i: numpy.sqrt(i/(n - 1.)) * L])

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(-1.) == tight_approx(-1./(n - 1.) * L)
        assert sfunc(0.) == tight_approx(0.)
        itest = 0.0001
        assert sfunc(itest) == pytest.approx(0.00001*3.*numpy.sqrt(2.), abs=1.e-11)
        assert sfunc(n-1.) == tight_approx(L)

    def test_combineSfuncsPoloidalSpacingRangeUpper(self, eqReg):
        n = len(eqReg)
        L = eqReg.totalDistance()

        eqReg.poloidalSpacingParameters.nonorthogonal_range_upper = .1
        eqReg.poloidalSpacingParameters.N_norm = 40.

        sfunc_orthogonal = lambda i: numpy.piecewise(i, [i > n - 1.],
                [100., lambda i: numpy.sqrt(i/(n - 1.)) * L])

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(0.) == tight_approx(0.)
        itest = n - 1. - 0.0001
        assert L - sfunc(itest) == pytest.approx(0.00001*3.*numpy.sqrt(2.), abs=1.e-11)
        assert sfunc(n - 1.) == tight_approx(L)
        assert sfunc(float(n)) == tight_approx(n/(n - 1.) * L)

    def test_combineSfuncsPoloidalSpacingRangeBoth(self, eqReg):
        n = len(eqReg)
        L = eqReg.totalDistance()

        eqReg.poloidalSpacingParameters.nonorthogonal_range_lower = .1
        eqReg.poloidalSpacingParameters.nonorthogonal_range_upper = .1
        eqReg.poloidalSpacingParameters.N_norm = 40.

        sfunc_orthogonal = lambda i: numpy.piecewise(i, [i < 0., i > n - 1.],
                [-100., 100., lambda i: numpy.sqrt(i/(n - 1.)) * L])

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(-1.) == tight_approx(-1./(n - 1.) * L)
        assert sfunc(0.) == tight_approx(0.)
        itest = 1.e-6
        assert sfunc(itest) == pytest.approx(3.e-7*numpy.sqrt(2.), rel=1.e-7)
        itest = n - 1. - 1.e-6
        assert L - sfunc(itest) == pytest.approx(3.e-7*numpy.sqrt(2.), rel=1.e-7)
        assert sfunc(n - 1.) == tight_approx(L)
        assert sfunc(float(n)) == tight_approx(n/(n - 1.) * L)

    def test_combineSfuncsPerpSpacingIntegrated(self, eqReg):
        # This test follows roughly the operations in
        # MeshRegion.distributePointsNonorthogonal for the default 'combined' option.
        n = len(eqReg)
        L = eqReg.totalDistance()

        eqReg.poloidalSpacingParameters.polynomial_d_lower = .1
        eqReg.poloidalSpacingParameters.polynomial_d_upper = .1
        eqReg.poloidalSpacingParameters.nonorthogonal_range_lower = .3
        eqReg.poloidalSpacingParameters.nonorthogonal_range_upper = .3
        eqReg.poloidalSpacingParameters.N_norm = 40

        sfunc_orthogonal_original = eqReg.contourSfunc()

        # as if lower_wall
        intersect_index = 2
        original_start = eqReg.startInd
        distance_at_original_start = eqReg.distance[original_start]

        d = 1.5*3. / (n - 1.)
        eqReg.points.insert(intersect_index, Point2D(d, d))

        if original_start >= intersect_index:
            original_start += 1

        distance_at_wall = eqReg.distance[intersect_index]

        sfunc_orthogonal = lambda i: sfunc_orthogonal_original(i) + distance_at_original_start - distance_at_wall

        eqReg.startInd = intersect_index

        d = eqReg.totalDistance()

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal, [0., 1.], [1., 0.])

        assert sfunc(0.) == tight_approx(0.)
        assert sfunc(n - 1.) == tight_approx(eqReg.totalDistance())

    def test_combineSfuncsPoloidalSpacingIntegrated(self, eqReg):
        # This test follows roughly the operations in
        # MeshRegion.distributePointsNonorthogonal, for 'poloidal_orthogonal_combined'
        # option
        n = len(eqReg)

        eqReg.poloidalSpacingParameters.polynomial_d_lower = .1
        eqReg.poloidalSpacingParameters.polynomial_d_upper = .1
        eqReg.poloidalSpacingParameters.nonorthogonal_range_lower = .2
        eqReg.poloidalSpacingParameters.nonorthogonal_range_upper = .2
        eqReg.poloidalSpacingParameters.N_norm = 40.

        sfunc_orthogonal_original = eqReg.contourSfunc()

        # as if lower_wall
        intersect_index = 2
        original_start = eqReg.startInd

        d = 1.5*3. / (n - 1.)
        eqReg.points.insert(intersect_index, Point2D(d, d))

        if original_start >= intersect_index:
            original_start += 1
        distance_at_original_start = eqReg.distance[original_start]

        distance_at_wall = eqReg.distance[intersect_index]

        sfunc_orthogonal = lambda i: sfunc_orthogonal_original(i) + distance_at_original_start - distance_at_wall

        eqReg.startInd = intersect_index

        d = eqReg.totalDistance()

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(0.) == tight_approx(0.)
        assert sfunc(n - 1.) == tight_approx(eqReg.totalDistance())
