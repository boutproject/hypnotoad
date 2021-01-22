# Copyright 2019 J.T. Omotani
#
# Contact John Omotani john.omotani@ukaea.uk
#
# This file is part of Hypnotoad 2.
#
# Hypnotoad 2 is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Hypnotoad 2 is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# Hypnotoad 2.  If not, see <http://www.gnu.org/licenses/>.

import numpy
import pytest
from copy import deepcopy
from hypnotoad.core.equilibrium import (
    calc_distance,
    find_intersections,
    Equilibrium,
    EquilibriumRegion,
    FineContour,
    Point2D,
    PsiContour,
)
from .utils_for_tests import tight_approx

PsiContour.user_options_factory = PsiContour.user_options_factory.add(
    finecontour_Nfine=1000, finecontour_atol=2.0e-8
)


class TestPoints:
    p0 = Point2D(1.0, 2.0)
    p1 = Point2D(3.0, 4.0)

    def test_add(self):
        p = self.p0 + self.p1
        assert p.R == tight_approx(4.0)
        assert p.Z == tight_approx(6.0)

    def test_sub(self):
        p = self.p0 - self.p1
        assert p.R == tight_approx(-2.0)
        assert p.Z == tight_approx(-2.0)

    def test_mul(self):
        p = self.p0 * 1.5
        assert p.R == tight_approx(1.5)
        assert p.Z == tight_approx(3.0)

    def test_rmul(self):
        p = 1.5 * self.p0
        assert p.R == tight_approx(1.5)
        assert p.Z == tight_approx(3.0)

    def test_div(self):
        p = self.p0 / 3.0
        assert p.R == tight_approx(1.0 / 3.0)
        assert p.Z == tight_approx(2.0 / 3.0)

    def test_iter(self):
        assert [x for x in self.p0] == tight_approx([1.0, 2.0])

    def test_repr(self):
        assert str(self.p0) == "Point2D(1.0,2.0)"

    def test_distance(self):
        assert calc_distance(self.p0, self.p1) == tight_approx(2.0 * numpy.sqrt(2.0))

    def test_as_ndarray(self):
        p = self.p0.as_ndarray()
        assert issubclass(type(p), numpy.ndarray)
        assert p == tight_approx(numpy.array([1.0, 2.0]))


def test_find_intersectionRR1():
    l1 = numpy.array([[-1.0, -0.1], [1.0, 0.1]])
    l2start = Point2D(-1.0, 0.1)
    l2end = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionRR2():
    l1 = numpy.array([[1.0, 0.1], [-1.0, -0.1]])
    l2start = Point2D(-1.0, 0.1)
    l2end = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionRR3():
    l1 = numpy.array([[-1.0, -0.1], [1.0, 0.1]])
    l2end = Point2D(-1.0, 0.1)
    l2start = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionRR4():
    l1 = numpy.array([[1.0, 0.1], [-1.0, -0.1]])
    l2end = Point2D(-1.0, 0.1)
    l2start = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionRRNone1():
    l1 = numpy.array([[2.0, -0.1], [4.0, 0.1]])
    l2start = Point2D(-1.0, 0.1)
    l2end = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionRRNone2():
    l1 = numpy.array([[4.0, 0.1], [2.0, -0.1]])
    l2start = Point2D(-1.0, 0.1)
    l2end = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionRRNone3():
    l1 = numpy.array([[2.0, -0.1], [4.0, 0.1]])
    l2end = Point2D(-1.0, 0.1)
    l2start = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionRRNone4():
    l1 = numpy.array([[4.0, 0.1], [2.0, -0.1]])
    l2end = Point2D(-1.0, 0.1)
    l2start = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionRZ1():
    l1 = numpy.array([[-1.0, -0.1], [1.0, 0.1]])
    l2start = Point2D(-0.1, 1.0)
    l2end = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionRZ2():
    l1 = numpy.array([[1.0, 0.1], [-1.0, -0.1]])
    l2start = Point2D(-0.1, 1.0)
    l2end = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionRZ3():
    l1 = numpy.array([[-1.0, -0.1], [1.0, 0.1]])
    l2end = Point2D(-0.1, 1.0)
    l2start = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionRZ4():
    l1 = numpy.array([[1.0, 0.1], [-1.0, -0.1]])
    l2end = Point2D(-0.1, 1.0)
    l2start = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionRZNone1():
    l1 = numpy.array([[2.0, -0.1], [4.0, 0.1]])
    l2start = Point2D(-0.1, 1.0)
    l2end = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionRZNone2():
    l1 = numpy.array([[4.0, 0.1], [2.0, -0.1]])
    l2start = Point2D(-0.1, 1.0)
    l2end = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionRZNone3():
    l1 = numpy.array([[2.0, -0.1], [4.0, 0.1]])
    l2end = Point2D(-0.1, 1.0)
    l2start = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionRZNone4():
    l1 = numpy.array([[4.0, 0.1], [2.0, -0.1]])
    l2end = Point2D(-0.1, 1.0)
    l2start = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionZR1():
    l1 = numpy.array([[-0.1, -1.0], [0.1, 1.0]])
    l2start = Point2D(-1.0, 0.1)
    l2end = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionZR2():
    l1 = numpy.array([[0.1, 1.0], [-0.1, -1.0]])
    l2start = Point2D(-1.0, 0.1)
    l2end = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionZR3():
    l1 = numpy.array([[-0.1, -1.0], [0.1, 1.0]])
    l2end = Point2D(-1.0, 0.1)
    l2start = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionZR4():
    l1 = numpy.array([[0.1, 1.0], [-0.1, -1.0]])
    l2end = Point2D(-1.0, 0.1)
    l2start = Point2D(1.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionZRNone1():
    l1 = numpy.array([[-0.1, -1.0], [0.1, 1.0]])
    l2start = Point2D(2.0, 0.1)
    l2end = Point2D(4.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionZRNone2():
    l1 = numpy.array([[0.1, 1.0], [-0.1, -1.0]])
    l2start = Point2D(2.0, 0.1)
    l2end = Point2D(4.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionZRNone3():
    l1 = numpy.array([[-0.1, -1.0], [0.1, 1.0]])
    l2end = Point2D(2.0, 0.1)
    l2start = Point2D(4.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionZRNone4():
    l1 = numpy.array([[0.1, 1.0], [-0.1, -1.0]])
    l2end = Point2D(2.0, 0.1)
    l2start = Point2D(4.0, -0.1)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionZZ1():
    l1 = numpy.array([[-0.1, -1.0], [0.1, 1.0]])
    l2start = Point2D(-0.1, 1.0)
    l2end = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionZZ2():
    l1 = numpy.array([[0.1, 1.0], [-0.1, -1.0]])
    l2start = Point2D(-0.1, 1.0)
    l2end = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionZZ3():
    l1 = numpy.array([[-0.1, -1.0], [0.1, 1.0]])
    l2end = Point2D(-0.1, 1.0)
    l2start = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionZZ4():
    l1 = numpy.array([[0.1, 1.0], [-0.1, -1.0]])
    l2end = Point2D(-0.1, 1.0)
    l2start = Point2D(0.1, -1.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect[0, :] == tight_approx([0.0, 0.0])


def test_find_intersectionZZNone1():
    l1 = numpy.array([[-0.1, -1.0], [0.1, 1.0]])
    l2start = Point2D(-0.1, 4.0)
    l2end = Point2D(0.1, 2.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionZZNone2():
    l1 = numpy.array([[0.1, 1.0], [-0.1, -1.0]])
    l2start = Point2D(-0.1, 4.0)
    l2end = Point2D(0.1, 2.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionZZNone3():
    l1 = numpy.array([[-0.1, -1.0], [0.1, 1.0]])
    l2end = Point2D(-0.1, 4.0)
    l2start = Point2D(0.1, 2.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


def test_find_intersectionZZNone4():
    l1 = numpy.array([[0.1, 1.0], [-0.1, -1.0]])
    l2end = Point2D(-0.1, 4.0)
    l2start = Point2D(0.1, 2.0)
    intersect = find_intersections(l1, l2start, l2end)
    assert intersect is None


class TestContour:
    @pytest.fixture
    def testcontour(self):
        # make a circle, not centred on origin
        class returnObject:
            npoints = 23
            r = 1.0
            R0 = 0.2
            Z0 = 0.3
            theta = numpy.linspace(0.0, numpy.pi, npoints)

            def __init__(self):
                def psifunc(R, Z):
                    return (R - self.R0) ** 2 + (Z - self.Z0) ** 2

                self.R = self.R0 + self.r * numpy.cos(self.theta)
                self.Z = self.Z0 + self.r * numpy.sin(self.theta)

                psi_xpoint = psifunc(self.R[0], self.Z[0])

                self.c = PsiContour(
                    points=[Point2D(R, Z) for R, Z in zip(self.R, self.Z)],
                    psi=psifunc,
                    psival=psi_xpoint,
                    settings={"refine_width": 1.0e-3, "refine_methods": "line"},
                )

        return returnObject()

    def test_distance(self, testcontour):
        segment_length = testcontour.r * numpy.pi / (testcontour.npoints - 1)
        assert testcontour.c.distance == pytest.approx(
            segment_length * numpy.arange(testcontour.npoints), abs=1.0e-4
        )

    def test_iter(self, testcontour):
        clist = list(testcontour.c)

        for i, item in enumerate(clist):
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
        assert c.distance[-1] == pytest.approx(expected_distance, abs=4.0e-5)

    def test_reverse(self, testcontour):
        c = testcontour.c
        orig = deepcopy(c)

        c.reverse()

        n = len(orig)
        total_d = orig.distance[-1]
        for i in range(n):
            assert orig[n - 1 - i].R == tight_approx(c[i].R)
            assert orig[n - 1 - i].Z == tight_approx(c[i].Z)
            assert total_d - orig.distance[n - 1 - i] == tight_approx(c.distance[i])

    def test_refine(self, testcontour):
        # PsiContour.refine just calls PsiContour.getRefined, so this tests both

        c = testcontour.c
        c.psival = 0.7
        # c does not start close to a contour of psi_func, so need to use a large width
        c.refine(width=2.0, atol=1.0e-13)
        for p in c:
            assert c.psi(p.R, p.Z) == tight_approx(0.7)

    def test_coarseInterp(self, testcontour):
        c = testcontour.c
        c.startInd = 2
        f, distance_estimate = c._coarseInterp()
        pstart = f(0.0)
        pend = f(distance_estimate[-1])
        assert pstart.R == pytest.approx(c[2].R, abs=1.0e-9)
        assert pstart.Z == pytest.approx(c[2].Z, abs=1.0e-5)
        assert pend.R == pytest.approx(testcontour.R0 - testcontour.r, abs=1.0e-9)
        assert pend.Z == pytest.approx(testcontour.Z0, abs=1.0e-5)

    def test_interpFunction(self, testcontour):
        f = testcontour.c.interpFunction()
        pstart = f(0.0)
        pend = f(numpy.pi * testcontour.r)
        assert pstart.R == pytest.approx(testcontour.R0 + testcontour.r, abs=1.0e-9)
        assert pstart.Z == pytest.approx(testcontour.Z0, abs=1.0e-5)
        assert pend.R == pytest.approx(testcontour.R0 - testcontour.r, abs=1.0e-9)
        assert pend.Z == pytest.approx(testcontour.Z0, abs=1.0e-5)

    def test_getRegridded(self, testcontour):
        orig = testcontour.c
        r = testcontour.r

        newNpoints = 97

        def sfunc_true(i):
            return numpy.sqrt(i / (newNpoints - 1)) * numpy.pi * r

        def sfunc(i):
            return numpy.sqrt(i / (newNpoints - 1)) * orig.distance[-1]

        newTheta = sfunc_true(numpy.arange(newNpoints)) / r
        newR = testcontour.R0 + r * numpy.cos(newTheta)
        newZ = testcontour.Z0 + r * numpy.sin(newTheta)

        new = orig.getRegridded(newNpoints, sfunc=sfunc, width=1.0e-3)

        assert [p.R for p in new] == pytest.approx(newR, abs=4.0e-4)
        assert [p.Z for p in new] == pytest.approx(newZ, abs=4.0e-4)

    def test_getRegridded_extend(self, testcontour):
        c = testcontour.c
        orig = c.newContourFromSelf()

        new = c.getRegridded(
            testcontour.npoints, width=0.1, extend_lower=1, extend_upper=2
        )

        assert numpy.array([[*p] for p in new[1:-2]]) == pytest.approx(
            numpy.array([[*p] for p in orig]), abs=1.0e-8
        )

        # test the extend_lower
        assert [*new[0]] == pytest.approx(
            [orig[1].R, 2.0 * testcontour.Z0 - orig[1].Z], abs=1.0e-10
        )

        # test the extend_upper
        assert [*new[-2]] == pytest.approx(
            [orig[-2].R, 2.0 * testcontour.Z0 - orig[-2].Z], abs=1.0e-10
        )
        assert [*new[-1]] == pytest.approx(
            [orig[-3].R, 2.0 * testcontour.Z0 - orig[-3].Z], abs=1.0e-10
        )

    def test_contourSfunc(self, testcontour):
        c = testcontour.c
        c.startInd = 2
        c.endInd = len(c) - 2
        n = c.endInd - c.startInd + 1
        c.extend_lower = 2
        c.extend_upper = 2

        f = c.contourSfunc()

        indices = numpy.arange(n, dtype=float)
        assert f(indices) == tight_approx(
            [d - c.distance[c.startInd] for d in c.distance[c.startInd : c.endInd + 1]]
        )
        assert f(-1.0) == tight_approx(0.0)
        assert f(n + 1.0) == tight_approx(c.distance[c.endInd] - c.distance[c.startInd])

    def test_contourSfunc_list(self, testcontour):
        c1 = testcontour.c
        c1.startInd = 2
        c1.endInd = len(c1) - 2
        n1 = c1.endInd - c1.startInd
        c1.extend_lower = 2
        c1.extend_upper = 2
        npoints1 = len(c1)

        r = testcontour.r
        R0 = testcontour.R0
        Z0 = testcontour.Z0
        theta = numpy.linspace(0.0, 1.5 * numpy.pi, 37)

        R = R0 + r * numpy.cos(theta)
        Z = Z0 + r * numpy.sin(theta)

        c2 = PsiContour(
            points=[Point2D(R, Z) for R, Z in zip(R, Z)],
            psi=c1.psi,
            psival=c1.psival,
            settings=dict(c1.user_options),
        )
        c2.startInd = 2
        c2.endInd = len(c2) - 1
        n2 = c2.endInd - c2.startInd
        c2.extend_lower = 2
        c2.extend_upper = 1
        npoints2 = len(c2)

        c_list = [c1, c2]
        sfunc_list = []

        # This does NOT work as desired - when the function added to 'sfunc_list' is
        # evaluated, it takes 'sfunc_orig' from the enclosing scope, which means it uses
        # the last value from the loop instead of the value when 'sfunc' was added to
        # 'sfunc_list'
        for c in c_list:
            sfunc_orig = c.contourSfunc()

            sfunc_list.append(lambda i: sfunc_orig(i) - 3.0)

        # notice we check that the first test *fails*
        assert not sfunc_list[0](float(n1)) == pytest.approx(
            n1 / (npoints1 - 1.0) * numpy.pi * r - 3.0, abs=1.0e-6
        )
        assert sfunc_list[1](float(n2)) == pytest.approx(
            n2 / (npoints2 - 1.0) * 1.5 * numpy.pi * r - 3.0, abs=4.0e-6
        )

        sfunc_list2 = []

        # This version does work, because when the lambda is evaluated it uses
        # 'sfunc_orig' from the scope of 'shift_sfunc' in which it was created.
        def shift_sfunc(c):
            sfunc_orig = c.contourSfunc()
            return lambda i: sfunc_orig(i) - 3.0

        for c in c_list:
            sfunc_list2.append(shift_sfunc(c))

        assert sfunc_list2[0](float(n1)) == pytest.approx(
            n1 / (npoints1 - 1.0) * numpy.pi * r - 3.0, abs=1.0e-6
        )
        assert sfunc_list2[1](float(n2)) == pytest.approx(
            n2 / (npoints2 - 1.0) * 1.5 * numpy.pi * r - 3.0, abs=4.0e-6
        )

    def test_interpSSperp(self, testcontour):
        c = testcontour.c

        # Make c.startInd > 0
        c.insert(0, Point2D(c[1].R, 2.0 * c[0].Z - c[1].Z))

        # 'vec' argument is in Z-direction, so 's_perp' is displacement in R-direction
        sfunc, s_perp_total = c.interpSSperp([0.0, 1.0])
        assert sfunc(0.0) == tight_approx(0.0)
        assert sfunc(1.0) == pytest.approx(numpy.pi / 2.0, abs=1.0e-6)
        assert sfunc(2.0) == pytest.approx(numpy.pi, abs=2.0e-6)
        assert s_perp_total == tight_approx(2.0)

    def test_FineContour(self, testcontour):
        testcontour.c.refine_width = 1.0e-2

        fc = FineContour(testcontour.c, dict(testcontour.c.user_options))

        assert fc.totalDistance() == pytest.approx(numpy.pi, abs=1.0e-5)

        interpFunc = fc.interpFunction()
        r = testcontour.r
        for theta in numpy.linspace(0.0, numpy.pi, 17):
            p = interpFunc(r * theta)
            assert p.R == pytest.approx(
                testcontour.R0 + r * numpy.cos(theta), abs=1.0e-4
            )
            assert p.Z == pytest.approx(
                testcontour.Z0 + r * numpy.sin(theta), abs=1.0e-4
            )


class ThisEquilibrium(Equilibrium):
    def __init__(self, settings=None):
        if settings is None:
            settings = {}
        self.user_options = Equilibrium.user_options_factory.add(
            refine_width=1.0e-5, refine_atol=2.0e-8
        ).create(settings)

        super().__init__({})


class TestEquilibrium:
    @pytest.fixture
    def eq(self):
        eq = ThisEquilibrium()
        eq.wall = [
            Point2D(-1.0, -1.0),
            Point2D(1.0, -1.0),
            Point2D(1.0, 1.0),
            Point2D(-1.0, 1.0),
        ]
        return eq

    def test_make1dGrid(self, eq):
        n = 4

        def f(i):
            return i ** 2

        r = eq.make1dGrid(n, f)
        assert r == tight_approx([0.0, 0.5, 1.0, 2.5, 4.0, 6.5, 9.0, 12.5, 16.0])

    def test_wallIntersection(self, eq):
        intersect = eq.wallIntersection(Point2D(0.0, 0.0), Point2D(2.0, 0.0))
        assert intersect.R == tight_approx(1.0)
        assert intersect.Z == tight_approx(0.0)

        # Line goes exactly through a corner of the wall, so intersects with two wall
        # segments
        intersect = eq.wallIntersection(Point2D(0.0, 0.0), Point2D(2.0, 2.0))
        assert intersect.R == tight_approx(1.0)
        assert intersect.Z == tight_approx(1.0)

    @pytest.mark.parametrize(
        ["grad_lower", "lower", "upper"], [[0.2, 0.4, 2.0], [-0.2, 2.0, 0.4]]
    )
    def test_getPolynomialGridFuncGradLowerDecreasing(
        self, eq, grad_lower, lower, upper
    ):
        # Test getPolynomialGridFunc() with grad_lower set so that it needs to decrease
        # the average spacing
        N = 10.0
        f = eq.getPolynomialGridFunc(N, lower, upper, grad_lower=grad_lower)
        # f(0) = lower
        assert f(0.0) == tight_approx(lower)
        # f(N) = upper
        assert f(N) == pytest.approx(upper, rel=8.0e-12, abs=1.0e-13)
        # for i<<1, df/di = grad_lower
        itest = 1.0e-3
        assert (f(itest) - f(0.0)) / itest == pytest.approx(grad_lower, abs=1.0e-5)
        # for i<<1, d2f/di2 = 0
        assert (f(0.0) - 2.0 * f(itest) + f(2.0 * itest)) / itest ** 2 == pytest.approx(
            0.0, abs=1.0e-5
        )

    @pytest.mark.parametrize(
        ["grad_lower", "lower", "upper"], [[0.02, 0.4, 2.0], [-0.02, 2.0, 0.4]]
    )
    def test_getPolynomialGridFuncGradLowerIncreasing(
        self, eq, grad_lower, lower, upper
    ):
        # Test getPolynomialGridFunc() with grad_lower set so that it needs to increase
        # the average spacing
        N = 10.0
        f = eq.getPolynomialGridFunc(N, lower, upper, grad_lower=grad_lower)
        # f(0) = lower
        assert f(0.0) == tight_approx(lower)
        # f(N) = upper
        assert f(N) == tight_approx(upper)
        # for i<<1, df/di = grad_lower
        itest = 1.0e-3
        assert (f(itest) - f(0.0)) / itest == pytest.approx(grad_lower, abs=1.0e-5)
        # for i<<1, d2f/di2 = 0
        assert (f(0.0) - 2.0 * f(itest) + f(2.0 * itest)) / itest ** 2 == pytest.approx(
            0.0, abs=1.0e-5
        )

    @pytest.mark.parametrize(
        ["grad_upper", "lower", "upper"], [[0.5, 0.4, 2.0], [-0.5, 2.0, 0.4]]
    )
    def test_getPolynomialGridFuncGradUpperDecreasing(
        self, eq, grad_upper, lower, upper
    ):
        # Test getPolynomialGridFunc() with grad_upper set so that it needs to decrease
        # the average spacing
        N = 10.0
        f = eq.getPolynomialGridFunc(N, lower, upper, grad_upper=grad_upper)
        # f(0) = lower
        assert f(0.0) == tight_approx(lower)
        # f(N) = upper
        assert f(N) == tight_approx(upper)
        # for N-i<<1, df/di = grad_upper
        itest = 1.0e-4
        assert (f(N) - f(N - itest)) / itest == pytest.approx(grad_upper, abs=1.0e-5)
        # for N-i<<1, d2f/di2 = 0
        assert (
            f(N) - 2.0 * f(N - itest) + f(N - 2.0 * itest)
        ) / itest ** 2 == pytest.approx(0.0, abs=1.0e-5)

    @pytest.mark.parametrize(
        ["grad_upper", "lower", "upper"], [[0.1, 0.4, 2.0], [-0.1, 2.0, 0.4]]
    )
    def test_getPolynomialGridFuncGradUpperIncreasing(
        self, eq, grad_upper, lower, upper
    ):
        # Test getPolynomialGridFunc() with grad_upper set so that it needs to increase
        # the average spacing
        N = 10.0
        f = eq.getPolynomialGridFunc(N, lower, upper, grad_upper=grad_upper)
        # f(0) = lower
        assert f(0.0) == tight_approx(lower)
        # f(N) = upper
        assert f(N) == tight_approx(upper)
        # for N-i<<1, df/di = grad_upper
        itest = 1.0e-3
        assert (f(N) - f(N - itest)) / itest == pytest.approx(grad_upper, abs=1.0e-5)
        # for N-i<<1, d2f/di2 = 0
        assert (
            f(N) - 2.0 * f(N - itest) + f(N - 2.0 * itest)
        ) / itest ** 2 == pytest.approx(0.0, abs=1.0e-5)

    @pytest.mark.parametrize(
        ["grad_lower", "grad_upper", "lower", "upper"],
        [[0.4, 0.2, 0.4, 2.0], [-0.4, -0.2, 2.0, 0.4]],
    )
    def test_getPolynomialGridFuncGradBothDecreasing(
        self, eq, grad_lower, grad_upper, lower, upper
    ):
        # Test getPolynomialGridFunc() with grad_lower and grad_upper set so that it
        # needs to decrease the average spacing
        N = 10.0
        f = eq.getPolynomialGridFunc(
            N, lower, upper, grad_lower=grad_lower, grad_upper=grad_upper
        )
        # f(0) = lower
        assert f(0.0) == tight_approx(lower)
        # f(N) = upper
        assert f(N) == pytest.approx(upper, rel=1.0e-10, abs=1.0e-13)
        # for i<<1, df/di = grad_lower
        itest = 1.0e-5
        assert (f(itest) - f(0.0)) / itest == pytest.approx(grad_lower, abs=1.0e-5)
        # for i<<1, d2f/di2 = 0
        assert (f(0.0) - 2.0 * f(itest) + f(2.0 * itest)) / itest ** 2 == pytest.approx(
            0.0, abs=1.0e-5
        )
        # for N-i<<1, df/di = grad_upper
        assert (f(N) - f(N - itest)) / itest == pytest.approx(grad_upper, abs=1.0e-5)
        # for N-i<<1, d2f/di2 = 0
        assert (
            f(N) - 2.0 * f(N - itest) + f(N - 2.0 * itest)
        ) / itest ** 2 == pytest.approx(0.0, abs=1.0e-5)

    @pytest.mark.parametrize(
        ["grad_lower", "grad_upper", "lower", "upper"],
        [[0.2, 0.1, 0.4, 2.0], [-0.2, -0.1, 2.0, 0.4]],
    )
    def test_getPolynomialGridFuncGradBothIncreasing(
        self, eq, grad_lower, grad_upper, lower, upper
    ):
        # Test getPolynomialGridFunc() with grad_lower and grad_upper set so that it
        # needs to increase the average spacing
        grad_lower = 0.2
        grad_upper = 0.1
        lower = 0.4
        upper = 2.0
        N = 10.0
        f = eq.getPolynomialGridFunc(
            N, lower, upper, grad_lower=grad_lower, grad_upper=grad_upper
        )
        # f(0) = lower
        assert f(0.0) == tight_approx(lower)
        # f(N) = upper
        assert f(N) == tight_approx(upper)
        # for i<<1, df/di = grad_lower
        itest = 5.0e-4
        assert (f(itest) - f(0.0)) / itest == pytest.approx(grad_lower, abs=1.0e-5)
        # for i<<1, d2f/di2 = 0
        assert (f(0.0) - 2.0 * f(itest) + f(2.0 * itest)) / itest ** 2 == pytest.approx(
            0.0, abs=1.0e-5
        )
        # for N-i<<1, df/di = grad_upper
        assert (f(N) - f(N - itest)) / itest == pytest.approx(grad_upper, abs=1.0e-5)
        # for N-i<<1, d2f/di2 = 0
        assert (
            f(N) - 2.0 * f(N - itest) + f(N - 2.0 * itest)
        ) / itest ** 2 == pytest.approx(0.0, abs=1.0e-5)


class TestEquilibriumRegion:
    @pytest.fixture
    def eqReg(self):
        equilib = ThisEquilibrium(settings={"y_boundary_guards": 1})
        equilib.psi = lambda R, Z: R - Z
        n = 11.0
        points = [
            Point2D(i * 3.0 / (n - 1.0), i * 3.0 / (n - 1.0)) for i in numpy.arange(n)
        ]
        eqReg = EquilibriumRegion(
            equilibrium=equilib,
            name="",
            nSegments=1,
            nx=[1],
            ny=5,
            kind="wall.wall",
            ny_total=5,
            points=points,
            psival=0.0,
        )
        return eqReg

    def test_getMonotonicPoloidalDistanceFunc(self, eqReg):
        d_lower = 0.02
        d_upper = 0.01
        L = 2.0
        N = 10.0
        N_norm = 40.0
        f = eqReg.getMonotonicPoloidalDistanceFunc(
            L, N, N_norm, d_lower=d_lower, d_upper=d_upper
        )
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f = d_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(d_lower * itest / N_norm, abs=1.0e-5)
        # for N-i<<1, f = L - d_upper*i/N_norm
        assert f(N - itest) == pytest.approx(L - d_upper * itest / N_norm, abs=1.0e-5)

    def test_getSqrtPoloidalDistanceFuncLinear(self, eqReg):
        L = 2.0
        N = 10.0
        N_norm = 1
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm)
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(10.0) == tight_approx(2.0)
        # f(i) = i/N*L
        assert f(3.0) == tight_approx(0.6)

        # test gradients at upper and lower bounds
        dfdi = L / N
        # i=0, interior
        itest = 1.0e-3
        assert (f(itest) - f(0.0)) / itest == tight_approx(dfdi)
        # i=0, extropolating
        itest = -1.0e-3
        assert (f(itest) - f(0.0)) / itest == tight_approx(dfdi)
        # i=N, interior
        itest = N - 1.0e-3
        assert (f(N) - f(itest)) / (N - itest) == tight_approx(dfdi)
        # i=N, extrapolating
        itest = N + 1.0e-3
        assert (f(N) - f(itest)) / (N - itest) == tight_approx(dfdi)

    def test_getSqrtPoloidalDistanceFuncBLower(self, eqReg):
        b_lower = 0.01
        L = 2.0
        N = 10.0
        N_norm = 40.0
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm, b_lower=b_lower)
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(b_lower * itest / N_norm, abs=1.0e-5)

        # test gradient at upper bound
        ##############################
        # Copied from `b_upper is None` case of getSqrtPoloidalDistanceFunc():
        # f(iN) = c + d*iN + e*iN**2
        # df/diN(N/N_norm) = d + 2*e*N/N_norm
        # c = 0
        d = b_lower
        e = (L - d * N / N_norm) / (N / N_norm) ** 2
        dfdi = (d + 2 * e * N / N_norm) / N_norm

        # test in interior
        delta = 1.0e-4
        itest = N - delta
        assert (f(N) - f(itest)) / delta == pytest.approx(dfdi, abs=1.0e-5)
        # test extrapolating
        itest = N + delta
        assert (f(itest) - f(N)) / delta == pytest.approx(dfdi, abs=1.0e-5)

    def test_getSqrtPoloidalDistanceFuncBUpper(self, eqReg):
        b_upper = 0.01
        L = 2.0
        N = 10.0
        N_norm = 40.0
        f = eqReg.getSqrtPoloidalDistanceFunc(L, N, N_norm, b_upper=b_upper)
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for (N-i)<<1, f ~ L - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(L - b_upper * (N - itest) / N_norm, abs=1.0e-5)

        # test gradient at lower bound
        ##############################
        # Copied from `b_lower is None` case of getSqrtPoloidalDistanceFunc():
        # f(iN) = c + d*iN + e*iN**2
        # df/diN(0) = d
        # c = 0
        e = (b_upper * N / N_norm - L) / (N / N_norm) ** 2
        d = b_upper - 2 * e * N / N_norm
        dfdi = d / N_norm

        # test in interior
        itest = 1.0e-4
        assert (f(itest) - f(0.0)) / itest == pytest.approx(dfdi, abs=1.0e-5)
        # test extrapolating
        itest = -1.0e-4
        assert (f(itest) - f(0.0)) / itest == pytest.approx(dfdi, abs=1.0e-5)

    def test_getSqrtPoloidalDistanceFuncBBoth(self, eqReg):
        b_lower = 0.1
        b_upper = 0.2
        L = 2.0
        N = 10.0
        N_norm = 40.0
        f = eqReg.getSqrtPoloidalDistanceFunc(
            L, N, N_norm, b_lower=b_lower, b_upper=b_upper
        )
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(b_lower * itest / N_norm, abs=1.0e-5)
        # Check we can extrapolate at lower end
        itest = -0.01
        assert f(itest) == pytest.approx(b_lower * itest / N_norm, abs=1.0e-5)
        # for (N-i)<<1, f ~ L - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(L - b_upper * (N - itest) / N_norm, abs=1.0e-5)
        # Check we can extrapolate at upper end
        itest = N + 0.01
        assert f(itest) == pytest.approx(L - b_upper * (N - itest) / N_norm, abs=1.0e-5)

    def test_getSqrtPoloidalDistanceFuncBothLower(self, eqReg):
        b_lower = 0.01
        a_lower = 0.05
        L = 2.0
        N = 10.0
        N_norm = 2.0
        f = eqReg.getSqrtPoloidalDistanceFunc(
            L, N, N_norm, b_lower=b_lower, a_lower=a_lower
        )
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ 2*a_lower*sqrt(i/N_norm) + b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(
            2.0 * a_lower * numpy.sqrt(itest / N_norm) + b_lower * itest / N_norm,
            abs=1.0e-5,
        )

        # test gradient at upper bound
        ##############################
        # Copied from `b_upper is None` case of getSqrtPoloidalDistanceFunc():
        # f(iN) = a*sqrt(iN) + c + d*iN + e*iN**2
        # dfdiN(N/N_norm) = a/2/sqrt(N/N_norm) + d + 2*e*N/N_norm
        # c = 0
        a = 2.0 * a_lower
        d = b_lower
        e = (L - a * numpy.sqrt(N / N_norm) - d * N / N_norm) / (N / N_norm) ** 2
        dfdi = (a / 2.0 / numpy.sqrt(N / N_norm) + d + 2.0 * e * N / N_norm) / N_norm

        # test in interior
        delta = 1.0e-4
        itest = N - delta
        assert (f(N) - f(itest)) / delta == pytest.approx(dfdi, abs=1.0e-5)
        # test extrapolating
        itest = N + delta
        assert (f(itest) - f(N)) / delta == pytest.approx(dfdi, abs=1.0e-5)

    def test_getSqrtPoloidalDistanceFuncBothUpper(self, eqReg):
        b_upper = 0.01
        a_upper = 0.05
        L = 2.0
        N = 10.0
        N_norm = 1
        f = eqReg.getSqrtPoloidalDistanceFunc(
            L, N, N_norm, b_upper=b_upper, a_upper=a_upper
        )
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for (N-i)<<1, f ~ L - 2*a_upper*sqrt((N-i)/N_norm) - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(
            L
            - 2.0 * a_upper * numpy.sqrt((N - itest) / N_norm)
            - b_upper * (N - itest) / N_norm,
            abs=1.0e-5,
        )

        # test gradient at lower bound
        ##############################
        # Copied from `b_lower is None` case of getSqrtPoloidalDistanceFunc():
        # f(iN) = -b*sqrt(N/N_norm-iN) + c + d*iN + e*iN**2
        # df/diN(0) = b/2/sqrt(N/N_norm) + d
        b = 2.0 * a_upper
        c = b * numpy.sqrt(N / N_norm)
        e = (c + b_upper * N / N_norm - L) / (N / N_norm) ** 2
        d = b_upper - 2 * e * N / N_norm
        dfdi = (b / 2.0 / numpy.sqrt(N / N_norm) + d) / N_norm

        # test in interior
        itest = 1.0e-4
        assert (f(itest) - f(0.0)) / itest == pytest.approx(dfdi, abs=1.0e-5)
        # test extrapolating
        itest = -1.0e-4
        assert (f(itest) - f(0.0)) / itest == pytest.approx(dfdi, abs=1.0e-5)

    def test_getSqrtPoloidalDistanceFuncALowerBBoth(self, eqReg):
        b_lower = 0.01
        a_lower = 0.05
        b_upper = 0.2
        L = 2.0
        L = 2.0
        N = 10.0
        N_norm = 1
        f = eqReg.getSqrtPoloidalDistanceFunc(
            L,
            N,
            N_norm,
            b_lower=b_lower,
            a_lower=a_lower,
            b_upper=b_upper,
        )
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ 2*a_lower*sqrt(i/N_norm) + b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(
            2.0 * a_lower * numpy.sqrt(itest / N_norm) + b_lower * itest / N_norm,
            abs=1.0e-5,
        )
        # for (N-i)<<1, f ~ L - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(
            L - b_upper * (N - itest) / N_norm,
            abs=1.0e-5,
        )
        # Check we can also extrapolate past the end
        itest = N + 0.01
        assert f(itest) == pytest.approx(
            L - b_upper * (N - itest) / N_norm,
            abs=1.0e-5,
        )

    def test_getSqrtPoloidalDistanceFuncAUpperBBoth(self, eqReg):
        b_lower = 0.01
        b_upper = 0.2
        a_upper = 0.07
        L = 2.0
        L = 2.0
        N = 10.0
        N_norm = 1
        f = eqReg.getSqrtPoloidalDistanceFunc(
            L,
            N,
            N_norm,
            b_lower=b_lower,
            b_upper=b_upper,
            a_upper=a_upper,
        )
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(
            b_lower * itest / N_norm,
            abs=1.0e-5,
        )
        # Check we can also extrapolate past the end
        itest = -0.01
        assert f(itest) == pytest.approx(
            b_lower * itest / N_norm,
            abs=1.0e-5,
        )
        # for (N-i)<<1, f ~ L - 2*a_upper*sqrt((N-i)/N_norm) - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(
            L
            - 2.0 * a_upper * numpy.sqrt((N - itest) / N_norm)
            - b_upper * (N - itest) / N_norm,
            abs=1.0e-5,
        )

    def test_getSqrtPoloidalDistanceFuncBothBoth(self, eqReg):
        b_lower = 0.01
        a_lower = 0.05
        b_upper = 0.2
        a_upper = 0.07
        L = 2.0
        L = 2.0
        N = 10.0
        N_norm = 1
        f = eqReg.getSqrtPoloidalDistanceFunc(
            L,
            N,
            N_norm,
            b_lower=b_lower,
            a_lower=a_lower,
            b_upper=b_upper,
            a_upper=a_upper,
        )
        # f(0) = 0
        assert f(0.0) == tight_approx(0.0)
        # f(N) = L
        assert f(N) == tight_approx(L)
        # for i<<1, f ~ 2*a_lower*sqrt(i/N_norm) + b_lower*i/N_norm
        itest = 0.01
        assert f(itest) == pytest.approx(
            2.0 * a_lower * numpy.sqrt(itest / N_norm) + b_lower * itest / N_norm,
            abs=1.0e-5,
        )
        # for (N-i)<<1, f ~ L - 2*a_upper*sqrt((N-i)/N_norm) - b_upper*(N-i)/N_norm
        itest = N - 0.01
        assert f(itest) == pytest.approx(
            L
            - 2.0 * a_upper * numpy.sqrt((N - itest) / N_norm)
            - b_upper * (N - itest) / N_norm,
            abs=1.0e-5,
        )

    def test_combineSfuncsPoloidalSpacing(self, eqReg):
        n = len(eqReg)
        L = eqReg.totalDistance()

        eqReg.resetNonorthogonalOptions(
            {"nonorthogonal_target_all_poloidal_spacing_length": L}
        )
        eqReg.ny_total = n - 1
        eqReg.name = "inner"

        def sfunc_orthogonal(i):
            return i / (n - 1.0) * L

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(0.0) == tight_approx(0.0)
        assert sfunc((n - 1.0) / 3.0) == tight_approx(L / 3.0)
        assert sfunc(n - 1.0) == tight_approx(L)

    def test_combineSfuncsPoloidalSpacingRangeLower(self, eqReg):
        n = len(eqReg)
        L = eqReg.totalDistance()

        new_settings = {
            "nonorthogonal_target_all_poloidal_spacing_length": L * 40.0 / (n - 1),
            "nonorthogonal_target_all_poloidal_spacing_range": 0.1,
        }
        eqReg.resetNonorthogonalOptions(new_settings)
        eqReg.ny_total = 40
        eqReg.name = "outer"

        def sfunc_orthogonal(i):
            return numpy.piecewise(
                i, [i < 0.0], [100.0, lambda i: numpy.sqrt(i / (n - 1.0)) * L]
            )

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(-1.0) == tight_approx(-1.0 / (n - 1.0) * L)
        assert sfunc(0.0) == tight_approx(0.0)
        itest = 0.0001
        assert sfunc(itest) == pytest.approx(
            0.00001 * 3.0 * numpy.sqrt(2.0), abs=1.0e-11
        )
        assert sfunc(n - 1.0) == tight_approx(L)

    def test_combineSfuncsPoloidalSpacingRangeUpper(self, eqReg):
        n = len(eqReg)
        L = eqReg.totalDistance()

        new_settings = {
            "nonorthogonal_target_all_poloidal_spacing_length": L * 40.0 / (n - 1),
            "nonorthogonal_target_all_poloidal_spacing_range": 0.1,
        }
        eqReg.resetNonorthogonalOptions(new_settings)
        eqReg.ny_total = 40
        eqReg.name = "inner_upper"

        def sfunc_orthogonal(i):
            return numpy.piecewise(
                i, [i > n - 1.0], [100.0, lambda i: numpy.sqrt(i / (n - 1.0)) * L]
            )

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(0.0) == tight_approx(0.0)
        itest = n - 1.0 - 0.0001
        assert L - sfunc(itest) == pytest.approx(
            0.00001 * 3.0 * numpy.sqrt(2.0), abs=1.0e-11
        )
        assert sfunc(n - 1.0) == tight_approx(L)
        assert sfunc(float(n)) == tight_approx(n / (n - 1.0) * L)

    def test_combineSfuncsPoloidalSpacingRangeBoth(self, eqReg):
        n = len(eqReg)
        L = eqReg.totalDistance()

        new_settings = {
            "nonorthogonal_target_all_poloidal_spacing_length": L * 40.0 / (n - 1),
            "nonorthogonal_target_all_poloidal_spacing_range": 0.1,
        }
        eqReg.resetNonorthogonalOptions(new_settings)
        eqReg.ny_total = 40
        eqReg.name = "outer_upper"

        def sfunc_orthogonal(i):
            return numpy.piecewise(
                i,
                [i < 0.0, i > n - 1.0],
                [-100.0, 100.0, lambda i: numpy.sqrt(i / (n - 1.0)) * L],
            )

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(-1.0) == tight_approx(-1.0 / (n - 1.0) * L)
        assert sfunc(0.0) == tight_approx(0.0)
        itest = 1.0e-6
        assert sfunc(itest) == pytest.approx(3.0e-7 * numpy.sqrt(2.0), rel=1.0e-7)
        itest = n - 1.0 - 1.0e-6
        assert L - sfunc(itest) == pytest.approx(3.0e-7 * numpy.sqrt(2.0), rel=1.0e-7)
        assert sfunc(n - 1.0) == tight_approx(L)
        assert sfunc(float(n)) == tight_approx(n / (n - 1.0) * L)

    def test_combineSfuncsPerpSpacingIntegrated(self, eqReg):
        # This test follows roughly the operations in
        # MeshRegion.distributePointsNonorthogonal for the default 'combined' option.
        n = len(eqReg)

        new_settings = {
            "nonorthogonal_target_all_poloidal_spacing_length": 0.1,
            "nonorthogonal_target_all_poloidal_spacing_range": 0.3,
        }
        eqReg.resetNonorthogonalOptions(new_settings)
        eqReg.ny_total = 40
        eqReg.sin_angle_at_start = 1.0
        eqReg.sin_angle_at_end = 1.0
        eqReg.name = "inner_lower"

        sfunc_orthogonal_original = eqReg.contourSfunc()

        # as if lower_wall
        intersect_index = 2
        original_start = eqReg.startInd
        distance_at_original_start = eqReg.distance[original_start]

        d = 1.5 * 3.0 / (n - 1.0)
        eqReg.points.insert(intersect_index, Point2D(d, d))

        if original_start >= intersect_index:
            original_start += 1

        distance_at_wall = eqReg.distance[intersect_index]

        sfunc_orthogonal = (
            lambda i: sfunc_orthogonal_original(i)
            + distance_at_original_start
            - distance_at_wall
        )

        eqReg.startInd = intersect_index
        eqReg.extend_lower = intersect_index
        eqReg.extend_upper = 1

        d = eqReg.totalDistance()

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal, [0.0, 1.0], [1.0, 0.0])

        assert sfunc(0.0) == tight_approx(0.0)
        assert sfunc(n - 1.0) == tight_approx(eqReg.totalDistance())

    def test_combineSfuncsPoloidalSpacingIntegrated(self, eqReg):
        # This test follows roughly the operations in
        # MeshRegion.distributePointsNonorthogonal, for 'poloidal_orthogonal_combined'
        # option
        n = len(eqReg)

        new_settings = {
            "nonorthogonal_target_all_poloidal_spacing_length": 0.1,
            "nonorthogonal_target_all_poloidal_spacing_range": 0.2,
        }
        eqReg.resetNonorthogonalOptions(new_settings)
        eqReg.ny_total = 40
        eqReg.name = "outer_lower"

        sfunc_orthogonal_original = eqReg.contourSfunc()

        # as if lower_wall
        intersect_index = 2
        original_start = eqReg.startInd

        d = 1.5 * 3.0 / (n - 1.0)
        eqReg.points.insert(intersect_index, Point2D(d, d))

        if original_start >= intersect_index:
            original_start += 1
        distance_at_original_start = eqReg.distance[original_start]

        distance_at_wall = eqReg.distance[intersect_index]

        sfunc_orthogonal = (
            lambda i: sfunc_orthogonal_original(i)
            + distance_at_original_start
            - distance_at_wall
        )

        eqReg.startInd = intersect_index
        eqReg.extend_lower = intersect_index
        eqReg.extend_upper = 1

        d = eqReg.totalDistance()

        sfunc = eqReg.combineSfuncs(eqReg, sfunc_orthogonal)

        assert sfunc(0.0) == tight_approx(0.0)
        assert sfunc(n - 1.0) == tight_approx(eqReg.totalDistance())
