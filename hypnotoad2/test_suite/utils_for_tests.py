import pytest


def tight_approx(val):
    return pytest.approx(val, rel=1.0e-12, abs=1.0e-13)
