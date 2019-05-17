import pytest

def tight_approx(val):
    return pytest.approx(val, rel=1.e-12, abs=1.e-13)
