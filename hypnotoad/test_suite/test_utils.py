from hypnotoad.utils.utils import with_default


class TestUtils:
    def test_with_default(self):
        assert with_default(1.0, 2.0) == 1.0
        assert type(with_default(1.0, 2.0)) is float
        assert with_default(None, 2.0) == 2.0
        assert type(with_default(None, 2.0)) is float
        assert with_default(-1.0, 2.0) == -1.0
        assert type(with_default(-1.0, 2.0)) is float
        assert with_default(1.0, -2.0) == 1.0
        assert type(with_default(1.0, -2.0)) is float
        assert with_default(None, -2.0) == -2.0
        assert type(with_default(None, -2.0)) is float
        assert with_default(0.0, 2.0) == 0.0
        assert type(with_default(0.0, 2.0)) is float
        assert with_default(None, 0.0) == 0.0
        assert type(with_default(None, 0.0)) is float
        assert with_default(0, 2.0) == 0
        assert type(with_default(0, 2.0)) is int
        assert with_default(None, 0) == 0
        assert type(with_default(None, 0)) is int
        assert with_default(False, 2.0) is False
        assert with_default(True, 2.0) is True
