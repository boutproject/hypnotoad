import pytest

from hypnotoad.utils.options import (
    OptionsFactory,
    WithMeta,
    NoneType,
    is_positive,
    is_positive_or_None,
    is_non_negative,
    is_non_negative_or_None,
)


class TestValueCheckUtilities:
    def test_is_positive(self):
        assert is_positive(1)
        assert not is_positive(-1)
        assert not is_positive(0.0)
        with pytest.raises(TypeError):
            is_positive(None)

    def test_is_positive_or_None(self):
        assert is_positive_or_None(1)
        assert not is_positive_or_None(-1)
        assert not is_positive_or_None(0.0)
        assert is_positive_or_None(None)

    def test_is_not_negative(self):
        assert is_non_negative(1)
        assert not is_non_negative(-1)
        assert is_non_negative(0.0)
        with pytest.raises(TypeError):
            is_non_negative(None)

    def test_is_not_negative_or_None(self):
        assert is_non_negative_or_None(1)
        assert not is_non_negative_or_None(-1)
        assert is_non_negative_or_None(0.0)
        assert is_non_negative_or_None(None)


class TestWithMeta:
    def test_init(self):
        x = WithMeta(1)
        assert x.value == 1

    def test_doc(self):
        x = WithMeta(3.0, doc="option x")
        assert x.doc == "option x"

    def test_value_type(self):
        x = WithMeta(3.0, value_type=float)
        assert x.value_type is float
        assert x.evaluate_expression({}) == 3.0

        x.value = 3
        with pytest.raises(TypeError):
            x.evaluate_expression({})

    def test_value_type_sequence(self):
        x = WithMeta(3.0, value_type=[float, NoneType])
        assert x.evaluate_expression({}) == 3.0

        x.value = None
        assert x.evaluate_expression({}) is None

        x.value = 3
        with pytest.raises(TypeError):
            x.evaluate_expression({})

    def test_allowed(self):
        x = WithMeta("foo", allowed="foo")
        assert x.value == "foo"
        assert x.allowed == ("foo",)
        assert x.evaluate_expression({}) == "foo"
        x.value = "baz"
        with pytest.raises(ValueError):
            x.evaluate_expression({})

    def test_allowed_sequence(self):
        x = WithMeta("foo", allowed=["foo", "bar"])
        assert x.value == "foo"
        assert x.allowed == ("foo", "bar")
        assert x.allowed != ["foo", "bar"]
        assert x.evaluate_expression({}) == "foo"
        x.value = "baz"
        with pytest.raises(ValueError):
            x.evaluate_expression({})

    def test_checks(self):
        x = WithMeta(4.0, checks=is_positive_or_None)
        assert x.value == 4.0
        assert x.evaluate_expression({}) == 4.0
        x.value = -2.0
        with pytest.raises(ValueError):
            x.evaluate_expression({})
        x.value = None
        assert x.evaluate_expression({}) is None

    def test_checks_sequence(self):
        x = WithMeta(5.0, checks=[is_positive, lambda x: x < 6.0])
        assert x.value == 5.0
        assert x.evaluate_expression({}) == 5.0
        x.value = -3.0
        with pytest.raises(ValueError):
            x.evaluate_expression({})
        x.value = 7.0
        with pytest.raises(ValueError):
            x.evaluate_expression({})

    def test_expression_allowed(self):
        x = WithMeta(lambda options: 2.0 * options["foo"], allowed=[4.0, 6.0])
        assert x.evaluate_expression({"foo": 2.0}) == 4.0
        assert x.evaluate_expression({"foo": 3.0}) == 6.0
        with pytest.raises(ValueError):
            x.evaluate_expression({"foo": 1.0})

    def test_expression_checks(self):
        x = WithMeta(
            lambda options: 2.0 + options["foo"],
            checks=[is_positive, lambda x: x < 10.0],
        )
        assert x.evaluate_expression({"foo": 2.0}) == 4.0
        assert x.evaluate_expression({"foo": 4.0}) == 6.0
        with pytest.raises(ValueError):
            x.evaluate_expression({"foo": -3.0})
        with pytest.raises(ValueError):
            x.evaluate_expression({"foo": 9.0})


class TestOptions:
    def test_defaults(self):
        factory = OptionsFactory(
            a=1,
            b=lambda options: options.a,
            c=lambda options: options["a"],
            d=lambda options: options.e + options.c,
            e=WithMeta(2.0, doc="option e", value_type=float, allowed=[2.0, 3.0]),
            f=WithMeta(
                11,
                doc="option f",
                value_type=int,
                checks=[is_positive, lambda x: x < 20],
            ),
            g=WithMeta(
                lambda options: options.a + 2,
                doc="option g",
                value_type=int,
                checks=[is_positive, lambda x: x < 20],
            ),
        )

        # test default values
        opts = factory.create()

        assert opts.a == 1
        assert opts.b == 1
        assert opts.c == 1
        assert opts.d == 3.0
        assert opts.e == 2.0
        assert opts.f == 11
        assert opts.g == 3

        assert opts["a"] == 1
        assert opts["b"] == 1
        assert opts["c"] == 1
        assert opts["d"] == 3.0
        assert opts["e"] == 2.0
        assert opts["f"] == 11
        assert opts["g"] == 3

        with pytest.raises(TypeError):
            opts.a = 2

        with pytest.raises(TypeError):
            opts["a"] = 2

        assert opts.is_default("a")
        assert opts.is_default("b")
        assert opts.is_default("c")
        assert opts.is_default("d")
        assert opts.is_default("e")
        assert opts.is_default("f")
        assert opts.is_default("g")
        with pytest.raises(KeyError):
            opts.is_default("x")

        assert "a" in opts
        assert "b" in opts
        assert "c" in opts
        assert "d" in opts
        assert "e" in opts
        assert "f" in opts
        assert "g" in opts
        assert not ("x" in opts)

        assert len(opts) == 7
        assert sorted([k for k in opts]) == sorted(["a", "b", "c", "d", "e", "f", "g"])
        assert sorted(opts.values()) == sorted([1, 1, 1, 3.0, 2.0, 11, 3])
        assert sorted(opts.items()) == sorted(
            [("a", 1), ("b", 1), ("c", 1), ("d", 3.0), ("e", 2.0), ("f", 11), ("g", 3)]
        )

    def test_initialise(self):
        factory = OptionsFactory(
            a=1,
            b=lambda options: options.a,
            c=lambda options: options["a"],
            d=lambda options: options.b + options.c,
            e=WithMeta(2.0, doc="option e", value_type=float, allowed=[2.0, 3.0]),
            f=WithMeta(
                11,
                doc="option f",
                value_type=int,
                checks=[is_positive, lambda x: x < 20],
            ),
            g=WithMeta(
                lambda options: options.a + 2,
                doc="option g",
                value_type=int,
                checks=[is_positive, lambda x: x < 20],
            ),
        )

        # test default values
        opts = factory.create({"a": 4, "b": 5, "e": 3.0, "f": 13, "z": 17})

        assert opts.a == 4
        assert opts.b == 5
        assert opts.c == 4
        assert opts.d == 9
        assert opts.e == 3.0
        assert opts.f == 13
        assert opts.g == 6

        # "z" should have been ignored
        with pytest.raises(AttributeError):
            opts.z

        assert opts["a"] == 4
        assert opts["b"] == 5
        assert opts["c"] == 4
        assert opts["d"] == 9
        assert opts["e"] == 3.0
        assert opts["f"] == 13
        assert opts["g"] == 6

        # "z" should have been ignored
        with pytest.raises(KeyError):
            opts["z"]

        with pytest.raises(TypeError):
            opts.a = 2

        with pytest.raises(TypeError):
            opts["a"] = 2

        assert not opts.is_default("a")
        assert not opts.is_default("b")
        assert opts.is_default("c")
        assert opts.is_default("d")
        assert not opts.is_default("e")
        assert not opts.is_default("f")
        assert opts.is_default("g")
        with pytest.raises(KeyError):
            opts.is_default("x")

        assert "a" in opts
        assert "b" in opts
        assert "c" in opts
        assert "d" in opts
        assert "e" in opts
        assert "f" in opts
        assert "g" in opts
        assert not ("x" in opts)

        assert len(opts) == 7
        assert sorted([k for k in opts]) == sorted(["a", "b", "c", "d", "e", "f", "g"])
        assert sorted(opts.values()) == sorted([4, 5, 4, 9, 3.0, 13, 6])
        assert sorted(opts.items()) == sorted(
            [("a", 4), ("b", 5), ("c", 4), ("d", 9), ("e", 3.0), ("f", 13), ("g", 6)]
        )

        with pytest.raises(ValueError):
            opts = factory.create({"e": 2.5})
        with pytest.raises(TypeError):
            opts = factory.create({"e": "2.0"})
        with pytest.raises(TypeError):
            opts = factory.create({"e": 2})
        with pytest.raises(ValueError):
            opts = factory.create({"f": -1})
        with pytest.raises(ValueError):
            opts = factory.create({"f": 30})
        with pytest.raises(TypeError):
            opts = factory.create({"f": 3.5})
        with pytest.raises(ValueError):
            opts = factory.create({"g": -7})
        with pytest.raises(ValueError):
            opts = factory.create({"g": 21})
        with pytest.raises(TypeError):
            opts = factory.create({"g": 3.5})
        with pytest.raises(ValueError):
            opts = factory.create({"a": -7})
        with pytest.raises(ValueError):
            opts = factory.create({"a": 21})
        with pytest.raises(TypeError):
            opts = factory.create({"a": 3.5})

    def test_circular(self):
        factory = OptionsFactory(
            a=lambda options: options.b, b=lambda options: options.a,
        )
        with pytest.raises(ValueError, match="Circular definition"):
            opts = factory.create()

        opts = factory.create({"b": 3})
        assert opts.a == 3
        assert opts.b == 3

        assert opts.is_default("a")
        assert not opts.is_default("b")
        with pytest.raises(KeyError):
            opts.is_default("x")

        assert "a" in opts
        assert "b" in opts
        assert not ("x" in opts)

        assert len(opts) == 2
        assert sorted([k for k in opts]) == sorted(["a", "b"])
        assert sorted(opts.values()) == sorted([3, 3])
        assert sorted(opts.items()) == sorted([("a", 3), ("b", 3)])
