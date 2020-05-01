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

from collections.abc import Sequence
from copy import deepcopy


# Some common tests that can be imported for clarity
def is_positive(x):
    return x > 0


def is_positive_or_None(x):
    if x is None:
        return True
    return is_positive(x)


def is_non_negative(x):
    return x >= 0


def is_non_negative_or_None(x):
    if x is None:
        return True
    return is_non_negative(x)


NoneType = type(None)


def with_default(value, default):

    if value is not None:
        return value

    return default


class WithMeta:
    """Type for passing metadata with options value or expression into *OptionsFactory

    """

    def __init__(self, value, *, doc=None, value_type=None, allowed=None, checks=None):
        """
        Parameters
        ----------
        value : expression, value, str, or WithMeta
            - If a callable expression is passed, evaluate expression(options) for the
              default value
            - If a value is passed, used as default value for this option
            - If (i) a str is passed and (ii) it is not in the allowed values for this
              option, and (iii) it is the name of another option, then set the default
              for this option as the value of the other option
            - If a WithMeta object is passed, check no other arguments were set and
              copy all attributes from value
        doc : str, optional
            Docstring for this option
        value_type : type, optional
            Type that this option should have
        allowed : value or sequence of values, optional
            When the option is set, it must have one of these values.
            Cannot be set if 'checks' is given.
        checks : expression or sequence of expressions, optional
            When a value is set for this option, all the expressions must return True
            when called with that value.
            Cannot be set if 'allowed' is given.
        """
        if isinstance(value, WithMeta):
            if (
                (doc is not None)
                or (value_type is not None)
                or (allowed is not None)
                or (checks is not None)
            ):
                raise ValueError(
                    f"doc={doc}, value_type={value_type}, allowed={allowed}, and "
                    f"checks={checks} should all be None when value is a WithMeta"
                )
            self.value = value.value
            self.doc = value.doc
            self.value_type = value.value_type
            self.allowed = value.allowed
            self.checks = value.checks
            return

        self.value = value
        self.doc = doc

        if isinstance(value_type, Sequence):
            value_type = tuple(value_type)
        self.value_type = value_type

        if (allowed is not None) and (checks is not None):
            raise ValueError("Cannot set both 'allowed' and 'checks'")

        if allowed is not None:
            if (not isinstance(allowed, Sequence)) or isinstance(allowed, str):
                # make allowed values a sequence
                allowed = (allowed,)
            self.allowed = tuple(allowed)
        else:
            self.allowed = None

        if checks is not None:
            if (not isinstance(checks, Sequence)) or isinstance(checks, str):
                # make checks expressions a sequence
                checks = (checks,)
            self.checks = tuple(checks)
            for check in self.checks:
                if not callable(check):
                    raise ValueError(
                        f"{check} is not callable, but was passed as a check"
                    )
        else:
            self.checks = None

    def __eq__(self, other):
        if not isinstance(other, WithMeta):
            return False
        return (
            self.value == other.value
            and self.doc == other.doc
            and self.allowed == other.allowed
            and self.checks == other.checks
        )

    def __str__(self):
        return (
            f"WithMeta({self.value}, doc={self.doc}, value_type={self.value_type}), "
            f"allowed={self.allowed}, checks={self.checks})"
        )

    def evaluate_expression(self, options, *, name=None):
        # Value may be expression or value. Try evaluating as an expression using options
        # first
        default_maybe_expression = self.value
        try:
            default_value = default_maybe_expression(options)
        except TypeError:
            # Try evaluating as name of another option
            if (
                isinstance(default_maybe_expression, str)
                and (
                    not isinstance(self.allowed, Sequence)
                    or default_maybe_expression not in self.allowed
                )
                and (default_maybe_expression in options)
            ):
                default_value = options[default_maybe_expression]
            else:
                default_value = default_maybe_expression

        return _checked(default_value, meta=self, name=name)


def _checked(value, *, meta=None, name=None):
    if (
        (meta is not None)
        and (meta.value_type is not None)
        and (not isinstance(value, meta.value_type))
    ):
        raise TypeError(
            f"{value} is not of type {meta.value_type}"
            f"{'' if name is None else ' for key=' + str(name)}"
        )

    if meta.allowed is not None:
        if value not in meta.allowed:
            raise ValueError(
                f"{value} is not in the allowed values {meta.allowed}"
                f"{'' if name is None else ' for key=' + str(name)}"
            )

    if meta.checks is not None:
        for check in meta.checks:
            if not check(value):
                raise ValueError(
                    f"The value {value}{'' if name is None else ' of key=' + str(name)} "
                    f"is not compatible with the checks"
                )

    return value


class OptionsFactory:
    """Factory to create Options instances

    """

    def __init__(self, *args, **kwargs):
        """Define the members of Options instances that this factory will create

        Parameters
        ----------
        *args : dicts of {key: [Withmeta, value or expression]}
            These dicts are combined with the kwargs to create the default values for
            this object.
            Intended to allow collecting defaults from contained objects. For example, if
            we have a class A, with members from classes B and C which each have an
            OptionsFactory, we could have something like:

                class A:
                    options_factory = OptionsFactory(
                        B.options_factory.defaults,
                        C.options_factory.defaults,
                        extra_option1 = 1,
                        extra_option2 = 2,
                    )

            It is an error for any keys in *args to be repeated or be the same as any in
            **kwargs.
        **kwargs : key=[WithMeta, value or expression]
            Keys are the names of the members of the Options that the factory will
            create.
            If a value is passed, that is used as the default for this key. If an
            expression is passed, it should take one argument, and can access values of
            other Options members from that argument. WithMeta allows values or
            expressions to be passed with extra metadata. For example,

                factory = OptionsFactory(
                    a=1,
                    b=lambda options: options.c + options.a
                    c=lambda options: 3*options["a"]
                    d=WithMeta(
                        4, doc="option d", value_type=int, allowed=[4, 5, 6]
                    )
                    e=WithMeta(
                        lambda options: options.a + options.b + options.c + options.d,
                        doc="option e",
                        value_type=float,
                        checks=lambda x: x > 0.0
                    )
                )
        """
        self.__default_values = {key: WithMeta(value) for key, value in kwargs.items()}

        # Add defaults from *args
        for a in args:
            for key, value in a.items():
                if key in self.__default_values:
                    if value != self.__default_values[key]:
                        raise ValueError(
                            f"{key} has been passed more than once with different values"
                        )
            self.__default_values.update(
                {key: WithMeta(value) for key, value in a.items()}
            )

        # Can be flipped (in create()) to True to allow values to be accessed using
        # __getitem__ or __getattr__
        self.__allow_access = False

    @property
    def defaults(self):
        """Get the default values defined for this OptionsFactory

        """
        return deepcopy(self.__default_values)

    def add(self, **kwargs):
        """Create a more specific version of the factory with extra options. For example,
        may be useful for a subclass like

            class Parent:
                options_factory = OptionsFactory(...)

            class Child:
                options_factory = Parent.options_factory.add(
                    an_extra_option="used only by Child"
                )

        Parameters
        ----------
        **kwargs : key=[WithMeta, value or expression]
            The new options to add, these override the ones in the parent factory if key
            already exists, but keep the doc, allowed and checks if the option is just a
            new value/expression (not a new WithMeta)
        """
        new_default_values = deepcopy(self.__default_values)
        for key, value in kwargs.items():
            if (not isinstance(value, WithMeta)) and key in new_default_values:
                # just update the default value or expression
                new_default_values[key].value = value
            else:
                new_default_values[key] = value

        return OptionsFactory(new_default_values)

    def __getitem__(self, key):
        if not self.__allow_access:
            raise AttributeError("Cannot access values from an OptionsFactory")

        # If key is already in __values, then it has a definite value
        try:
            return self.__values[key]
        except KeyError:
            pass

        # When setting default values, detect circular definitions
        try:
            self.__key_chain
        except AttributeError:
            chain_start = True
            self.__key_chain = [key]
        else:
            if key in self.__key_chain:
                # Found a circular definition

                # Tidy up factory state
                key_chain = self.__key_chain
                del self.__key_chain
                self.__allow_access = False

                # Tell the user where the circular definition was
                index = key_chain.index(key)
                raise ValueError(
                    f"Circular definition of default values. At least one of "
                    f"{key_chain[index:]} must have a definite value"
                )

            chain_start = False
            self.__key_chain.append(key)

        self.__values[key] = self.__default_values[key].evaluate_expression(
            self, name=key
        )

        if chain_start:
            # Tidy up temporary member variable
            del self.__key_chain

        return self.__values[key]

    def __getattr__(self, key):
        if self.__allow_access:
            if key in self.__values or key in self.__default_values:
                return self.__getitem__(key)
        return super().__getattr__(key)

    def __contains__(self, key):
        return key in self.__default_values

    def create(self, values={}):
        """Create an Options instance

        The members of the created Options are defined by this
        OptionsFactory instance. Any values passed in the values dict are used,
        and the rest are set from defaults, which can be expressions depending on other
        members.

        Parameters
        ----------
        values : dict, optional
            Non-default values to be used
        """

        # do not modify passed-in values
        values = deepcopy(dict(values))

        # ignore values for keys not in the list of keys defined in the factory
        for key in list(values.keys()):
            if key not in self.__default_values:
                del values[key]

        # check passed-in values
        self.__values = {
            key: _checked(value, meta=self.__default_values[key], name=key,)
            for key, value in values.items()
        }

        # make a list of the explicitly-set (non-default) values
        non_default = [k for k in self.__values]

        # Return new Options instance
        # Now create all other values from default expressions. Default value evalutation
        # is implemented in the __getitem__ method so we can pass 'self' to the default
        # expressions.
        # Set __allow_access to True to allow internal use of __getitem__, which is
        # forbidden externally
        self.__allow_access = True
        for key in self.__default_values:
            # Note self.__values[key] is updated in self.__getitem__(key)
            self[key]

        # Create Options to return
        new_instance = OptionsFactory.Options(
            self.__values,
            {key: self.__default_values[key].doc for key in self.__default_values},
            non_default,
        )

        self.__allow_access = False

        # Clean up state of the factory
        del self.__values

        return new_instance

    class Options:
        """Provide access to a pre-defined set of options, with values fixed when the
        instance is created

        """

        __frozen = False

        def __init__(self, data, doc, non_default):
            self.__data = data
            self.__doc = doc
            self.__non_default = non_default

            # Set self.__frozen to True to prevent attributes being changed
            self.__frozen = True

        @property
        def doc(self):
            return deepcopy(self.__doc)

        def __getitem__(self, key):
            return deepcopy(self.__data.__getitem__(key))

        def __setitem__(self, key, value):
            raise TypeError("Options does not allow assigning to keys")

        def __getattr__(self, key):
            if key == "_Options__data":
                # need to treat __data specially, as we use it for the next test
                return super.__getattr__(key)
            if key in self.__data:
                return self.__getitem__(key)
            return super.__getattr__(key)

        def __setattr__(self, key, value):
            if self.__frozen:
                raise TypeError("Options does not allow assigning to attributes")
            super().__setattr__(key, value)

        def is_default(self, key):
            if key not in self.__data:
                raise KeyError(f"{key} is not in this Options")
            return key not in self.__non_default

        def __contains__(self, key):
            return key in self.__data

        def __len__(self):
            return len(self.__data)

        def __iter__(self):
            return iter(deepcopy(self.__data))

        def keys(self):
            return [key for key in self]

        def values(self):
            return [deepcopy(v) for v in self.__data.values()]

        def items(self):
            return zip(self.keys(), self.values())

        def __str__(self):
            return str(self.__data)


# Helper function to convert options to string
def optionsTableString(options):
    """Return a string containing a table of options set"""
    formatstring = "{:<50}|  {:<30}\n"

    # Header
    result = (
        "\nOptions\n=======\n" + formatstring.format("Name", "Value") + "-" * 80 + "\n"
    )

    # Row for each value
    for name, value in sorted(options.items()):
        valuestring = str(value)
        if options.is_default(name):
            valuestring += "\t(default)"
        result += formatstring.format(name, valuestring)
    return result
