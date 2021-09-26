import numbers
import numpy


class MultiLocationArray(numpy.lib.mixins.NDArrayOperatorsMixin):
    """
    Container for arrays representing points at different cell locations
    Not all have to be filled.
    """

    _centre_array = None
    _xlow_array = None
    _ylow_array = None
    _corners_array = None

    def __init__(self, nx, ny):
        self.nx = nx
        self.ny = ny
        # Attributes that will be saved to output files along with the array
        self.attributes = {}

    @property
    def centre(self):
        if self._centre_array is None:
            self._centre_array = numpy.zeros([self.nx, self.ny])
        return self._centre_array

    @centre.setter
    def centre(self, value):
        if self._centre_array is None:
            self._centre_array = numpy.zeros([self.nx, self.ny])
        self._centre_array[...] = value

    @property
    def xlow(self):
        if self._xlow_array is None:
            self._xlow_array = numpy.zeros([self.nx + 1, self.ny])
        return self._xlow_array

    @xlow.setter
    def xlow(self, value):
        if self._xlow_array is None:
            self._xlow_array = numpy.zeros([self.nx + 1, self.ny])
        self._xlow_array[...] = value

    @property
    def ylow(self):
        if self._ylow_array is None:
            self._ylow_array = numpy.zeros([self.nx, self.ny + 1])
        return self._ylow_array

    @ylow.setter
    def ylow(self, value):
        if self._ylow_array is None:
            self._ylow_array = numpy.zeros([self.nx, self.ny + 1])
        self._ylow_array[...] = value

    @property
    def corners(self):
        if self._corners_array is None:
            self._corners_array = numpy.zeros([self.nx + 1, self.ny + 1])
        return self._corners_array

    @corners.setter
    def corners(self, value):
        if self._corners_array is None:
            self._corners_array = numpy.zeros([self.nx + 1, self.ny + 1])
        self._corners_array[...] = value

    def copy(self):
        new_multilocationarray = MultiLocationArray(self.nx, self.ny)
        if self.centre is not None:
            new_multilocationarray.centre = self.centre.copy()
        if self.xlow is not None:
            new_multilocationarray.xlow = self.xlow.copy()
        if self.ylow is not None:
            new_multilocationarray.ylow = self.ylow.copy()
        if self.corners is not None:
            new_multilocationarray.corners = self.corners.copy()

        return new_multilocationarray

    # The following __array_ufunc__ implementation allows the MultiLocationArray class to
    # be handled by Numpy functions, and add, subtract, etc. like an ndarray.
    # The implementation is mostly copied from the example in
    # https://docs.scipy.org/doc/numpy/reference/generated/numpy.lib.mixins.NDArrayOperatorsMixin.html#numpy.lib.mixins.NDArrayOperatorsMixin

    # One might also consider adding the built-in list type to this
    # list, to support operations like np.add(array_like, list)
    _HANDLED_TYPES = (numpy.ndarray, numbers.Number)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.get("out", ())
        for x in inputs + out:
            # Only support operations with instances of _HANDLED_TYPES.
            # Use MultiLocationArray instead of type(self) for isinstance to
            # allow subclasses that don't override __array_ufunc__ to
            # handle MultiLocationArray objects.
            if not isinstance(x, self._HANDLED_TYPES + (MultiLocationArray,)):
                return NotImplemented

        MLArrays = [self] + [x for x in inputs if isinstance(x, MultiLocationArray)]

        result = MultiLocationArray(self.nx, self.ny)

        # Defer to the implementation of the ufunc on unwrapped values.
        if all(x._centre_array is not None for x in MLArrays):
            this_inputs = tuple(
                x._centre_array if isinstance(x, MultiLocationArray) else x
                for x in inputs
            )
            if out:
                kwargs["out"] = tuple(
                    x.centre if isinstance(x, MultiLocationArray) else x for x in out
                )
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(
                        MultiLocationArray(self.nx, self.ny) for x in this_result
                    )
                for i, x in enumerate(this_result):
                    result[i].centre = x

            elif method == "at":
                # no return value
                result = None
            else:
                # one return value
                result.centre = this_result

        # Defer to the implementation of the ufunc on unwrapped values.
        if all(x._xlow_array is not None for x in MLArrays):
            this_inputs = tuple(
                x._xlow_array if isinstance(x, MultiLocationArray) else x
                for x in inputs
            )
            if out:
                kwargs["out"] = tuple(
                    x.xlow if isinstance(x, MultiLocationArray) else x for x in out
                )
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(
                        MultiLocationArray(self.nx, self.ny) for x in this_result
                    )
                for i, x in enumerate(this_result):
                    result[i].xlow = x

            elif method == "at":
                # no return value
                result = None
            else:
                # one return value
                result.xlow = this_result

        # Defer to the implementation of the ufunc on unwrapped values.
        if all(x._ylow_array is not None for x in MLArrays):
            this_inputs = tuple(
                x._ylow_array if isinstance(x, MultiLocationArray) else x
                for x in inputs
            )
            if out:
                kwargs["out"] = tuple(
                    x.ylow if isinstance(x, MultiLocationArray) else x for x in out
                )
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(
                        MultiLocationArray(self.nx, self.ny) for x in this_result
                    )
                for i, x in enumerate(this_result):
                    result[i].ylow = x

            elif method == "at":
                # no return value
                result = None
            else:
                # one return value
                result.ylow = this_result

        # Defer to the implementation of the ufunc on unwrapped values.
        if all(x._corners_array is not None for x in MLArrays):
            this_inputs = tuple(
                x._corners_array if isinstance(x, MultiLocationArray) else x
                for x in inputs
            )
            if out:
                kwargs["out"] = tuple(
                    x.corners if isinstance(x, MultiLocationArray) else x for x in out
                )
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(
                        MultiLocationArray(self.nx, self.ny) for x in this_result
                    )
                for i, x in enumerate(this_result):
                    result[i].corners = x

            elif method == "at":
                # no return value
                result = None
            else:
                # one return value
                result.corners = this_result

        return result

    def zero(self):
        # Initialise all locations, set them to zero and return the result
        self.centre = 0.0
        self.xlow = 0.0
        self.ylow = 0.0
        self.corners = 0.0
        return self
