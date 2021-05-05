# Parallel map using multiprocessing, but using dill for pickling so that lambdas can be
# pickled. Inspired by example here https://stackoverflow.com/a/57190433

import dill
import itertools
import multiprocessing


def _run_dill_encoded(payload):
    fun, args = dill.loads(payload)
    res = fun(args)
    res = dill.dumps(res)
    return res


class ParallelMap:
    """
    Apply functions in parallel, using dill for pickling, inspired by example here
    https://stackoverflow.com/a/57190433

    Note: could not figure out a nice way of unpacking function arguments
    without breaking the case when a PsiContour is the only argument. In that
    case the PsiContour would be unpacked into many arguments, each of which is
    a Point2D. To work around this, all functions passed to
    ParallelMap.__call__() must take a single positional argument (which they
    can unpack inside the function).

    Parameters
    ----------
    np : int
        Number of processors to use
    """

    def __init__(self, np):
        if np > 1:
            self.pool = multiprocessing.Pool(np)
        else:
            self.pool = None

    def __del__(self):
        # Not sure why this is necessary, but without it get an `OSError:
        # [Errno 9] Bad file descriptor` thrown by multiprocessing
        del self.pool

    def __call__(self, function, args, **kwargs):
        if self.pool is None:
            return [function(arg, **kwargs) for arg in args]

        it = map(
            dill.dumps,
            zip(itertools.cycle([lambda x: function(x, **kwargs)]), args),
        )

        map_result = self.pool.map(_run_dill_encoded, it)

        return list(map(dill.loads, map_result))
