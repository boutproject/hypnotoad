# Parallel map using multiprocessing, but using dill for pickling so that lambdas can be
# pickled. Inspired by example here https://stackoverflow.com/a/57190433

import dill
import multiprocessing


def _run_dill_encoded(payload):
    fun, args, kwargs = dill.loads(payload)
    res = fun(args, **kwargs)
    res = dill.dumps(res)
    return res


class ParallelMap:
    """
    Apply functions in parallel, using dill for pickling, inspired by example here
    https://stackoverflow.com/a/57190433

    Note: could not figure out a nice way of unpacking function arguments without
    breaking the case when a PsiContour is the only argument. In that case the
    PsiContour would be unpacked into many arguments, each of which is a Point2D. To
    work around this, all functions passed to ParallelMap.__call__() must take a single
    positional argument (which they can unpack inside the function).

    Parameters
    ----------
    np : int
        Number of processors to use
    """

    def __init__(self, np, *, equilibrium):
        if np > 1:
            self.task_queue = multiprocessing.Queue()
            self.result_queue = multiprocessing.Queue()
            equilibrium = dill.dumps(equilibrium)
            self.workers = [
                multiprocessing.Process(
                    target=ParallelMap.worker_run,
                    args=(self.task_queue, self.result_queue, equilibrium),
                )
                for i in range(np)
            ]
            for worker in self.workers:
                worker.start()
        else:
            self.workers = None
            self.equilibrium = equilibrium
            self.psi = equilibrium.psi
            self.f_R = equilibrium.f_R
            self.f_Z = equilibrium.f_Z

    def __del__(self):
        if self.workers is not None:
            for worker in self.workers:
                worker.terminate()
                worker.join()

    def worker_run(task_queue, result_queue, equilibrium):
        """
        'main' function for workers. Gets tasks from task_queue and puts results back in
        result_queue. Results need to be re-assembled into a list, so all tasks come
        with a list index, and results are returned with one. Re-assembly is done in
        __call__().
        """
        equilibrium = dill.loads(equilibrium)
        psi = equilibrium.psi
        f_R = equilibrium.f_R
        f_Z = equilibrium.f_Z
        while True:
            i, function, args, kwargs = task_queue.get()
            result = function(
                *args, equilibrium=equilibrium, psi=psi, f_R=f_R, f_Z=f_Z, **kwargs
            )
            result_queue.put((i, result))

    def __call__(self, function, args_list, **kwargs):
        if self.workers is None:
            return [
                function(
                    *args,
                    equilibrium=self.equilibrium,
                    psi=self.psi,
                    f_R=self.f_R,
                    f_Z=self.f_Z,
                    **kwargs
                )
                for args in args_list
            ]

        args_list = tuple(args_list)
        n_tasks = len(args_list)

        for i, args in enumerate(args_list):
            self.task_queue.put((i, function, args, kwargs))

        result = [None for i in range(n_tasks)]
        for count in range(n_tasks):
            i, this_result = self.result_queue.get()
            result[i] = this_result

        if not self.task_queue.empty():
            raise ValueError("Some tasks not finished")
        if not self.result_queue.empty():
            raise ValueError("Some results not handled")

        return result
