import time


def time_it(method):
    # TODO: docstring
    def timed(*args, **kw):
        start_time = time.time()
        result = method(*args, **kw)
        end_time = time.time()
        time_elapsed = end_time - start_time
        return result, time_elapsed

    return timed
