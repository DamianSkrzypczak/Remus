import time


def time_it(method):
    def timed(*args, **kw):
        start_time = time.time()
        result = method(*args, **kw)
        end_time = (time.time() - start_time)
        return result, round(end_time, 3)

    return timed
