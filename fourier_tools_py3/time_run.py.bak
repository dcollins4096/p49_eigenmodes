import time
from functools import wraps

def time_run(func):
    @wraps(func)
    def run(*args, **kwargs):
        t1 = time.time()
        retval = func(*args, **kwargs)
        t2 = time.time()
        print "Took %0.3e for %s" % (
                t2-t1, func.__name__)
        return retval
    return run

