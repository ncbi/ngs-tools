#!/usr/bin/python
import time

def timed(fn):
    
    def decorated(*args, **kwargs):
        before = time.time()
        try:
#			print "start", fn.__name__
			return fn(*args, **kwargs)
        finally:
            after = time.time()
            print '%s execution took %.2f s' % (fn.__name__, after-before)
    
    return decorated

