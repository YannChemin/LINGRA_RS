import sys

class Foo(object):
    pass

def print_classes():
    current_module = sys.modules[__name__]
    for key in dir(current_module):
        if isinstance( getattr(current_module, key), type ):
            print(key)
