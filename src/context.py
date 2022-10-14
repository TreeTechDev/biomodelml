import sys

MAX_RECURSION = 2**(32-1)-1


class RecursionContext:
    def __init__(self, limit=MAX_RECURSION):
        self.limit = limit
        self.default_limit = sys.getrecursionlimit()

    def __enter__(self):
        sys.setrecursionlimit(self.limit)

    def __exit__(self, type, value, traceback):
        sys.setrecursionlimit(self.default_limit)
