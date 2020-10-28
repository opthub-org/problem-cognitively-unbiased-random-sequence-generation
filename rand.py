#!/bin/env python
import random
import sys


def main(n):
    for _ in range(n):
        print(random.randint(1, 6), end="")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: %s SEQLEN" % sys.argv[0])
        sys.exit(1)
    n = int(sys.argv[1])
    main(n)
