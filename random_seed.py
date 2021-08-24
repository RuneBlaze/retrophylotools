import random
import os
import sys


def main():
    random.seed(None)
    sys.stdout.write("%d" % random.randint(1, 1000000000))
    sys.stdout.flush()
    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE


if __name__ == "__main__":
    main()
