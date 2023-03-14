#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE
import pprint


def main(args):
    """ run  anaDir ... ana on directory  """
    myEnv = os.environ.copy()
    print(myEnv)
    files = []
    p = os.listdir('rootData')
    print(p)
    for i in p:
        files.append(i)

    n = len(files)
    if (n < 1):
        return

    print(" files %i ", len(p), " files %i ", len(files))
    if (len(sys.argv) > 1):
        n = int(args[0])

    print(" args ", args, " number of files to run  ", n)
    for i in range(0, n):
        print(" run job %i ", i, " file %d", files[i])
        os.environ['LD_LIBRARY_PATH'] = os.getcwd()  # 
        process = Popen(['compiled/ana1', files[i]],
                        stdout=PIPE, stderr=PIPE, env=myEnv)
        stdout, stderr = process.communicate()
        process.wait()
        print(stdout)
        print(stderr)


if __name__ == '__main__':
    main(sys.argv[1:])
