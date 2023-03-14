#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE
import pprint


def main(args):
    """ run  anaDir ... ana on directory  """
    print(sys.argv)
    if (len(sys.argv) < 2):
        print("usage: anaDir <directory>")
        return;
        
    myEnv = os.environ.copy()
    #print(myEnv)
    files = []
    theDir = "rootData/"+sys.argv[1]
    print("dir = %s ",theDir)
    p = os.listdir('rootData')
    #print(p)
    for i in p:
        files.append(i)

    n = len(files)
    if (n < 1):
        return

    #print(" files %i ", len(p), " files %i ", len(files))
    if (len(sys.argv) > 2):
        n = int(sys.argv[2])

    print(" number of files to run  %i ", n)
    for i in range(0, n):
        print(" file %i ", i, " file %d", files[i])

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
