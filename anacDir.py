#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE
import pprint


def main(args):
    """ run  anaDir ... ana on directory  """
    print(sys.argv)
    if (len(sys.argv) < 2):
        print("usage: anaDir <directory> <optional number of files> ")
        return;
        
    myEnv = os.environ.copy()
    #print(myEnv)
    files = []
    tag =sys.argv[1] 
    theDir = 'rootData/'
    print("dir =  ",theDir, " file date  ", tag)
    p = os.listdir(theDir)
    #print(p)
    for i in p:
        if i.find(tag) != -1:
            files.append(i)

    n = len(files)
    if (n < 1):
        return

    #print(" files %i ", len(p), " files %i ", len(files))
    if (len(sys.argv) > 2):
        n = int(sys.argv[2])

    print(" number of files to run  %i ", n)
    for i in range(0, n):
        print(" file ", i, " file ", files[i])


    for i in range(0, n):
        print(" run job %i ", i, " file %d", files[i])
        os.environ['LD_LIBRARY_PATH'] = os.getcwd()  # 
        process = Popen(['compiled/anac1', files[i]],
                        stdout=PIPE, stderr=PIPE, env=myEnv)
        stdout, stderr = process.communicate()
        process.wait()
        print(stdout)
        print(stderr)


if __name__ == '__main__':
    main(sys.argv[1:])
