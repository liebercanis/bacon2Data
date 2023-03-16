#!/usr/bin/env python3
import sys
import os
import scp
from subprocess import Popen, PIPE
import pprint


def main(args):
    """ copy to nupacs1  """
    print(sys.argv)
    if (len(sys.argv) < 2):
        print("usage: copyDir.py  <date> ")
        return;
        
    myEnv = os.environ.copy()
    # init scp
    client = SCPClient(host="64.106.62.27", user="gold")
    client.use_system_keys()
    #print(myEnv)
    files = []
    tag =sys.argv[1] 
    theDir = 'rootData/'
    print("dir =  ",theDir, " file date  ", tag)
    p =os.listdir(theDir)
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
       # scp rootData/files[i]  gold@64.106.62.27:/data3/bacon/rootData/


if __name__ == '__main__':
    main(sys.argv[1:])
