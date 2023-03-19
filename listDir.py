#!/usr/bin/env python3
import sys
import array as arr 
import os
import subprocess
from subprocess import Popen, PIPE
import pprint


def main(args):
    """ list runs in data  """
    dates = []
    theDir = 'data/'
    print("dir =  ",theDir, " date  ")
    p =os.listdir(theDir)

    for i in range(0,len(p)):
        s = p[i]
        c=s.count('_')
        if c>1 : 
            print(p[i])
    

if __name__ == '__main__':
    main(sys.argv[1:])
