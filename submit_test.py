import os
import sys


def get_matching_files(date_tag, rootdata_dir='rootData'):
    """Get list of files matching the date tag"""
    files = []
    try:
        for filename in os.listdir(rootdata_dir):
            if date_tag in filename:
                files.append(os.path.join(rootdata_dir, filename))
    except FileNotFoundError:
        print(f"Error: Directory '{rootdata_dir}' not found", file=sys.stderr)
        return []
    
    return sorted(files)


def main(args):

    print(sys.argv)
    if (len(sys.argv) < 2):
        print("usage: anaDir <date tag> <optional number of files> ")
        return;
        
    myEnv = os.environ.copy()
    #print(myEnv)
    tag =sys.argv[1] 
     
    theDir = 'rootData/'
    print("dir =  ",theDir, " file date  ", tag)

    files = get_matching_files(tag)
    n = len(files)
    ntot = n

     #print(" files %i ", len(p), " files %i ", len(files))
    if (len(sys.argv) > 2):
        n = int(sys.argv[2])



    print(" number of files to run  %i of %i  " % (n, ntot))
    for i in range(0, n):
        print(" file ", i, " file ", files[i]) 
    

    
if __name__ == '__main__':
    main(sys.argv[1:])

