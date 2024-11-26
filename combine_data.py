import shutil
import time, glob

outfilename = "output.txt"

with open(outfilename, 'wb') as outfile:
    for filename in glob.glob('/home/weinfub/CSCI597/project/DIMHPA/DATA/*'):

        with open(filename, 'rb') as readfile:
            shutil.copyfileobj(readfile, outfile)