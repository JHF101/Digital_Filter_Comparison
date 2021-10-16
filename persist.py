import numpy as np
import matplotlib.pyplot as plt
import glob
import os, time
from pylab import *
from parameters import fs

continuousMode = True

def plotGraphs():
    # Get the names of all the textfiles in the current working directory
    names = glob.glob('*.txt')
    counter = 0
    if len(names)>0:
        for name in names:
            fileOutputArr = []
            f = open(name)
            for line in f:
                fields = line.split("\n")
                fileOutputArr.append(float(fields[0]))
            f.close()
            sampleSize = len(fileOutputArr)
            
            # Detect frequency
            if str.lower(name[0])=="f":
                f = [(i/sampleSize)*fs for i in range(0,len(fileOutputArr))]
                plt.plot(f,fileOutputArr)
            else:
                print("Something went wrong")
        plt.xlabel("Frequency(Hz)")
        plt.ylabel("Magnitude")
        plt.title("Persist Function")
        plt.grid()
        plt.show()

if continuousMode == True:
    path_to_watch = str(os.getcwd()) # Watch the current directory
    before = dict ([(f, None) for f in os.listdir (path_to_watch)])
    while 1:
        time.sleep(5) # Give time to check
        after = dict ([(f, None) for f in os.listdir (path_to_watch)]) 
        added = [f for f in after if not f in before] # Checks for files added
        removed = [f for f in before if not f in after] # Checks for files removed
        if added: 
            print("Added: ", ", ".join (added))
        if removed: 
            print ("Removed: ", ", ".join (removed))

        plotGraphs() # Run the function

        before = after

        # Just deletes old textfiles
        nametodel = glob.glob('*.txt')
        if nametodel != None:
            for name in nametodel:
                os.remove(name)
else: 
    plotGraphs()

