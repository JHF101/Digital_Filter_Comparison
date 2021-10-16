import numpy as np
import matplotlib.pyplot as plt
import glob
import os, time
from pylab import *
from matplotlib.pyplot import figure
from parameters import fs

continuousMode = True
multiPlot = True 

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 7 }

matplotlib.rc('font', **font)
fig = figure(figsize=(12, 12), dpi=60)
fig.tight_layout(pad=10)

def plotGraphs():
    # Get the names of all the textfiles in the current working directory
    names = glob.glob('*.txt')
    counter = 0
    for name in names:
        fileOutputArr = []
        f = open(name)
        for line in f:
            fields = line.split("\n")
            fileOutputArr.append(float(fields[0]))
        f.close()
        sampleSize = len(fileOutputArr)
        if multiPlot == False:
            # Detect frequency
            if str.lower(name[0])=="f":
                f = [(i/sampleSize)*fs for i in range(0,len(fileOutputArr))]
                plt.plot(f,fileOutputArr)
            
            # Detect time
            elif str.lower(name[0])=="t":
                t = np.arange(0,len(fileOutputArr)/fs,1/fs)
                plt.plot(t,fileOutputArr)

            else:
                print("Something went wrong")

            plt.title("{}".format(name[2:-4]))
            plt.grid()
            plt.show()
        
        else:
            # If you want multiple plots on one
            counter += 1
            # Detect frequency
            if str.lower(name[0])=="f":
                ax1 = subplot(len(names),1,counter)
                f = [(i/sampleSize)*fs for i in range(0,len(fileOutputArr))]
                ax1.plot(f,fileOutputArr)
                ax1.title.set_text("{}".format(name[2:-4]))
                # ax1.set_xlabel("Frequency(Hz)") -- Removed because it caused overlap in images

            # Detect time
            if str.lower(name[0])=="t":
                ax1 = subplot(len(names),1,counter)
                t = np.arange(0,len(fileOutputArr)/fs,1/fs)
                ax1.plot(t,fileOutputArr)
                ax1.title.set_text("{}".format(name[2:-4]))
                # ax1.set_xlabel("Time(s)") -- Removed because it caused overlap in images
            
            
    if multiPlot == True and len(names)>0:    
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

