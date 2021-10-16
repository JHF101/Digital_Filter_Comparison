import numpy as np
import math
from scipy import signal #To plot frequency response
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os, time
from pylab import *
from parameters import fs

continuousMode = True

logVal = True # Plot on logarithmic scale

def designFilter(zfiltNum,zfiltDen,stageName):
    freq, hz = signal.freqz(zfiltNum,zfiltDen, fs= fs, worN = 100000)
    fig, ax1 = plt.subplots()
    ax1.set_title(f'Digital filter frequency response {stageName}')
    if logVal == True:
        ax1.plot(freq, 20 * np.log10(abs(hz)), 'b')
    else:
        ax1.plot(freq, (abs(hz)), 'b')
    ax1.set_ylabel('Amplitude [dB]', color='b')
    ax1.set_xlabel('Frequency [Hz]')

    ax2 = ax1.twinx()
    angles = np.unwrap(np.angle(hz))
    ax2.plot(freq, angles, 'g')
    ax2.set_ylabel('Angle (radians)', color='g')
    ax2.grid()
    ax2.axis('tight')
    # plt.savefig(f"Filter Coefficencts of {stageName} stage")
    plt.show()


def plotGraphs():
    # Get the names of all the textfiles in the current working directory
    z_filter_num1 = []
    z_filter_den1 = []

    z_filter_num2 = []
    z_filter_den2 = []
    
    z_filter_num3 = []
    z_filter_den3 = []
    
    z_filter_num4 = []
    z_filter_den4 = []

    z_filter_num5 = []
    z_filter_den5 = []

    z_filter_num6 = []
    z_filter_den6 = []
    
    z_filter_num7 = []
    z_filter_den7 = []

    z_filter_num8 = []
    z_filter_den8 = []

    names = glob.glob('*.txt')
    counter = 0
    for name in names:
        fileOutputArr = []
        f = open(name)
        for line in f:
            fields = line.split("\n")
            fileOutputArr.append(float(fields[0]))
        f.close()

        print(str.lower(name[0:len(name)-5]))
        if str.lower(name[0:len(name)-5]) == "z_zfilternum":
            # Numerators
            if str.lower(name[len(name)-5]) == "1":
                print("At first Part")
                z_filter_num1 = fileOutputArr
                print(z_filter_num1)

            elif str.lower(name[len(name)-5]) == "2":
                z_filter_num2 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "3":
                z_filter_num3 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "4":
                z_filter_num4 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "5":
                z_filter_num5 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "6":
                z_filter_num6 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "7":
                z_filter_num7 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "8":
                z_filter_num8 = fileOutputArr

        if str.lower(name[0:len(name)-5]) == "z_zfilterden":
            # Denominators
            if str.lower(name[len(name)-5]) == "1":
                z_filter_den1 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "2":
                z_filter_den2 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "3":
                z_filter_den3 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "4":
                z_filter_den4 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "5":
                z_filter_den5 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "6":
                z_filter_den6 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "7":
                z_filter_den7 = fileOutputArr
            elif str.lower(name[len(name)-5]) == "8":
                z_filter_den8 = fileOutputArr

    if len(z_filter_den1)>0 and len(z_filter_num1)>0:
        designFilter(zfiltNum = z_filter_num1,
                        zfiltDen = z_filter_den1,
                        stageName = "First")
        print(z_filter_num1)
        print(z_filter_den1)
        
        
    if len(z_filter_den2)>0 and len(z_filter_num2)>0:
        designFilter(zfiltNum = z_filter_num2,
                        zfiltDen = z_filter_den2,
                        stageName = "Second")
        print(z_filter_num2)
        print(z_filter_den2)

    if len(z_filter_den3)>0 and len(z_filter_num3)>0:
        designFilter(zfiltNum = z_filter_num3,
                        zfiltDen = z_filter_den3,
                        stageName = "Third")
        print(z_filter_num3)
        print(z_filter_den3)

    if len(z_filter_den4)>0 and len(z_filter_num4)>0:
        designFilter(zfiltNum = z_filter_num4,
                     zfiltDen = z_filter_den4,
                     stageName = "Fourth")

    if len(z_filter_den5)>0 and len(z_filter_num5)>0:
        designFilter(zfiltNum = z_filter_num5,
                     zfiltDen = z_filter_den5,
                     stageName = "Fifth")

    if len(z_filter_den6)>0 and len(z_filter_num6)>0:
        designFilter(zfiltNum = z_filter_num6,
                    zfiltDen = z_filter_den6,
                    stageName = "Sixth")
    
    if len(z_filter_den7)>0 and len(z_filter_num7)>0:
        designFilter(zfiltNum = z_filter_num7,
                     zfiltDen = z_filter_den7,
                     stageName = "Seven")

    if len(z_filter_den8)>0 and len(z_filter_num8)>0:
        designFilter(zfiltNum = z_filter_num8,
                     zfiltDen = z_filter_den8,
                     stageName = "Eigth")
        
                
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

