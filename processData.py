import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
from bisect import bisect

#Load the data
data = np.loadtxt("database/SimonMurphyCatalog.csv", dtype=np.str_)

dataKIC = []
for i in range(len(data)):
    if "a" in data[i][0]:
        dataKIC.append(int(data[i][0][:-1]))
    elif "b" in data[i][0]:
        dataKIC.append(int(data[i][0][:-1]))
    else:
        dataKIC.append(int(data[i][0]))

dataKIC = np.array(dataKIC)
AllKIC = data[:,0]
AllOrbitalPeriod = data[:,1]
AllOrbitalPeriodErr = data[:,2]
AllFileNames = glob.glob("data/*.txt")


for FileItem in AllFileNames:
    data= pd.read_csv(FileItem)
    
    #Lomb Scargle Periodogram
    frequency, power = LombScargle(data['Time'], data[' PDC Flux']).autopower()
    LS_Period = 1/frequency
    ArrangedIndex = np.argsort(LS_Period)
    LS_Period = LS_Period[ArrangedIndex]
    power = power[ArrangedIndex]

    SelectRange = (LS_Period > 0.05) & (LS_Period < 100)
    LS_Period = LS_Period[SelectRange]
    power = power[SelectRange]
    

    MedianPower = np.median(power)
    AllPeaks = find_peaks(power, distance=15,prominence=10.*MedianPower)
    PeakLocations = AllPeaks[0] 
    PeakProminence = AllPeaks[1]['prominences']
    ArrangeIndex = np.argsort(PeakProminence)[::-1]
    PeakLocations = PeakLocations[ArrangeIndex]
    PeakProminence = PeakProminence[ArrangeIndex]

    NumPoints = 100

    #Work on this tomorrow
    PeriodPeaks = LS_Period[PeakLocations][:100]
    PowerPeaks = power[PeakLocations][:100]
    ArrangeIndex = np.argsort(PeriodPeaks)

    PeriodPeaks = PeriodPeaks[ArrangeIndex]
    PowerPeaks = PowerPeaks[ArrangeIndex]

  

    KIC_ID = int(FileItem.split("/")[1][:9])
    SelectIndex = KIC_ID == dataKIC
    
    if np.sum(SelectIndex) == 0:
        print("No data for KIC_ID: ", KIC_ID)
        continue

    CurrentPeriod = float(AllOrbitalPeriod[SelectIndex][0])
    CurrentPeriodErr = float(AllOrbitalPeriodErr[SelectIndex][0])
    print("KIC_ID: ", KIC_ID, "Orbital Period: ", CurrentPeriod, "Orbital Period Error: ", CurrentPeriodErr)    
    input("Wait here.e.re.") 
    
    PeriodRange = np.linspace(CurrentPeriod-4*CurrentPeriodErr, CurrentPeriod+4*CurrentPeriodErr, 10000)
    PowerValues = np.zeros(len(PeriodRange))
    
    for p in PeriodRange:
        #Check harmonics in the peal
        #HighestHarmonic = int(0.01/
        Harmonics = np.arange(1, 10000)

        print("The value of h in harmonics is: ", Harmonics)
        for h in Harmonics:
            Index = bisect(PeriodPeaks, h)

            print("The value of h is: ", h)
            #print(h-LS_Period[Index-1], LS_Period[Index]-h)

            #PowerValues.append(power[Index]) 
            plt.figure(figsize=(10, 5))
            plt.plot(LS_Period, power, "k-")    
            for i in range(len(PeriodPeaks)):
                plt.plot(PeriodPeaks[i], PowerPeaks[i], "r+")
            plt.axvline(x=h, color="b", lw=3)
            plt.axvline(x=PeriodPeaks[Index], color="g", lw=3)
            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("Period (Days)")
            plt.ylabel("Power")
            plt.title("Lomb Scargle Periodogram")
            plt.show()


        





    input("Wait here...")