#!/usr/bin/env python

#make compatible 2.7 and 3
from __future__ import print_function
#get root
# from ROOT import TFile, TH1F, TDirectoryFile
import ROOT as root
#get the OS features
import os, sys, re, math
import numpy as np


def getClusterDistributions(file):

    disks,rings = 4,5
    h = [[root.TH1D() for j in range(rings)] for i in range(disks)]

    pileupstring = re.findall('summary_PU_(.*).root', file)
    pileup = float(pileupstring[0])
    print("Found a root file for pileup", pileup, "in file", file, "Objects: Clusters")

    rootfile = root.TFile.Open(file)
    rootfile.cd('BRIL_IT_Analysis/TEPX/Clusters')

    #build the histogram names
    histname = "Number of clusters for Disk "

    #outfile = root.TFile("Results_Clusters.root","RECREATE")

    #loop the disks
    for disk in range(1,5):
        histminusz = root.gDirectory.Get(histname +"-"+str(disk))
        histplusz = root.gDirectory.Get(histname+str(disk))
        #add plus and minus Z histograms
        histminusz.Add(histplusz)

        #now loop the rings
        for ring in range(5):
            h[disk-1][ring] = histminusz.ProjectionY("Ring", ring+1, ring+1)


    outfile = root.TFile("nClusters_"+str(pileup)+".root","RECREATE")
    for i in range(4):
        for j in range (5):
            #savecanvas = root.TCanvas("Clusters Disk "+str(disk)+" Ring "+str(ring+1),"Clusters Disk "+str(disk)+" Ring "+str(ring+1))
            #savecanvas.cd()
            h[i][j].Draw()
            outfile.WriteTObject(h[i][j],"Clusters Disk "+str(i+1)+" Ring "+str(j+1))
            #savecanvas.Write("Number of Clusters for Disk "+str(i+1)+" Ring "+str(j+1))

    outfile.Close()
    return


# def main():
args = len(sys.argv)
if args==3:
    path = sys.argv[1]
    observable = sys.argv[2]
else:
    print("Error, call with command line arguments: [1] path and [2] observable; the options for the latter are Clusters or 2x or 3x")
    path = "/afs/cern.ch/user/g/gauzinge/ITsim/mySummaryPlots/"
    observable ="Clusters"
    print("default values are:")
    print(path)
    print(observable)

print("Filepath", path, "Observable:",observable)
#now get the files in the path
files = os.listdir(path)
files = [item for item in files if not (item.find("summary") and (item.find(".root")))]
files.sort()
print(files)

if observable == "Clusters":

    for file in files:
        if file.find(".root"):
            filename = path+file
            #fill the actual graph for all available PU steps
            getClusterDistributions(filename)
        else:
            print("Not a root file, skipping")
            continue

