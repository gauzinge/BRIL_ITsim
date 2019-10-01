#!/usr/bin/env python

#make compatible 2.7 and 3
from __future__ import print_function
#get root
# from ROOT import TFile, TH1F, TDirectoryFile
import ROOT as root
#get the OS features
import os, sys, re, math

#extract mean and sigma from 1D projections of # of Clusters histograms
def getParams(hist, ring):
    ringhist= hist.ProjectionY("Ring", ring+1, ring+1)
    # (fit,status)=fitPoisson(ringhist)
    # if status != 4000:
    mean = ringhist.GetMean()
        # print("Problem with the fit, using simple mean - fit status:",status)
    # else:
        # mean = fit.GetParameter(1)
    if mean == 0.5:
        mean = 0
        sigma = 0
    else:
        sigma = math.sqrt(mean)
        # sigma = ringhist.GetRMS()
    # if mean == 0.5:
        # mean = 0

    print("Ring",ring," wiht mean",mean,"RMS",sigma)
    return (mean,sigma)

#extract mean and sigma from 1D projection of global histograms
def getGlobalParams(hist):
    ringhist= hist.ProjectionY("Global Projection", 0, -1)
    mean = ringhist.GetMean()
    if mean == 0.5:
        mean = 0
        sigma = 0
    else:
        sigma = math.sqrt(mean)
    return (mean,sigma)

#get the linearity graph for clusters
def getLinearityClusters(file, graphs=[], globalgraph=[]):

    pileupstring = re.findall('summary_PU_(.*).root', file)
    pileup = float(pileupstring[0])
    print("Found a root file for pileup", pileup, "in file", file, "Objects: Clusters")

    rootfile = root.TFile.Open(file)
    rootfile.cd('BRIL_IT_Analysis/TEPX/Clusters')

    #build the histogram names
    histname = "Number of clusters for Disk "

    globalhist = root.gDirectory.Get(histname+"-1")
    globalhist.Add(root.gDirectory.Get(histname+"1"))
    for disk in range(2,5):
        hmz = root.gDirectory.Get(histname +"-"+str(disk))
        hpz = root.gDirectory.Get(histname+str(disk))
        hmz.Add(hpz)
        globalhist.Add(hmz)

    #trigger = 0.075
    #(globalmean, globalsigma) = getGlobalParams(globalhist)
    #if pileup==0.0:
    #    globalgraph.SetPoint(globalgraph.GetN(),pileup, 0.0)
    #    globalgraph.SetPointError(globalgraph.GetN()-1,0, 0.0)
    #else:
    #    globalgraph.SetPoint(globalgraph.GetN(), pileup, 
    #                         (math.sqrt(globalmean*pow(2,14))/(globalmean*pow(2,14)))*(math.sqrt(trigger/40)/(trigger/40)) )
    #    globalgraph.SetPointError(globalgraph.GetN()-1,0,0) #1/(2*globalmean*globalsigma))

    #loop the disks
    sumofmeans = 0
    sumofstats = 0
    trigger = 0.075
    for disk in range(1,5):
        histminusz = root.gDirectory.Get(histname +"-"+str(disk))
        histplusz = root.gDirectory.Get(histname+str(disk))
        #add plus and minus Z histograms
        histminusz.Add(histplusz)

        #now loop the rings
        for ring in range(5):
            if ring==0:
                if disk==4:
                    trigger = 0.825
            else:
                trigger = 0.075
            (mean, sigma) = getParams(histminusz, ring)
            sumofmeans = sumofmeans + mean
            if pileup==0.0:
                graphs[disk-1][ring].SetPoint(graphs[disk-1][ring].GetN(),pileup, 0.0)
                graphs[disk-1][ring].SetPointError(graphs[disk-1][ring].GetN()-1,0, 0.0)
                sumofstats = sumofstats + 0
            else:    
                graphs[disk-1][ring].SetPoint(graphs[disk-1][ring].GetN(),pileup,
                                              (math.sqrt(mean*pow(2,14))/(mean*pow(2,14)))*(math.sqrt(trigger/40)/(trigger/40)) )
                graphs[disk-1][ring].SetPointError(graphs[disk-1][ring].GetN()-1,0,0) # 1/(2*mean*sigma))
                sumofstats = sumofstats + (math.sqrt(mean*pow(2,14))/(mean*pow(2,14)))*(math.sqrt(trigger/40)/(trigger/40))

    meanofmeans = sumofmeans/20
    meanofstats = sumofstats/20
    if pileup==0.0:
        globalgraph.SetPoint(globalgraph.GetN(),pileup, 0.0)
        globalgraph.SetPointError(globalgraph.GetN()-1,0, 0.0)
    else:
        globalgraph.SetPoint(globalgraph.GetN(), pileup,
                             (math.sqrt(meanofmeans*pow(2,14))/(meanofmeans*pow(2,14)))*(math.sqrt(trigger/40)/(trigger/40)) )
        globalgraph.SetPointError(globalgraph.GetN()-1,0,0)

    rootfile.Close()
    return

#get the linearity graph for hits
def getLinearityHits(file, graphs=[], globalgraph=[]):

    pileupstring = re.findall('summary_PU_(.*).root', file)
    pileup = float(pileupstring[0])
    print("Found a root file for pileup", pileup, "in file", file, "Objects: Hits")

    rootfile = root.TFile.Open(file)
    rootfile.cd('BRIL_IT_Analysis/TEPX/Hits')

    #build the histogram names
    histname = "Number of hits for Disk "

    globalhist = root.gDirectory.Get(histname+"-1")
    globalhist.Add(root.gDirectory.Get(histname+"1"))
    for disk in range(2,5):
        hmz = root.gDirectory.Get(histname +"-"+str(disk))
        hpz = root.gDirectory.Get(histname+str(disk))
        hmz.Add(hpz)
        globalhist.Add(hmz)

    (globalmean, globalsigma) = getGlobalParams(globalhist)
    if pileup==0.0:
        globalgraph.SetPoint(globalgraph.GetN(),pileup, 0.0)
        globalgraph.SetPointError(globalgraph.GetN()-1,0, 0.0)
    else:
        globalgraph.SetPoint(globalgraph.GetN(), pileup, math.sqrt(globalmean)/globalmean)
        globalgraph.SetPointError(globalgraph.GetN()-1, 0, 1/(2*globalmean*globalsigma))

    #loop the disks
    for disk in range(1,5):
        histminusz = root.gDirectory.Get(histname +"-"+str(disk))
        histplusz = root.gDirectory.Get(histname+str(disk))
        #add plus and minus Z histograms
        histminusz.Add(histplusz)

        #now loop the rings
        for ring in range(5):
            (mean, sigma) = getParams(histminusz, ring)
            if pileup==0.0:
                graphs[disk-1][ring].SetPoint(graphs[disk-1][ring].GetN(),pileup, 0.0)
                graphs[disk-1][ring].SetPointError(graphs[disk-1][ring].GetN()-1,0, 0.0)
            else:
                graphs[disk-1][ring].SetPoint(graphs[disk-1][ring].GetN(),pileup, math.sqrt(mean)/mean)
                graphs[disk-1][ring].SetPointError(graphs[disk-1][ring].GetN()-1,0, 1/(2*mean*sigma))

    rootfile.Close()
    return


# def main():
disks,rings = 4,5
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

# a TCanvas
c_canvas = root.TCanvas("Summary","Summary")
c_canvas.Divide(5,4)

# a 2D array of TGraph Errors
#first index is disk and second is ring
#first, let's only deal with the case for Clusters
if observable == "Clusters":
    graphs = [[root.TGraphErrors() for j in range(rings)] for i in range(disks)]
    extrapolated = [[root.TF1() for j in range(rings)] for i in range(disks)]
    errors = [[root.TGraphErrors() for j in range(rings)] for i in range(disks)]

    globalgraph = root.TGraphErrors()

    for file in files:
        if file.find(".root"):
            filename = path+file
            #fill the actual graph for all available PU steps
            getLinearityClusters(filename,graphs,globalgraph)
        else:
            print("Not a root file, skipping")
            continue

    rootfile = root.TFile("Results_Clusters_StatError.root","RECREATE")
    globalgraph.SetLineColor(1)
    globalgraph.SetMarkerStyle(20)
    globalgraph.SetMarkerColor(1)
    globalgraph.SetMarkerSize(1)
    globalgraph.SetTitle("Stat error size - All disks;Pileup;Statistical error")
    saveglobalcanvas = root.TCanvas("Stat error size - All disks","Stat error size - All disks")
    saveglobalcanvas.cd()
    globalgraph.Draw("ap")
    saveglobalcanvas.Write("StatErrorSize_AllDisks")
    index = 1
    for i in range(disks):
        for j in range(rings):
            #Cosmetics
            graphs[i][j].SetLineColor(1)
            graphs[i][j].SetMarkerStyle(20)
            graphs[i][j].SetMarkerColor(1)
            graphs[i][j].SetMarkerSize(1)
            graphs[i][j].SetTitle("Stat error size - Disk "+str(i+1)+" Ring "+str(j+1)+";Pileup;Statistical error")
            c_canvas.cd(index)
            graphs[i][j].Draw("ap")

            #fit and extrapolate
            #(extrapolated[i][j],errors[i][j]) = extrapolateLinear(graphs[i][j],2)
            #errors[i][j].Draw("e3 same")
            #extrapolated[i][j].Draw("same")

            #calculate relative nonlinearity
            #deviation = relativeNonlinearity(graphs[i][j], extrapolated[i][j])
            #deviation.Write("Deviation Clusters Disk"+str(i+1)+"Ring"+str(j+1))

            #save canvases for the individual disk/ring combos
            savecanvas = root.TCanvas("Stat error size - Disk "+str(i+1)+" Ring "+str(j+1),"Stat error size - Disk "+str(i+1)+" Ring "+str(j+1))
            savecanvas.cd()
            graphs[i][j].Draw("ap")
            #errors[i][j].Draw("e3 same")
            #extrapolated[i][j].Draw("same")
            savecanvas.Write("Stat error size - Disk "+str(i+1)+" Ring "+str(j+1))
            index = index+1

    #Write out the summary as well
    c_canvas.Write("StatErrorSize_SummaryClusters")
    rootfile.Close()

#second, let's only deal with the case for Hits
elif observable == "Hits":
    graphs = [[root.TGraphErrors() for j in range(rings)] for i in range(disks)]
    extrapolated = [[root.TF1() for j in range(rings)] for i in range(disks)]
    errors = [[root.TGraphErrors() for j in range(rings)] for i in range(disks)]

    globalgraph = root.TGraphErrors()

    for file in files:
        if file.find(".root"):
            filename = path+file
            #fill the actual graph for all available PU steps
            getLinearityHits(filename,graphs,globalgraph)
        else:
            print("Not a root file, skipping")
            continue

    rootfile = root.TFile("Results_Hits_StatError.root","RECREATE")
    globalgraph.SetLineColor(1)
    globalgraph.SetMarkerStyle(20)
    globalgraph.SetMarkerColor(1)
    globalgraph.SetMarkerSize(1)
    globalgraph.SetTitle("Stat error size - All disks;Pileup;#sqrt{n_hits}/n_hits")
    saveglobalcanvas = root.TCanvas("Stat error size - All disks","Stat error size - All disks")
    saveglobalcanvas.cd()
    globalgraph.Draw("ap")
    saveglobalcanvas.Write("StatErrorSize_AllDisks")
    index = 1
    for i in range(disks):
        for j in range(rings):
            #Cosmetics
            graphs[i][j].SetLineColor(1)
            graphs[i][j].SetMarkerStyle(20)
            graphs[i][j].SetMarkerColor(1)
            graphs[i][j].SetMarkerSize(1)
            graphs[i][j].SetTitle("Stat error size - Disk "+str(i+1)+" Ring "+str(j+1)+";Pileup;#sqrt{n_hits}/n_hits")
            c_canvas.cd(index)
            graphs[i][j].Draw("ap")

            #fit and extrapolate
            #(extrapolated[i][j],errors[i][j]) = extrapolateLinear(graphs[i][j],2)
            #errors[i][j].Draw("e3 same")
            #extrapolated[i][j].Draw("same")

            #calculate relative nonlinearity
            #deviation = relativeNonlinearity(graphs[i][j], extrapolated[i][j])
            #deviation.Write("Deviation Hits Disk"+str(i+1)+"Ring"+str(j+1))

            #save canvases for the individual disk/ring combos
            savecanvas = root.TCanvas("Stat error size - Disk "+str(i+1)+" Ring "+str(j+1),"Stat error size - Disk "+str(i+1)+" Ring "+str(j+1))
            savecanvas.cd()
            graphs[i][j].Draw("ap")
            #errors[i][j].Draw("e3 same")
            #extrapolated[i][j].Draw("same")
            savecanvas.Write("Stat error size - Disk "+str(i+1)+" Ring "+str(j+1))
            index = index+1

    #Write out the summary as well
    c_canvas.Write("StatErrorSize_SummaryHits")
    rootfile.Close()

# if __name__ == '__main__':
        # main()
