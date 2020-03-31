#!/usr/bin/env python

#make compatible 2.7 and 3
from __future__ import print_function
#get root
# from ROOT import TFile, TH1F, TDirectoryFile
import ROOT as root
#get the OS features
import os, sys, re, math
import numpy as np

#fit a poisson distribution
def fitPoisson(hist):
    fit = root.TF1("fit","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)", 1, 2500)
    fit.SetParameters(1,1,1)
    if hist.GetMean() > 0.51:
        status = int(hist.Fit("fit","RMQN"))
    else:
        status = 0
        # if status != 4000:
            # print(status)
            # hist.Draw()
            # raw_input("Press Enter to Continue...")
    return fit, status

#extract mean/ and sigma from 1D projections of # of Clusters histograms
def getParams(hist, ring):
    ringhist= hist.ProjectionY("Ring", ring+1, ring+1)
    # (fit,status)=fitPoisson(ringhist)
    # if status != 4000:
    prob = np.array([0.5])
    x = 0.0
    q = np.array([0.])
    x = ringhist.GetQuantiles(1,q,prob)
    median = q
    mean = ringhist.GetMean()
    print("mean ",mean," median ", median)
        # print("Problem with the fit, using simple mean - fit status:",status)
    # else:
        # mean = fit.GetParameter(1)
    #if mean == 0.5:
    #    mean = 0
    #    sigma = 0
    #else:
        #sigma = math.sqrt(mean)
    sigma = ringhist.GetRMS()
    # if mean == 0.5:
        # mean = 0

    print("Ring",ring," wiht m",mean,"RMS",sigma)
    return (mean,sigma)

#get the linearity graph for clusters
def getLinearityClusters(file, graphs=[]):

    pileupstring = re.findall('summary_PU_(.*).root', file)
    pileup = float(pileupstring[0])
    print("Found a root file for pileup", pileup, "in file", file, "Objects: Clusters")

    rootfile = root.TFile.Open(file)
    rootfile.cd('BRIL_IT_Analysis/TEPX/Clusters')

    #build the histogram names
    histname = "Number of clusters for Disk "

    #loop the disks
    for disk in range(1,5):
        histminusz = root.gDirectory.Get(histname +"-"+str(disk))
        histplusz = root.gDirectory.Get(histname+str(disk))
        #add plus and minus Z histograms
        histminusz.Add(histplusz)

        #now loop the rings
        for ring in range(5):
            (m, sigma) = getParams(histminusz, ring)
            graphs[disk-1][ring].SetPoint(graphs[disk-1][ring].GetN(),pileup, m)
            graphs[disk-1][ring].SetPointError(graphs[disk-1][ring].GetN()-1,0, sigma)

    rootfile.Close()
    return

#get the stat error graph for clusters
def getStatErrorClusters(file,TEPXgraph, graphs=[]):
    ln4=16384
    trgkhz=0.075
    pileupstring = re.findall('summary_PU_(.*).root', file)
    pileup = float(pileupstring[0])
    print("Found a root file for pileup", pileup, "in file", file, "Objects: Clusters")

    rootfile = root.TFile.Open(file)
    rootfile.cd('BRIL_IT_Analysis/Clusters')

    #build the histogram names
    histname = "Number of clusters for Disk "
    nClustersTEPX=0

    #loop the disks
    for disk in range(1,5):
        histminusz = root.gDirectory.Get(histname +"-"+str(disk))
        histplusz = root.gDirectory.Get(histname+str(disk))
        #add plus and minus Z histograms
        histminusz.Add(histplusz)

        #now loop the rings
        for ring in range(5):
            (mean, sigma) = getParams(histminusz, ring)
            nClustersTEPX += mean
            if (pileup==0):
                value=0
            else:
                value=(math.sqrt(mean*ln4)/(mean*ln4))*(math.sqrt(trgkhz/40)/(trgkhz/40))
            graphs[disk-1][ring].SetPoint(graphs[disk-1][ring].GetN(),pileup, value)
            # graphs[disk-1][ring].SetMarkerStyle(8)
            # graphs[disk-1][ring].SetPointError(graphs[disk-1][ring].GetN()-1,0, sigma)
    
    if(pileup==0):
        TEPXvalue=0
    else:
        TEPXvalue=(math.sqrt(nClustersTEPX*ln4)/(nClustersTEPX*ln4))*(math.sqrt(trgkhz/40)/(trgkhz/40))
    TEPXgraph.SetPoint(TEPXgraph.GetN(),pileup,TEPXvalue)
    rootfile.Close()
    return

#get the linearity graph for hits
def getLinearityHits(file, graphs=[]):

    pileupstring = re.findall('summary_PU_(.*).root', file)
    pileup = float(pileupstring[0])
    print("Found a root file for pileup", pileup, "in file", file, "Objects: Hits")

    rootfile = root.TFile.Open(file)
    rootfile.cd('BRIL_IT_Analysis/TEPX/Hits')

    #build the histogram names
    histname = "Number of hits for Disk "

    #loop the disks
    for disk in range(1,5):
        histminusz = root.gDirectory.Get(histname +"-"+str(disk))
        histplusz = root.gDirectory.Get(histname+str(disk))
        #add plus and minus Z histograms
        histminusz.Add(histplusz)

        #now loop the rings
        for ring in range(5):
            (m, sigma) = getParams(histminusz, ring)
            graphs[disk-1][ring].SetPoint(graphs[disk-1][ring].GetN(),pileup, m)
            graphs[disk-1][ring].SetPointError(graphs[disk-1][ring].GetN()-1,0, sigma)

    rootfile.Close()
    return

# extract the mean/median, sigma and pileup from the folder in the rootfile - this is where the magic happens
def getLinearityCoincidences(file,nCoincidences, graphssum=[],graphsreal=[]):

    pileupstring = re.findall('summary_PU_(.*).root', file)
    pileup = float(pileupstring[0])
    print("Found a root file for pileup", pileup, "in file", file, "Objects:",nCoincidences,"Coincidences")

    rootfile = root.TFile.Open(file)
    rootfile.cd('BRIL_IT_Analysis/TEPX/'+str(nCoincidences)+'xCoincidences')

    #build the histogram names
    if nCoincidences == 2:
        histname = "Number of 2x Coincidences for Disk "
        realhistname = "Number of real 2x Coincidences for Disk "
    elif nCoincidences == 3:
        histname = "Number of 3x Coincidences for Disk "
        realhistname = "Number of real 3x Coincidences for Disk "

    #loop the disks
    for disk in range(1,5):
        histminusz = root.gDirectory.Get(histname +"-"+str(disk))
        histplusz = root.gDirectory.Get(histname+str(disk))

        realhistminusz = root.gDirectory.Get(realhistname +"-"+str(disk))
        realhistplusz = root.gDirectory.Get(realhistname +str(disk))

        #add plus and minus Z histograms
        histminusz.Add(histplusz)
        realhistminusz.Add(realhistplusz)


        #now loop the rings
        for ring in range(5):
            (msum, sigmasum) = getParams(histminusz, ring)
            (mreal, sigmareal) = getParams(realhistminusz, ring)

            graphssum[disk-1][ring].SetPoint(graphssum[disk-1][ring].GetN(),pileup, msum)
            graphssum[disk-1][ring].SetPointError(graphssum[disk-1][ring].GetN()-1,0, sigmasum)

            graphsreal[disk-1][ring].SetPoint(graphsreal[disk-1][ring].GetN(),pileup, mreal)
            graphsreal[disk-1][ring].SetPointError(graphsreal[disk-1][ring].GetN()-1,0, sigmareal)

    rootfile.Close()
    return

def extrapolateLinear(graph, basis=1):
    graph.Sort()
    #fit up to desired range, default is PU1
    pol1 = root.TF1("pol1",
                   "[0]*x+[1]",0,basis)
    graph.Fit("pol1","R")
    #draw the fit extrapolated to 200
    par0 = pol1.GetParameter(0)
    par1 = pol1.GetParameter(1)
    print(par0, par1)
    extrapolated = root.TF1("extra","[0]*x+[1]",0,200)
    extrapolated.SetParameter(0,par0)
    extrapolated.SetParameter(1,par1)
    extrapolated.SetLineColor(6)

    #draw a TGRaphErrors with the 1sigma confidence interval
    errors = root.TGraphErrors()
    for point in range(graph.GetN()):
        errors.SetPoint(point,graph.GetX()[point],0)
        #print("x:",graph.GetX()[point],"point:",point)

    root.TVirtualFitter.GetFitter().GetConfidenceIntervals(errors,0.95)
    errors.SetFillColor(6)
    errors.SetFillStyle(3003)

    return(extrapolated,errors)

def relativeNonlinearity(graph, fit):
    deviation = root.TGraph()
    deviation.SetTitle("Deviation from Linearity;PU;Diff[%]")
    deviation.SetMarkerStyle(8)
    graph.Sort()
    x=graph.GetX()
    y=graph.GetY()

    for point in range(graph.GetN()):
        diff = y[point]-fit.Eval(x[point])
        expect = fit.Eval(x[point])
        if expect == 0:
            relative_diff = 0
        else:
            relative_diff = diff/expect
        deviation.SetPoint(point, x[point],relative_diff*100)
        print("x:",x[point],"y:",y[point],"eval y:",expect,"point:",point,"diff:",relative_diff*100,"%")

    return deviation

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
    TEPXgraph=root.TGraph()
    TEPXgraph.SetMarkerStyle(8)
    TEPXgraph.SetTitle("Stat. Error TEPX [NB4]")

    for file in files:
        if file.find(".root"):
            filename = path+file
            #fill the actual graph for all available PU steps
            getStatErrorClusters(filename,TEPXgraph,graphs)
        else:
            print("Not a root file, skipping")
            continue

    rootfile = root.TFile("Results_StatError.root","RECREATE")
    index = 1
    for i in range(disks):
        for j in range(rings):
            #Cosmetics
            graphs[i][j].SetLineColor(1)
            graphs[i][j].SetTitle("Statistical Error [NB4] Disk"+str(i+1)+"Ring"+str(j+1)+";Pileup;Statistical Error [%]")
            graphs[i][j].SetMarkerStyle(8)
            c_canvas.cd(index)
            graphs[i][j].Draw("ap")

            #fit and extrapolate
            # (extrapolated[i][j],errors[i][j]) = extrapolateLinear(graphs[i][j],2)
            # errors[i][j].Draw("e3 same")
            # extrapolated[i][j].Draw("same")

            #calculate relative nonlinearity
            # deviation = relativeNonlinearity(graphs[i][j], extrapolated[i][j])
            # deviation.Write("Deviation Clusters Disk"+str(i+1)+"Ring"+str(j+1))

            #save canvases for the individual disk/ring combos
            savecanvas = root.TCanvas("Clusters Disk"+str(i+1)+"Ring"+str(j+1),"Clusters Disk"+str(i+1)+"Ring"+str(j+1))
            savecanvas.cd()
            graphs[i][j].GetYaxis().SetTitleOffset(1.2)
            graphs[i][j].GetXaxis().SetTitleSize(0.04)
            graphs[i][j].GetXaxis().SetLabelSize(0.04)
            graphs[i][j].GetYaxis().SetTitleSize(0.04)
            graphs[i][j].GetYaxis().SetLabelSize(0.04)
            graphs[i][j].Draw("ap")
            # errors[i][j].Draw("e3 same")
            # extrapolated[i][j].Draw("same")
            savecanvas.Write("Clusters Disk"+str(i+1)+"Ring"+str(j+1))
            index = index+1

    #Write out the summary as well
    c_canvas.Write("SummaryStatError")
    aCanvas = root.TCanvas("TEPX","TEPX")
    aCanvas.cd()
    TEPXgraph.Draw("ap")
    aCanvas.Write("StatErrorTEPX")
    rootfile.Close()

#second, let's only deal with the case for Hits
elif observable == "Hits":
    graphs = [[root.TGraphErrors() for j in range(rings)] for i in range(disks)]
    extrapolated = [[root.TF1() for j in range(rings)] for i in range(disks)]
    errors = [[root.TGraphErrors() for j in range(rings)] for i in range(disks)]

    for file in files:
        if file.find(".root"):
            filename = path+file
            #fill the actual graph for all available PU steps
            getLinearityHits(filename,graphs)
        else:
            print("Not a root file, skipping")
            continue

    rootfile = root.TFile("Results_Hits.root","RECREATE")
    index = 1
    for i in range(disks):
        for j in range(rings):
            #Cosmetics
            graphs[i][j].SetLineColor(1)
            graphs[i][j].SetTitle("Linearity Disk"+str(i+1)+"Ring"+str(j+1)+";Pileup;# of Hits")
            c_canvas.cd(index)
            graphs[i][j].Draw("ap")

            #fit and extrapolate
            (extrapolated[i][j],errors[i][j]) = extrapolateLinear(graphs[i][j],2)
            errors[i][j].Draw("e3 same")
            extrapolated[i][j].Draw("same")

            #calculate relative nonlinearity
            deviation = relativeNonlinearity(graphs[i][j], extrapolated[i][j])
            deviation.Write("Deviation Hits Disk"+str(i+1)+"Ring"+str(j+1))

            #save canvases for the individual disk/ring combos
            savecanvas = root.TCanvas("Hits Disk"+str(i+1)+"Ring"+str(j+1),"Hits Disk"+str(i+1)+"Ring"+str(j+1))
            savecanvas.cd()
            graphs[i][j].Draw("ap")
            errors[i][j].Draw("e3 same")
            extrapolated[i][j].Draw("same")
            savecanvas.Write("Hits Disk"+str(i+1)+"Ring"+str(j+1))
            index = index+1

    #Write out the summary as well
    c_canvas.Write("SummaryHits")
    rootfile.Close()

#for 2 and 3x coincidences
else:
    nCoincidences = int(observable[0])
    print("Working on", nCoincidences,"Coincidences")
    graphsreal = [[root.TGraphErrors() for j in range(rings)] for i in range(disks)]
    graphssum = [[root.TGraphErrors() for j in range(rings)] for i in range(disks)]
    extrapolated = [[root.TF1() for j in range(rings)] for i in range(disks)]
    errors = [[root.TGraphErrors() for j in range(rings)] for i in range(disks)]

    for file in files:
        if file.find(".root"):
            filename = path+file
            #fill the actual graph for all available PU steps
            getLinearityCoincidences(filename,nCoincidences,graphssum, graphsreal)
        else:
            print("Not a root file, skipping")
            continue

    rootfile = root.TFile("Results_Coincidences_"+str(nCoincidences)+"x.root","RECREATE")
    index = 1
    for i in range(disks):
        for j in range(rings):

            #Cosmetics
            graphssum[i][j].SetLineColor(4)
            graphsreal[i][j].SetLineColor(8)
            graphssum[i][j].SetMarkerColor(4)
            graphsreal[i][j].SetMarkerColor(8)
            graphssum[i][j].SetMarkerStyle(8)
            graphsreal[i][j].SetMarkerStyle(8)
            graphssum[i][j].SetTitle("Linearity Disk"+str(i+1)+"Ring"+str(j+1)+";<N_{PU}>;<Number of "+observable+" Coincidences")
            c_canvas.cd(index)
            graphssum[i][j].Draw("ap")
            graphsreal[i][j].Draw("p same")

            #fit and extrapolate
            (extrapolated[i][j],errors[i][j]) = extrapolateLinear(graphssum[i][j],2)
            errors[i][j].Draw("e3 same")
            extrapolated[i][j].Draw("same")

            #calculate relative nonlinearity
            deviation = relativeNonlinearity(graphssum[i][j], extrapolated[i][j])
            deviation.Write("Deviation "+observable+" Disk"+str(i+1)+"Ring"+str(j+1))

            #save canvases for the individual disk/ring combos
            savecanvas = root.TCanvas(observable+"Coincidences Disk"+str(i+1)+"Ring"+str(j+1),observable+"Coincidences Disk"+str(i+1)+"Ring"+str(j+1))
            savecanvas.cd()
            graphssum[i][j].Draw("ap")
            graphsreal[i][j].Draw("p same")
            errors[i][j].Draw("e3 same")
            extrapolated[i][j].Draw("same")
            savecanvas.Write(observable+"Coincidences Disk"+str(i+1)+"Ring"+str(j+1))
            index = index+1

    #Write out the summary as well
    c_canvas.Write("SummaryClusters")
    rootfile.Close()


# if __name__ == '__main__':
        # main()
