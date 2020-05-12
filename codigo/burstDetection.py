import scipy as sc
import scipy.io as io
import matplotlib.pyplot as gr
import csv
import time
import os
import string
import networkx as nx

def slopeAB(A,B):
    """
    slopeAB calculates the slope between two points A and B
    Example:
    sAB=slopeAB(a,b)
    """
    mAB = (B[1]-A[1]) / (B[0]-A[0])
    return mAB


def calcFwdSlopes(x,y): 
    """
    calcFwdSlopes calculates the forward-slopes between pairs of points from graph {(x,y)}
    Example:
    s=calcFwdSlopes(x,y)
    """
    nY=len(y)
    s=list()
    for m in sc.arange(nY-1):
        A=[x[m], y[m]]
        s.append(list())
        for n in sc.arange(m+1,nY):
            B=[x[n], y[n]]
            s[m].append(slopeAB(A,B))
    return s 

def fwdVisibilityInds(x,y):
    """
    fwdVisibilityInds obtains the indices of the forward-visible neighbours within a graph {(x,y)}
    neighs=fwdVisibilityInds(x,y)
    """
    nodes=list()
    neighs=list()
    nY=len(y)
    s=calcFwdSlopes(x,y)
    #nodes = sc.arange(nY)
    for m in sc.arange(nY-2):
        a=[x[m], y[m]]
        b=[x[m+1], y[m+1]]
        sAB=slopeAB(a,b)
        neighs.append(list())
        neighs[m].append(m+1)
        for n in sc.arange(m+2,len(s[m])):
            if sAB < s[m][n]:
                sAB=s[m][n]
                neighs[m].append(n)
    return neighs

def addFwdNeighs(g,neighInds):
    """
    addFwdNeighs adds new edges to the graph g from a list of neighbour indices in neighInds
    Example:
    g1=addFwdNeighs(g1,n1)

    """
    n=0
    for neigh in neighInds:
        for j in neigh:
            g.add_edge(g.nodes()[n], g.nodes()[j])
        n=n+1
    return g

def visibilityGraph(xVals,yVals):
    """
    Create a visibility graph from a series of sample points stored in the arrays xVals and yVals.
    """
    neighs=fwdVisibilityInds(x=xVals,y=yVals)
    g=nx.Graph()
    g.add_nodes_from(xVals)
    g=addFwdNeighs(g,neighs)
    g.remove_node(g.nodes()[-1])
    return g


def visibilityFromSpikeTrains(trains, ifrs, a=0, b=-1):
    """
    createVisibilityGraphs takes two lists with spike trains and their instantaneous firing rates and calculates visibility graphs
    """
    # Note the length of xList must be equal to the length of yList
    g=list()
    nSeries= len(trains)
    for n in sc.arange(nSeries):
        g.append(visibilityGraph(trains[n][a:b],ifrs[n][a:b]))
    return g


# Function that plots the graphs form a list into a set of axes (also specified from a list). 

def showGraphList(gList, gLabels, layout="spring", nodeSize=50, nodeColor1='blue',transp=0.5):
    nAx = len(gList)
    f=gr.figure(figsize=(17,11))
    rows=2; cols= sc.ceil(nAx/rows)
    ax=list(); gr.ioff();
    for n in sc.arange(nAx):
        ax.append(f.add_subplot(rows,cols,n+1))
        if layout=="spring":
            nx.draw(gList[n],ax=ax[n],pos=nx.spring_layout(gList[n]), node_size=nodeSize,node_color=nodeColor1,alpha=transp)
        if layout=="shel":
            nx.draw(gList[n],ax=ax[n],pos=nx.shell_layout(gList[n]), node_size=nodeSize,node_color=nodeColor1,alpha=transp)
        if layout=="random":
            nx.draw(gList[n],ax=ax[n],pos=nx.random_layout(gList[n]), node_size=nodeSize,node_color=nodeColor1,alpha=transp)
        if layout=="circular":
            nx.draw(gList[n],ax=ax[n],pos=nx.circular_layout(gList[n]), node_size=nodeSize,node_color=nodeColor1,alpha=transp)
        if layout=="spectral":
            nx.draw(gList[n],ax=ax[n],pos=nx.spectral_layout(gList[n]), node_size=nodeSize,node_color=nodeColor1,alpha=transp)
        ax[n].set_title(gLabels[n])
    gr.ion(); gr.draw()
    return f

# -----------------------------------------------------------------
# Trenes y estimulos
# -----------------------------------------------------------------
# Plot raster from trains
def plotRaster(ax, trains):
    nTrains=len(trains)
    for n in sc.arange(nTrains):
        train= trains[n]
        y = n * sc.ones(len(train))
        ax.plot(trains[n], y, "|")
    return ax

# --------------------------------------------------
# Spike trains with gamma interspike interval distributions
# --------------------------------------------------
def createNGammaTrainsNPulses(nPulses=500, nTrains=60, meanRate=10.0, msRefrac=3.0, graph=0.0):
    """
    Create N spike trains with gamma interspike interval distribution
    Example:
    tr, isis, ifrs= createNGammaTrains(nPulses=500, nTrains=60, meanRate=20.0,graph=0.0)
    """
    isis=sc.random.gamma(msRefrac,1000.0/meanRate,(nTrains,nPulses))
    trains= isis.cumsum(1)/1000.0
    ifrs = 1000.0/isis
    if graph:
        gr.figure()
        gr.ioff()
        for nn in sc.arange(nTrains):
            gr.plot(trains[nn], nn*sc.ones(nPulses), '|')
        gr.ion()
        gr.draw()
    return trains, isis, ifrs

# ............................
# Needs fixing
# ............................
def createNGammaTrainsMaxTime(maxTimeSecs=1.0, nTrains=10, meanRate=10.0, msRefrac=3.0, graph=0.0):
    """ Example:
    tr, isis, ifrs= createNGammaTrains(maxTime=500.0, nTrains=30, meanRate=20.0,graph=0.0)
    """
    kHz = 1000/meanRate
    maxNSpikes = 2*maxTimeSecs*kHz
    print(maxNSpikes)
    isiA=sc.random.gamma(msRefrac, kHz, (nTrains,maxNSpikes))
    trains=isiA.cumsum(1)/1000.0
    print(trains, trains.shape)
    spikeTrains=list()
    isis=list()
    ifrs=list()
    for n in sc.arange(nTrains):
        print(trains[n])
        keepInds = sc.where(trains[n]<maxTimeSecs)[1]
        spikeTrains.append(trains[n][keepInds])
        isis.append(isiA[keepInds])
        ifrs.append(1000.0/isis[n])
    print(keepInds,spikeTrains)
    if graph:
        gr.figure()
        gr.ioff()
        for nn in sc.arange(nTrains):
            gr.plot(trains[nn], nn*sc.ones(len(trains[nn])), '|')
        gr.ion()
        gr.draw()
    return trains, isis, ifrs

def createBurstingTrain(maxTime, burstRate=10.0, interBurstRate=30.0, nonBurstProp=0.1):
    nSpikes= sc.int32(maxTime/burstRate)
    nBursts = (1-nonBurstProp) * nSpikes 
    b, ibis, ibrs= createNGammaTrains(nPulses=nSpikes, nTrains=1, meanRate=burstRate,graph=0.0)
    
    for n in sc.arange(nBursts):
        tr, isi, ifr= createNGammaTrains(nPulses=nSpikes, nTrains=1, meanRate=burstRate,graph=0.0)
        
    return train

def calculateISIs(train):
    y=sc.zeros(len(train))
    y[1:]= train[1:]- train[:-1]
    return y


# ---------------------------
# Cross-correlations. No resolution.
# ---------------------------
def calcMultipleSpikeTrainCCs(trains):
    """
    Calculate the cross-correlograms of multiple spike trains. The results are put in an array of nx (n+1)/2 elements. Example:
    ccs, taus= calcMultipleSpikeTrainCCs(trains)
    """
    ccs = list()
    taus = list()
    nTrains = len(trains)
    for m in sc.arange(nTrains):
        for n in sc.arange(nTrains):
            cc=sc.correlate(trains[m],trains[n], mode='full')
            ccs.append(cc)
            nPts=len(cc)
            tau=sc.arange(-nPts/2.0,nPts/2.0)
            taus.append(tau)
    return ccs,taus

def alphaFunction(x, A=1.0, tau=1.0, downAccel=1.0):
    aa= sc.int32(x>0)
    xovertau = x/tau
    return A* xovertau * sc.exp( downAccel*(1 - xovertau))

def trainAlpha(samplingTimes, pulseTimes, tauTrain=1.0,downAccel=1.0):
    """
    Obtain a train of stimuli centered around a spike train. each stimulus is a
    a= trainAlpha(samplingTimes, pulseTimes, tauTrain=1.0,downAccel=1.0)
    """
    nPts= len(samplingTimes)
    train = sc.zeros(nPts)
    alpha=alphaFunction(samplingTimes,A=1.0,tau=tauTrain,downAccel=downAccel)
    for n in range(len(pulseTimes)):
        #nn=gr.find(samplingTimes<pulseTimes[n]).max()
        nn=sc.where(samplingTimes<pulseTimes[n])[1].max()
        train[nn:] = train[nn:] + alpha[0: nPts-nn]
    return train


# -----------------------------------------------------------------
# I/O tools to get data from different kinds of files
# -----------------------------------------------------------------
# -----------------------------------------------------------------
# Make a list of files with a common string in their name
# -----------------------------------------------------------------
def getFileList(dataDir, prefix, suffix, showNames=0):
    """
    Example:
    files=getFileList(dataDir, prefix, suffix)
    """
    files=list()
    for f in os.listdir(dataDir):
        a= sc.int32(os.path.isfile(os.path.join(dataDir,f)))
        b= sc.int32(str.find(f,prefix)>-1)
        c= sc.int32(str.find(f,suffix)>0)
        #print(a,b,c)
        if (a*b*c):
            files.append(f)
    nFiles = len(files)
    if showNames:
        print("Found %d files with the indicated string"%nFiles)
        print(files)
    return files


def preprocessTrainsFromMatFiles(fileNameList):
    """
    preprocessTrainsFromFiles extracts spike trains from a list of .mat files containing spike times, one train of time-stamps per file.
    The spike trains are preprocessed by calculating interspike intervals and instantaneous firing rates. 
    The output is a dictionary containing the preprocessed data. 
    Example of use:
    ctrFiles=getFileList(dataDir, prefix="CntSt", suffix=".mat")
    ctrData= preprocessTrainsFromMatFiles(ctrFiles)
    """
    dataDir='./trenesBursts/aleSpikeTimes/'
    spikes=list(); isis=list(); ifrs=list()
    n=0
    for f in fileNameList:
        ff=dataDir + f
        s=sc.squeeze(io.loadmat(ff)["savevar"][0][0])
        spikes.append(s)
        print("Found %d spikes in %s"%(len(s),ff))
        isi = sc.zeros(len(s)); ifr = sc.zeros(len(s)) #Note time stamps are in milliseconds
        isi[1:]= s[1:]-s[:-1]
        ifr[1:]= 1/isi[1:]
        isis.append(isi); ifrs.append(ifr)
        n=n+1
    return {"spikeTrains":spikes, "ISIs":isis, "IFRs":ifrs}


def readMatlabArray(dataDir,files,varName):
    nFiles=len(files)
    data= list()
    for n in sc.arange(nFiles):
        fName = files[n]
        pigas = sc.squeeze(io.loadmat(dataDir + fName)[varName][0][0])
        data.append(pigas)
    return data

def readColumnCSV(dataDir='./',  fileName ='Neuron.csv', delimiter=',', nHeaderRows=0):
    data=dict()
    data['name'] = fileName
    with open(dataDir+fileName,'rb') as f:
        reader=csv.reader(f, delimiter=delimiter)
        dataList=list()
        nRows=0
        for row in reader:
            if nRows>nHeaderRows-1:
                aa=sc.float32(row)
                dataList.append(aa)
            nRows=nRows+1
            
        #for row in reader:
        #    aa=sc.float32(row)
        #    dataList.append(aa)
        #    nRows=nRows+1
            
        data['values']= sc.array(dataList).squeeze()
        data['nRows']= nRows
        f.close()
    return data

def readMultipleColumnsCSV(delimiter=',', fileName ='.csv'):
    data=dict()
    data['name'] = fileName
    with open(dataDir+fileName,'rb') as f:
        reader=csv.reader(f, delimiter=delimiter)
        dataList=list()
        nPts=0
        for row in reader:
            aa=sc.float64(row)
            dataList.append(aa)
            nPts=nPts+1
            
        data['values']= sc.array(dataList).squeeze()
        data['nPts']= nPts
        f.close()
    return data

def entropy(relFreq):
    goodInds= relFreq.nonzero()[0]
    r = relFreq[goodInds]
    h = -(r * sc.log2(r)).sum()
    return h


# Calculate a cumulative distribution function from a sample, using a specified set of bins. 
def calcCDF(sample,graph=0):
    """Example:
    Calculate a cumulative distribution function from a sample, using a specified set of bins. 
    sSize= 100
    mySample=sc.randn(sSize)
    Cdf = calcCDF(bins=myBins, sample=mySample)
    """
    mySample=sc.copy(sample)
    mySample.sort()
    bins= sc.unique(mySample)
    nBins = len(bins)
    sSize= sc.float32(len(mySample))
    #print(bins,nBins)
    #print(mySample,sSize)
    cdf = sc.zeros(nBins)
    for n in sc.arange(nBins):
        A= sc.where(mySample<bins[n])[0]
        cdf[n]=len(A)/sSize
    #print(cdf,cdf.sum())
    f = sc.vectorize(lambda x : sc.interp(x, bins,cdf))
    fInverse = sc.vectorize(lambda y: sc.interp(y,cdf, bins))
    if graph:
        gr.figure()
        gr.plot(bins,cdf,'wo')
        #gr.plot(sample[30:-20],f(sample[30:-20]),'.')
    return f,fInverse


def calcThresholdsFromCDF(cdfInverse, quants=(0.95,)):
    """
    Calculate thresholds using the inverse of the cumulative distribution function given a list of quants. 
    Examples:
    calcThresholdsFromCDF(cdfInverse, quants=(0.90,))
    calcThresholdsFromCDF(cdfInverse, quants=0.95)
    """
    if len(quants)>0:
        nValues = len(quants)
        alphaX=sc.zeros(nValues)
        for n in sc.arange(nValues):
            alphaX[n]= cdfInverse(quants[n])
    else: 
        alphaX= cdfInverse(quants)
    return alphaX

def calcISI(train):
    y= sc.zeros(len(train))
    y[1:]= train[1:]-train[:-1]    
    return i

def calcIFR(train):
    ifr= sc.zeros(len(train))
    ifr[1:]= 1/(train[1:]-train[:-1])
    return ifr

def calcISICDF(train):
    y= calculateISIs(train)
    f,fInverse=calcCDF(y,graph=0)
    return f,fInverse

def calcIFRCDF(train):
    y= calculateIFRs(train)
    f,fInverse=calcCDF(y,graph=0)
    return f,fInverse

def findBurstsInTrainIFR(train, alpha=0.95):
    """
    Find bursts of high activity (high instantaneous firing rate) in a train of pulses with respect to a statistical threshold upQuant. Example:
    bursts= findBurstsInTrainIFR(train, alpha=0.05)
    """
    nSpikes=len(train)
    ifr= sc.zeros(nSpikes)
    isi= sc.zeros(nSpikes)
    dIFR= sc.zeros(nSpikes)
    isi[1:]= train[1:]-train[:-1]
    #dISI[1:] = isi[1:] /  
    ifr[1:]= 1/(train[1:]-train[:-1])
    dIFR[1:]= (ifr[1:]-ifr[:-1])/isi[1:]
    ifrCDF,ifrCDFInverse= calcCDF(ifr,graph=0)
    isiCDF,isiCDFInverse= calcCDF(isi,graph=0)
    # Find spikes during high activity
    ifrUpThresh= calcThresholdsFromCDF(ifrCDFInverse, (alpha,))
    ifrDnThresh= calcThresholdsFromCDF(ifrCDFInverse, (1-alpha,))
    isiThresh= calcThresholdsFromCDF(isiCDFInverse, (alpha,))
    #rHighInds= sc.where( ifr>ifrThresh)[0]
    rHighInds= sc.where( ifr>ifrUpThresh)[0]
    lHighInds= rHighInds-1
    highInds= sc.union1d(rHighInds,lHighInds)
    highSpikeTimes= train[highInds]
    #lowSpikeTimes= train[lowInds]
    aa= sc.zeros(len(highInds))
    aa[1:]= sc.diff(highInds)
    #bb= sc.zeros(len(lowInds))
    #bb[1:]= sc.diff(lowInds)
    startInds=sc.where(aa!=1)[0]
    burstStarts= highSpikeTimes[startInds]
    burstEnds= highSpikeTimes[startInds-1]
    nBursts= len(burstStarts)
    burstStarts=burstStarts[:-1]
    burstEnds=burstEnds[1:]    
    pBurst = sc.float32(nBursts)/nSpikes
    burstDurs = burstEnds-burstStarts
    bursts={"train": train, 
        "highInds":highInds, "highSpikeTimes":highSpikeTimes, 
        "burstStarts":burstStarts, "burstEnds":burstEnds, 
        "nBursts":nBursts, "nSpikes": nSpikes, "pBurst": pBurst, "burstDurs": burstDurs,
        "ifr":ifr, "ifrCDF":ifrCDF, "ifrCDFInverse":ifrCDFInverse,"alpha":alpha,"ifrThresh":ifrThresh,
        "isi":isi, "isiCDF":isiCDF, "isiCDFInverse":isiCDFInverse,"isiThresh":isiThresh,
        "dIFR":dIFR}
    return bursts


# --------------------------------------------------
# Detection of bursts in neurons using v and dv/dt
# --------------------------------------------------
def voltageBurstDetection(volts, vBeta=0.99, dvdtBeta=0.99, samplingPeriod= 100e-3):
    """Detects action potentials from an array of voltages and returns the indices where action potentials occur. Example:
    spikeInds= neuronBurstDetection(volts, alpha=0.01, samplingPeriod= 100e-3)
    """
    dvdt= sc.zeros(len(volts))
    dvdt[1:] = (volts[1:] - volts[0:-1])/samplingPeriod
    vCDF, vCDFInverse=calcCDF(volts)
    dvdtCDF, dvdtCDFInverse=calcCDF(dvdt)
    vThresh=dvdtCDFInverse(vBeta)
    dvdtThresh=dvdtCDFInverse(dvdtBeta)
    spikeInds = sc.where( (dvdt>dvdtThresh)&(volts>vThresh))[1]
    spikeInds2= sc.hstack([0, 1+ sc.where((spikeInds[1:]- spikeInds[:-1])>1)[1] ])
    spikeIndices =spikeInds[spikeInds2]
    bd={'v_cdf': vCDF, 'v_cdfInv': vCDFInverse, 'dvdt_cdf': dvdtCDF, 'dvdt_cdfInv': dvdtCDFInverse, 'dvdt': dvdt, 'spikeInds':spikeIndices,'nSpikes':len(spikeIndices),'vThreshold': vThresh, 'dvdtThreshold':dvdtThresh, 'vBeta': vBeta, 'dvdtBeta':dvdtBeta }
    return bd





# --------------------------------------------------
# Detection of bursts
# --------------------------------------------------
if 0:
    statThresh=0.99
    dataDir='./trenesBursts/'
    data= readColumnCSV()
    v= data['values'].squeeze()
    samplingStep= 25e-3
    sampTimes=sc.arange(0, len(v)*samplingStep,samplingStep)/1000.0
    bursts= voltageBurstDetection(vBeta=statThresh, dvdtBeta=statThresh, samplingPeriod= 100e-3, volts= v)
    spikeTimes=sampTimes[bursts['spikeInds']]
    dvdt=bursts['dvdt']
    dvdtArray=sc.arange(dvdt.min(), dvdt.max(), 0.1)
    thresh=bursts['dvdtThreshold']

# --------------------------------------------------
# Bursts triggering Ca-like signals
# --------------------------------------------------
if 0:
    nTrains= 60
    tr, isis, ifrs= createNGammaTrains(nPulses=100, nTrains=nTrains, meanRate=10.0, graph=0.0)
    sampTimes= sc.arange(0,tr.max(), 1e-3)

    transients= list()
    for n in sc.arange(nTrains):
        a= trainAlpha(samplingTimes=sampTimes, pulseTimes=tr[n], tauTrain=0.01,downAccel=0.2) + sc.rand(len(sampTimes))
        a = (a- a.min())/ (a.max()- a.min())
        transients.append(a)

    figTrains= gr.figure(figsize=(17,11))
    rows=1; cols=1;
    gr.ioff()
    ax=figTrains.add_subplot(rows,cols,1)
    for n in sc.arange(nTrains):
        ax.plot(tr[n], sc.ones(len(tr[n]))*(n+1),'|')
    ax.set_xlabel('seconds')
    ax.set_ylabel('neurons')
    gr.ion(); gr.draw()

    figTrans= gr.figure(figsize=(17,11))
    rows=1; cols=1;
    gr.ioff()
    ax=figTrans.add_subplot(rows,cols,1)
    for n in sc.arange(nTrains):
        ax.plot(sampTimes, n+transients[n], alpha=0.4)
    ax.set_xlabel('seconds')
    ax.set_ylabel('neurons')
    gr.ion(); gr.draw()


#
if 0:
    figBurst= gr.figure(figsize=(17,11))
    rows=1; cols=3; ax=list()
    gr.ioff()
    ax.append(figBurst.add_subplot(2,1,1))
    ax.append(figBurst.add_subplot(2,2,3))
    ax.append(figBurst.add_subplot(2,2,4))
    ax[0].plot(sampTimes,v, 'k', ms=1, alpha=1, label=r'$(t,v(t))$')
    ax[0].plot(spikeTimes, -90*sc.ones(len(spikeTimes)),'k|',ms=10)
    ax[1].plot(v, dvdt,'k.', label=r'$(v, \partial_t v)$')
    ax[2].plot([0,1], thresh*sc.ones(2), 'k', ms=3, label=r'(%g,%g)'%(statThresh, thresh) )
    ax[2].plot( bursts['dvdt_cdf'](dvdtArray), dvdtArray, 'wo', label=r'$(f_{\partial_t v}, \partial_t v)$')
    ax[0].set_xlim(sampTimes.max()/3.0, sampTimes.max()*2/3)
    for n in sc.arange(rows*cols):
        ax[n].legend()

    gr.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.98,
                       wspace=0.1, hspace=0.1)
    gr.ion(); gr.draw()



if 0: 
    fileName='T-Esther.csv'
    data=readColumnCSV(dataDir='./', delimiter=',', fileName =fileName, nHeaderRows=4)
    tStamps= data['values'].transpose()
    s1 = tStamps[10]
    s2 = tStamps[30]
    cc=sc.correlate(s1,s2, mode='full')
    nPts=len(cc)
    tau = sc.arange(-nPts/2, nPts/2)
    figg= gr.figure()
    gr.plot(tau,cc)


# Example: Create a spike train, threshold it, show the thresholds
if 0: 
    tr, isis, ifrs= createNGammaTrains(nPulses=500, nTrains=1, graph=0.0)
    train=tr[0];ifr=ifrs[0]; isi=isis[0]
    myBins= sc.unique(sc.sort(ifr))
    cdf,cdfInverse = calcCDF(bins=myBins, sample=ifr)
    alphaX= calcThresholdsFromCDF(cdf=cdf, quantValues=[0.1,0.9])
    fig=gr.figure()
    ax1= fig.add_subplot(121)
    ax2= fig.add_subplot(122)
    gr.ioff()
    ax1.plot(train, ifr.max()*sc.ones(len(train)),'|')
    ax1.plot(train, ifr, 'o')
    ax1.plot([train.min(), train.max()], [alphaX[0],alphaX[0]], 'r')
    ax1.plot([train.min(), train.max()], [alphaX[1],alphaX[1]], 'b')
    ax1.set_xlabel('secs')
    ax1.set_ylabel('Hz')
    ax2.plot([0,1], [alphaX[0],alphaX[0]],'r')
    ax2.plot([0,1], [alphaX[1],alphaX[1]],'b')
    ax2.plot(cdf, myBins)
    ax1.set_xlabel('cdf')
    gr.ion(); gr.draw()




if 0: 
    # Now create a coded representation for an ensemble: -1 for significantly silent, 1 for significantly firing, 0 for unsignifcantly active.
    nNeurons=20; nPulsos=500
    trains, isis, ifrs= createNGammaTrains(nPulses=nPulsos, nTrains=nNeurons, graph=0.0)
    myBins= sc.arange(ifrs.max())
    thrsDn,thrsUp=calcDnUpEnsembleThresholds(ensembleSample=ifrs, bins=myBins, dn=0.1, up=0.9)
    # Create a code based on low or high firing
    windowSize= 0.075
    stepSize=0.1
    nSteps= 1+ trains.max()/stepSize
    newN= sc.zeros((nNeurons,nPulsos))
    indsDn=list(); indsUp=list()
    for mm in sc.arange(nNeurons):
        upInds=sc.where(ifrs[mm]> thrsUp[mm])[1]
        dnInds=sc.where(ifrs[mm]< thrsDn[mm])[1]
        nUp = len(upInds)
        nDn = len(dnInds)
        if nUp:
            newN[mm][upInds]= 1
        if nDn:
            newN[mm][dnInds]= -1
        indsDn.append(dnInds)
        indsUp.append(upInds)

    f3=gr.figure();
    gr.ioff(); r=1;c=2
    ax1=f3.add_subplot(r,c,1)
    ax2=f3.add_subplot(r,c,2)
    for mm in sc.arange(nNeurons):
        up=indsUp[mm]
        dn=indsDn[mm]
        ax1.plot(trains[mm][dn], mm*sc.ones(len(dn)),'ro')
        ax1.plot(trains[mm][up], mm*sc.ones(len(up)),'go')
        # Show the thresholds for the different neurons
        ax2.plot(thsDn,'bo');
        ax2.plot(thsUp,'ro');
        ax2.set_xlabel('neuron');
        ax2.set_xlabel('Threshold (Hz)')

    gr.ion(); gr.draw()


if 0: 
    pref="Cnt"
    cell= 1
    st= trains[pref][cell]
    print(st.fName)
    print("% bursts", 100*len(st.bInds)/sc.float64(len(st.spikeTimes)))

    print(st.bStarts)
    print(st.bEnds)
    print(st.bEnds-st.bStarts)






# ----------------------------------------
# New calcThresholdsFromCDF...
# ----------------------------------------
if 0: 
    def calcThresholdsFromCDF(cdfInverse, quants):
        """
        Calculate thresholds using the inverse of the cumulative distribution function given a list of quants. 
        Examples:
        calcThresholdsFromCDF(cdfInverse, quants=(0.90,))
        calcThresholdsFromCDF(cdfInverse, quants=0.95)
        """
        if len(quants)>1:
            nValues = len(quants)
            alphaX=list()
            for a in quants:
                alphaX.append(cdfInverse(a))
        else: 
            alphaX= cdfInverse(quants)
        return alphaX

    #
    threshs=sc.arange(0.1,1.0001,0.1)
    cell=0
    train=stAnt[cell]
    if 1:
        bp=dict()
        nSpikes=len(train)
        isi= sc.zeros(nSpikes); 
        dIFR= sc.zeros(nSpikes)
        isi[1:]= train[1:]-train[:-1]
        ifr=1/isi
        dIFR_Right= (ifr[1:]-ifr[:-1])/isi[1:] 
        dIFR_Left= (ifr[1:]-ifr[:-1])/isi[:-1]
        dIFR[1:]= (dIFR_Right+dIFR_Left)/2.0 
        cdf,cdfInverse= calcCDF(dIFR,graph=0)
        dIFRThresh= sc.squeeze(calcThresholdsFromCDF(cdfInverse, threshs))
    #
    f=gr.figure()
    ax=list()
    rows=2
    for r in sc.arange(rows):
        ax.append(f.add_subplot(rows,1,r+1))
    gr.ioff()
    ax[0].plot(train,ifr,'bo', label=r"avg")
    ax[0].plot(train,dIFR, "ko", label=r"avg")
    ax[1].plot(train[:-1],dIFR_Right,alpha=0.3)
    ax[1].plot(train[1:],dIFR_Left,alpha=0.3)
    ax[1].plot(train,dIFR, "k", label=r"avg")
    for r in sc.arange(rows):
        ax[r].legend()
        ax[r].set_xlim(150,300)
    gr.ion(); gr.draw()

