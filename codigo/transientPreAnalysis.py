# Suite of functions to pre-analyze and process data sets obtained from simultaneously recorded calcium transients 
# Created: MAHV 20150801
# Last modified: MAHV 20160216. Extraction from IGOR files. 

import scipy as sc
from scipy import fftpack
import os
import matplotlib.pylab as gr
import igor.igorpy as igor
import string as st
import pickle

# Frequency-dependent smoothing
# FFT low-high cut
def lowHighCutFreqSmoothing(x,lowCut=400,hiCut=600):
	rft= fftpack.rfft(x); 
	rft[:lowCut]=0; 
	rft[hiCut:]=0; 
	#y_smooth = sc.ifft(rft, N)
	y=fftpack.irfft(rft); 
	return y #,y_smooth

# Approximate the derivative from a signal
def calcTwoSideSecantSlope(t,signal):
    nPts=len(f)
    dfFwd= sc.zeros(nPts)
    dfBwd= sc.zeros(nPts)
    df= sc.zeros(nPts)
    dfFwd[1:] = (signal[1:]-signal[:-1])/(t[1:]-t[:-1])
    dfBwd[:-1] = dfFwd[1:]
    df = (dfFwd + dfBwd )/2.0
    return df


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
    fInverse = sc.vectorize(lambda y: sc.interp(y,cdf,bins))
    if graph:
        gr.figure()
        gr.plot(bins,cdf,'wo')
        #gr.plot(sample[30:-20],f(sample[30:-20]),'.')
    return f,fInverse


def calcThresholdFromCDF(cdfInverse, quant=0.95):
    """
    Calculate thresholds using the inverse of the cumulative distribution function given a list of quants. 
    Example:
    alphaX=calcThresholdFromCDF(cdfInverse, quant=0.95)
    """
    alphaX= cdfInverse(quant)
    return alphaX

def calcThresholdsFromCDF(cdfInverse, quants=(0.95,)):
    """
    Calculate thresholds using the inverse of the cumulative distribution function given a list of quants. 
    """
    if len(quants)>1:
        nValues = len(quants)
        alphaX=sc.zeros(nValues)
        for n in sc.arange(nValues):
            alphaX[n]= cdfInverse(quants[n])
    else: 
        print('calcThresholdsFromCDF needs more than one quantile to work')
    return alphaX


def calcISICDF(train):
    y= calculateISIs(train)
    f,fInverse=calcCDF(y,graph=0)
    return f,fInverse


def fMinNorm(dataArray):
    """
    Nomalization of a time series by taking the minimum as a baseline. 
    Example:
    fmd = fMinNorm(dataArray)     
    where dataArray is a one-dimensional time series.  
    """
    fMinData=sc.zeros(dataArray.shape,'float64')
    for n in sc.arange(nCells):
        fMin= dataArray[n].min()
        fMinData[n]= (dataArray[n] - fMin)/(fMin)
    return fMinData

def fMaxMinNorm(dataArray):
    """
    Nomalization of a time series by taking the minimum as a baseline. 
    Example:
    fmd = fMinNorm(dataArray)     
    where dataArray is a one-dimensional time series.  
    """
    fMinData=sc.zeros(dataArray.shape,'float64')
    for n in sc.arange(nCells):
        fMin= dataArray[n].min()
        fMax= dataArray[n].max()
        fMinData[n]= (dataArray[n] - fMin)/(fMax-fMin)
    return fMinData


def lowHighCutFreqSmoothing(x,lowCut=400,hiCut=600):
    """
    Cut low and high frequency oscilations from a signal x. Example:
    lowHighCutFreqSmoothing(x,lowCut=400,hiCut=600)
    """
    rft= fftpack.rfft(x); 
    rft[:lowCut]=0; 
    rft[hiCut:]=0; 
    #y_smooth = sc.ifft(rft, N)
    y=fftpack.irfft(rft); 
    return y #,y_smooth

# --------------------------------------------------
# Data extraction from IGOR files
# --------------------------------------------------

def getFileList(dataDir, prefix, suffix, includeDataDir=1):
    """
    getFileList look for files with specific prefix and suffix within the directory dataDir
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
            if includeDataDir:
                files.append(dataDir+f)
            else:
                files.append(f)
    nFiles = len(files)
    print("Found %d files with the indicated string"%nFiles)
    print(files)
    return files


def extractPXPData(fName):
    """
    extractPXPData
    Example:
    allData=extractPXPData("microcircuitsNetworks/cort76dp1c.pxp")
    """
    allData=igor.load(fName)
    # Dictionary containing the names of recorded variables during the experiment
    bbb=allData.children[0].userstr[b"S_waveNames"]
    aaa=bbb.decode("UTF-8")
    dataNames= aaa.split(";")[:-1]
    extractedData1=list()
    extractedData2=list()
    for nam in dataNames:
        if len(nam)>3:
            myInd=sc.int32(nam[4:]) -1
            str2= "w2=allData." + nam + ".data"%(myInd)
            exec(str2)
            extractedData2.append(w2)
        else:
            print("Found empty string")

    rawData=sc.array(extractedData2)
    return rawData


def extractDebleachedFminData(filePath):
    """
    extractDebleachedFminData takes a path to a pxp data file as only argument and extracts all the wave data found there, 
    provided the data has the suffix F_min
    
    Example:
    dataDir="./microcircuitsNetworks/"
    fileName= "cort76dp1c.pxp"
    waveData, timeStamps=extractDebleachedFminData(dataDir+fileName)
    """
    allData= igor.load(filePath)
    #dataNames= st.digits.split(allData.children[0].userstr["S_waveNames"], ";")[:-1]
    bbb=allData.children[0].userstr[b"S_waveNames"]
    aaa=bbb.decode("UTF-8")
    dataNames= aaa.split(";")[:-1]
    waves= list()
    for m in sc.arange(len(dataNames)):
        waveNum=sc.int32(dataNames[m][1+str.rfind(dataNames[0],"e"):])
        str1= "waves.append(allData."+dataNames[m] + "_%dF_min.data)"%(waveNum-1)
        #print(dataNames[m],str1)
        exec(str1)
    return sc.array(waves), sc.array(allData.sec.data)



def extractDebleachedData(filePath):
    allData= igor.load(filePath)
    bbb=allData.children[0].userstr[b"S_waveNames"]
    aaa=bbb.decode("UTF-8")
    dataNames= aaa.split(";")[:-1]
    waves= list()
    for m in sc.arange(len(dataNames)):
        waveNum=sc.int32(dataNames[m][1+str.rfind(dataNames[0],"e"):])
        str1= "waves.append(allData."+dataNames[m] + "_%d.data)"%(waveNum-1)
        #print(dataNames[m],str1)
        exec(str1)
    return sc.array(waves)




# --------------------------------------------------
# Calculate activity packets from calcium rasters
# NEEDS WORK.... 
# --------------------------------------------------
def calcDnUpEnsembleThresholds(ensembleSample,bins, dn=0.1, up=0.9):
    """
    Example:
    thDn,thUp=calcDnUpEnsembleThresholds(ensembleSample, bins, dn=0.1, up=0.9)
    """
    cdfs =list()
    nNeurons= len(ensembleSample)
    threshDn=sc.zeros(nNeurons)
    threshUp=sc.zeros(nNeurons)
    for nn in sc.arange(nNeurons):
        mySample= ifrs[nn]
        myBins= sc.sort(sc.unique(ensembleSample))
        cdf = calcCDF(sample=mySample)
        threshDn[nn]= intCDF(dn, cdf, myBins)
        threshUp[nn]= intCDF(up, cdf, myBins)
        cdfs.append(cdf)
    return threshDn, threshUp



