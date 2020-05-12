import csv
import os
import string
import scipy as sc
import matplotlib.pyplot as gr
gr.ion()
import scipy.io as sio
import scipy.signal as signal

dataDir = "microcircuitsNetworks/"
dataFile = 'duhne20140816-1.mat'
data=sio.loadmat(dataDir+dataFile);
spikes= data['Spikes'].transpose()
nCells = len(spikes)



def delayedCounts(train1, train2, nBins, normalize=1):
    n1=len(train1)
    n2=len(train2)
    nTBins= 2*nBins + 1
    nR=nBins; nL=nBins
    counts = sc.zeros(nTBins)
    iOnes1=sc.where(train1)[0]
    iOnes2=sc.where(train2)[0]
    nSpikes1=len(iOnes1)
    nSpikes2=len(iOnes2)
    
    for k in sc.arange(nSpikes1):
        dummy = sc.zeros(nTBins)
        myInd= iOnes1[k]
        #print(myInd)
        b= iOnes1[k]+nR
        if myInd-nL>=0:
            a= iOnes1[k]-nL
            aa=0
        else:
            a= 0;
            aa= nL-myInd
        if nR+myInd < n2:
            b= 1+nR+iOnes1[k]
            bb= nTBins 
        else:
            b= n2
            bb= nTBins - (myInd+nR-b)-1
        #print(a,b,aa,bb)
        #print(len(spikes[n][a:b]))
        #print(len(dummy[aa:bb]))
        dummy[aa:bb]= train2[a:b]
        counts= counts + dummy;
    if sc.product(train1==train2)>0:
        counts[nL]=0.0
    if ((normalize) & (nSpikes2>0)):
        counts= counts/ sc.float32(nSpikes2);
    return counts


xCorrs= list()
nBins=100
for m in sc.arange(1):
    for n in sc.arange(nCells):
        xCorrs.append(delayedCounts(spikes[m], spikes[n],nBins=nBins))

        
bins=sc.arange(-nBins,nBins+1)
start=0; cols=10; rows=1
fig= gr.figure()
ax=list(); 
gr.ioff()
for m in sc.arange(start,start+rows):
    for n in sc.arange(start,start+cols):
        k= m*cols + n
        h = 100 * xCorrs[k] 
        ax.append(fig.add_subplot(rows,cols, k+1))
        ax[k].bar(bins-0.5, h)
        ax[k].set_ylim(-0.1, h.max())
gr.ion(); gr.draw()


