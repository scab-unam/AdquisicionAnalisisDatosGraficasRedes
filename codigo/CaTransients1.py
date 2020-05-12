import csv
import os
import string
import scipy as sc
import matplotlib.pyplot as gr
gr.ion()
import scipy.io as sio

dataDir = 'microcircuitsNetworks/'
dataFile = 'duhne20140816-1.mat'
data=sio.loadmat(dataDir+dataFile);
spikes= data['Spikes'].transpose()

myarray = sc.fromfile('BinaryData.dat',dtype=float)

def templatesFromCSVs(dataDir='microcircuitsNetworks/CalciumOrange/', dataIndicator='T_', graphData=1):
    # Read file names in the directory
    timeFiles = [ f for f in os.listdir(dataDir) if os.path.isfile(os.path.join(dataDir,f)) & string.find(f,dataIndicator)<1 ]
    nFiles = len(timeFiles)
    Data=dict()
    data['names']=list()
    m=0
    data['nTotFrames']=0
    if graphData: 
        fig=gr.figure()
        gr.ioff()
        rows=nFiles; cols=1
    for fileName in timeFiles:
        delimiter=','
        thisName =fileName[len(dataIndicator):-4]
        data['names'].append(thisName)
        data[thisName]=dict()
        with open(dataDir+fileName,'rb') as f:
            reader=csv.reader(f, delimiter=delimiter)
            data[thisName]['samplingRate']=reader.next()
            data[thisName]['acquisitionDate']=reader.next()
            data[thisName]['row3']=reader.next()
            data[thisName]['row4']=reader.next()

            dataList=list()
            for row in reader:
                aa=sc.float32(row)
                dataList.append(aa)
            data[thisName]['dF/F']=sc.array(dataList).transpose()
            data[thisName]['nFrames']= data[thisName]['dF/F'].shape[1]
            data[thisName]['nCells']= data[thisName]['dF/F'].shape[0]
            data['nTotFrames']=data['nTotFrames']+data[thisName]['nFrames']

        if graphData:
            gr.imshow(data[thisName]['dF/F'],aspect='auto')
            gr.title(thisName)
            #ax.colorbar()
        m=m+1
        gr.ion(); gr.draw()
                
    data['nFiles']=nFiles
    data['nCells']=data[thisName]['nCells']
    return data

data = templatesFromCSVs(graphData=1)

# Quantify the distributions for the different conditions
sampRate=6.0
quantValues=[0.5,0.95]
cell=50
cdfs= list()
sampTimes=list()
bins=list()
samples=list()
thresholds=list()
allFrames=sc.zeros(data['nTotFrames'],'float64')
allSampTimes=sc.arange(0,data['nTotFrames']/sampRate, 1/sampRate)
a=0; b=0
for mm in sc.arange(data['nFiles']):
    f1 = data[data['names'][mm]]
    sampTimes.append(sc.arange(0,f1['nFrames']/sampRate, 1/sampRate) )
    sample=f1['dF/F'][cell]
    samples.append(sample)
    b = a+ f1['nFrames']
    allFrames[a:b]=f1['dF/F'][cell]
    bins.append(sc.arange(sample.min(), sample.max(),  (sample.max()-sample.min())/100.0))
    cdf =calcCDF( sample)
    cdfs.append( cdf)
    thresholds.append(calcThresholdsFromCDF(cdf, quantValues))
    a = b

fig=gr.figure(); gr.ioff()
rows=data['nFiles']; cols =2
ax=list()
for n in sc.arange(rows*cols):
    ax.append(fig.add_subplot(rows,cols,n +1))

for n in sc.arange(rows):
    print(n)
    ax[2*n].plot(sampTimes[n], samples[n])
    ax[2*n + 1].plot( cdfs[n](bins[n]), bins[n])
gr.ion(); gr.draw(); 


    

