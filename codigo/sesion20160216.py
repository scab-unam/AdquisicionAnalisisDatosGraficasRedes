from transientPreAnalysis import * 
#import transientPreAnalysis PA

# Extraction of one data set
dirName="microcircuitsNetworks/"
fName="cort76dp1c.pxp"
#rawData=extractPXPData(dirName+fName)
fMinData=extractDebleachedFminData(dirName+fName)
# Calculate areas under the curve
areas=fMinData.cumsum(1)[:,-1]
# Look for cells that may have faulty data 
goodCells=sc.where(sc.isnan(areas)==False)[0]
# Keep the good data
fMinData=fMinData[goodCells]
nCells,nFrames= fMinData.shape
frameMS=150.00
# areas = (areas - areas.min() )/ ()
# Save the data to a pickle file
savePickle=1
if savePickle:
    outFile = open(dirName+fName[:-4]+".pkl", 'wb')
    pickle.dump(fMinData, outFile)
    outFile.close()
# Normalize the calcium traces by the peak 
if 0:
    normAreas= sc.zeros(nCells)
    for m in range(len(areas)):
        normAreas[m]=areas[m]/fMinData[m].max()
        
normAreas= (areas - areas.min())/(areas.max() - areas.min())
# Get the cumulative distribution of the normalized areas
cdfArea,cdfInvArea= calcCDF(normAreas)
# Set significance values to calculate thresholds for the areas
signif=sc.arange(0.90,0.96,0.05)
# Uncomment this line if more than one threshold is used
threshArea= calcThresholdsFromCDF(cdfInvArea, quants=signif)
#threshArea= calcThresholdFromCDF(cdfInvArea, signif[0])
areaArray= sc.arange(0,1,0.01)
sortedAreas= sc.sort(normAreas)
nThresh=len(threshArea)
cellsSignif=list()
for m in sc.arange(nThresh):
    cellsSignif.append(sc.where(normAreas>threshArea[m])[0])
    print(cellsSignif[m])

# Histogram of areas under the curve
nBars=100
figAreas=gr.figure()
rows=1; cols=2
ax1= figAreas.add_subplot(rows,cols,1)
ax2= figAreas.add_subplot(rows,cols,2)
ax1.hist(areas, nBars, label=r"$A=\int_0^T f(t) dt$")
ax2.hist(normAreas, nBars,label=r"$A=\frac{1}{\max{f(t)}} \int_0^T f(t) dt$")
ax1.set_xlabel(r"Area ($A$, arbitrary units)")
ax2.set_xlabel(r"Area ($A$, arbitrary units)")
ax1.legend()
gr.show()
# Normalize the Ca-transients
# Envelope
maxMinF=rawData.max(0)- rawData.min(0)
print(maxMinF[0],maxMinF[-1])
avgF=rawData.mean(0)
print(avgF[0],avgF[-1])
# Normalization by fMin
#fMinData= fMinNorm(rawData)
#fMaxMinData= fMaxMinNorm(rawData)
d1= fMaxMinData[100]
d2= fMaxMinData[130]
gr.figure(); gr.ioff()
gr.plot(d1)
gr.plot(d2)
gr.ion(); gr.draw()

# Calculation for the areas under the curves of the Ca-Transients
# Keep only the last element of the cumulative sum over time to get the total activity and normalize


bins= sc.arange(0, 1, 1/20.0)
counts, binEdges= sc.histogram(normAreas,bins)


# Graficas exploratoria de todas las trazas     
#if 1:
figRaw=gr.figure(figsize=(13,10))
gr.ioff(); axRaw=list(); axFMin=list(); axFMaxMin=list()
rows=3; cols=2
k= sc.int32(0.05*nCells/rows)
for m in sc.arange(rows):
    a= m*k; b= sc.minimum((m+1)*k, nCells)
    #axRaw.append(figRaw.add_subplot(rows,cols,m*cols +1))
    axFMin.append(figRaw.add_subplot(rows,cols,m*cols +1))
    axFMaxMin.append(figRaw.add_subplot(rows,cols,m*cols +2))
    print(a,b)
    for n in sc.arange(a,b):
        #axRaw[m].plot(rawData[n],alpha=0.3)
        axFMin[m].plot(fMinData[n],alpha=0.3)
        axFMaxMin[m].plot(fMinData[n]/fMinData[n].max(),alpha=0.3)
        #axRaw[m].set_ylabel("F")
        axFMin[m].set_ylabel(r"$\frac{F-F_{m}}{F_{m}}$")
        axFMaxMin[m].set_ylabel(r"$\frac{F-F_{m}}{F_{M}F_{m}}$")
    if m>rows-2:
        #axRaw[m].set_xlabel("Frames")
        axFMin[m].set_xlabel("Frames")

gr.ion(); gr.draw()


for n in sc.arange(nCells):
    gr.plot(rawData[n],aspect='auto')
gr.ion(); gr.draw()



# Graphs of the total area distributions
rows=3; cols=1; barW=0.85*(binEdges[1]-binEdges[0])
f=gr.figure()
gr.ioff()
ax=list()
for n in sc.arange(rows*cols):
    ax.append(f.add_subplot(rows,cols,n+1))

ax[0].plot(normAreas,sc.arange(nCells),'o')
ax[1].bar(left=binEdges[:-1],width=barW,height=counts)
ax[2].plot(areaArray, cdfArea(areaArray),'k', label='CDF')
ax[2].plot(sortedAreas, cdfArea(sortedAreas), 'b.', label='CDF')

if nThresh>1:
    for mm in sc.arange(nThresh):
        ax[2].plot([threshArea[mm], threshArea[mm]], [0,signif[mm]], 'k',  alpha=0.35, label='CDF')
        ax[2].plot([0, threshArea[mm]], [signif[mm],signif[mm]], 'k',  alpha=0.35, label='CDF')
        ax[1].plot([threshArea[mm], threshArea[mm]], [0,counts.max()], 'k',  alpha=0.35, label='CDF')
        ax[0].plot([threshArea[mm], threshArea[mm]], [0,nCells+1], 'k',  alpha=0.35, label='CDF')
else:
    ax[2].plot([threshArea[0], threshArea[0]], [0,signif[0]], 'k', alpha=0.35, label='CDF')

ax[0].set_ylabel("Cells")
ax[2].set_xlabel("Normalized area")
ax[1].set_ylabel("Counts")
ax[2].set_ylabel(r"$P$(area $\leq a$")
ax[0].set_xlim(0,1)
ax[0].set_ylim(0,nCells+1)
ax[1].set_xlim(0,1)
ax[2].set_xlim(0,1)
gr.ion(); gr.draw()

if 0:
    gr.figure()
    gr.ioff(); 
    gr.imshow(rawData,aspect='auto')
    gr.ion(); gr.draw()



if 0:
    # Smoothing attempt
    # Distribution of the time-dependent rates of change of the normalized fluorescence data
    dFdt= sc.zeros((nCells,nFrames),'float64')
    smoothdFdt= sc.zeros((nCells,nFrames),'float64')
    s3F= sc.zeros((nCells,nFrames),'float64')
    for m in sc.arange(nCells):
        dFdt[m][1:]= (fMinData[m][1:]- fMinData[m][:-1])
        s3F[m][1:-1]= (F[2:]+F[1:-1]+F[0:-2])/3.0
        smoothdFdt[m][1:]= (s3F[m][1:]- s3F[m][:-1])

    dFdt= dFdt/frameMS
    smoothdFdt= smoothdFdt/frameMS

    q=0.50
    threshDFs=list()
    cdfDFs=list()
    cdfInvDFs=list()
    indHighDFs=list()
    for m in sc.arange(nCells):
        #cdfdFdt,cdfInvdFdt= calcCDF(smoothdFdt[m])
        cdfdFdt,cdfInvdFdt= calcCDF(dFdt[m])
        cdfDFs.append(cdfdFdt)
        cdfInvDFs.append(cdfInvdFdt)
        threshdFdt= calcThresholdFromCDF(cdfInvdFdt, quant=q)
        threshDFs.append(threshdFdt)
        indHighDFs.append(sc.where(dFdt[m]>threshdFdt)[0])
        print(m, threshdFdt)


    # Smoothing data to obtain events 
    m=200
    dF= dFdt[m]
    F= fMinData[m]
    hI= sc.hstack([0,indHighDFs[m],len(F)-1])
    secs= frameMS * sc.arange(0,nFrames)/1000.0
    s = sc.interpolate.interp1d(secs[hI], F[hI], kind='cubic')
    interpF = s(secs)
    #
    fig=gr.figure()
    gr.ioff()
    rows=1; cols=2
    ax1=fig.add_subplot(rows,cols,1)
    ax2=fig.add_subplot(rows,cols,2)
    ax1.plot(secs,F,'k',alpha=0.35)
    ax1.plot(secs[hI],F[hI],'wo')
    ax1.plot(secs,interpF,'r')
    ax1.plot(secs[hI],interpF[hI],'ro')
    ax2.plot(secs,dF,'k')
    ax2.plot(secs[hI],dF[hI],'wo')
    gr.ion(); gr.draw()


    # Plot all the distributions of the derivatives to see where the thresholds for spiking should be
    gr.figure(); 
    gr.ioff()
    rangeDF=0.0002*sc.arange(-1,1,0.01)
    for m in sc.arange(nCells):
        myFunc=cdfDFs[m]
        gr.plot(rangeDF, myFunc(rangeDF))

    gr.ion(); gr.draw(); 
