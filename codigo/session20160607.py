from transientPreAnalysis import * 

# Extraction of one data set
dirName="./hipophysis/"
# Read the data directory to get a list of pxp files
pxpFiles=getFileList(dataDir=dirName, prefix="Fmin", suffix=".pxp",includeDataDir=0)
# Read the data from one file
fName = pxpFiles[3]
fMinData,sampTimes=extractDebleachedFminData(dirName+fName)
# Calculate all the areas under the curves before normalization
areas=fMinData.cumsum(1)[:,-1]
# Look for cells that may have faulty data 
goodCells=sc.where(sc.isnan(areas)==False)[0]
# Keep the good data
fMinData=fMinData[goodCells]
# Count cells and frames
nCells,nFrames= fMinData.shape
# Smooth the signals
smoothCa= sc.zeros(fMinData.shape)
for n in sc.arange(nCells):
	dd=fMinData[n]; 
	smoothCa[n]=lowHighCutFreqSmoothing(dd, lowCut=0, hiCut=40)

# Plot all the data for check
nPieces=9; q= nCells/nPieces
fig= gr.figure(figsize=(19, 11));
fig.canvas.set_window_title("%s"%(fName))
gr.ioff(); r=nPieces; c=2; 
ax1=list(); ax2=list()
for m in sc.arange(nPieces):
	ax1.append(fig.add_subplot(r,c,m*c +1))
	ax2.append(fig.add_subplot(r,c,m*c +2))
for m in sc.arange(nPieces):
	a= sc.int32(m * q)
	b= sc.int32(sc.minimum(((m+1) * q), nCells))
	dd = fMinData[a:b]
	ff = smoothCa[a:b]
	for n in sc.arange(len(dd)):
		ax1[m].plot(sampTimes, dd[n])
		ax2[m].plot(sampTimes, ff[n])
	ax1[m].set_ylim(fMinData.min()*0.97, fMinData.max()*1.01)
	ax2[m].set_ylim(fMinData.min()*0.97, fMinData.max()*1.01)
fig.subplots_adjust(left=0.05, right=0.98, top=0.98, bottom=0.05)
gr.ion(); gr.draw()


# -------------------------------
# Now let us build a function that extracts all the data contained in a .pxp file within a given directory, 
# smooths the data and normalizes all amplitudes in different ways
# -------------------------------
def getSmoothAllPXPData(dirName,fName,norm01=1,shift0=1,shift0Max=1):
	# Read the data from one file
	print("Extracting data from file %s"%(dirName+fName))
	fMinData,sampTimes=extractDebleachedFminData(dirName+fName)
	# Calculate all the areas under the curves before normalization
	areas=fMinData.cumsum(1)[:,-1]
	# Look for cells that may have faulty data 
	goodCells=sc.where(sc.isnan(areas)==False)[0]
	areas= areas[goodCells]
	# Keep the good data
	transCa=fMinData[goodCells]
	# Count cells and frames
	nCells,nFrames= transCa.shape
	# Smooth the signals
	smoothCa= sc.zeros(fMinData.shape)
	smoothCa01= sc.zeros(fMinData.shape)
	smoothCa0= sc.zeros(fMinData.shape)
	smoothCa0Max= sc.zeros(fMinData.shape)
	for n in sc.arange(nCells):
		smoothCa[n]=lowHighCutFreqSmoothing(transCa[n], lowCut=0, hiCut=40)
		if norm01:
			smoothCa01[n]= (smoothCa[n] - smoothCa[n].min())/(smoothCa[n].max() - smoothCa[n].min())
		if shift0:
			smoothCa0[n]= (smoothCa[n] - smoothCa[n].min())
	mCa=smoothCa.min()
	MCa=smoothCa.max()
	for n in sc.arange(nCells):
		if shift0Max:
			smoothCa0Max[n]= (smoothCa[n] - mCa)/(MCa - mCa)
	data={"smoothCa": smoothCa, "transCa": transCa, 
		"fName":fName, "nCells": nCells,
		"sampTimes":sampTimes, "areas":areas}
	if norm01:
		data["smoothCa01"]=smoothCa01
	if shift0:
		data["smoothCa0"]=smoothCa0
	if shift0Max:
		data["smoothCa0Max"]=smoothCa0Max
	return data

pericitosBasal= getSmoothAllPXPData(dirName=dirName, fName="FminBasalPericitos.pxp")
endocrinasBasal= getSmoothAllPXPData(dirName=dirName, fName="FminBasalEndocrinas.pxp")

# Plot all the data for check
recs=pericitosBasal
recs=endocrinasBasal
nPieces=9; q= nCells/nPieces
fig= gr.figure(figsize=(19, 11));
fig.canvas.set_window_title("%s"%(recs["fName"]))
gr.ioff(); r=nPieces; c=5; 
ax1=list(); ax2=list(); ax3=list(); ax4=list(); ax5=list()
for m in sc.arange(nPieces):
	ax1.append(fig.add_subplot(r,c,m*c +1))
	ax2.append(fig.add_subplot(r,c,m*c +2))
	ax3.append(fig.add_subplot(r,c,m*c +3))
	ax4.append(fig.add_subplot(r,c,m*c +4))
	ax5.append(fig.add_subplot(r,c,m*c +5))
for m in sc.arange(nPieces):
	a= sc.int32(m * q)
	b= sc.int32(sc.minimum(((m+1) * q), nCells))
	dd = recs["transCa"][a:b]
	ff = recs["smoothCa"][a:b]
	gg = recs["smoothCa01"][a:b]
	hh = recs["smoothCa0"][a:b]
	ii = recs["smoothCa0Max"][a:b]
	#print(len(recs["sampTimes"]),len(dd[n]))
	for n in sc.arange(len(dd)):
		ax1[m].plot(recs["sampTimes"], dd[n])
		ax2[m].plot(recs["sampTimes"], ff[n])
		ax3[m].plot(recs["sampTimes"], gg[n])
		ax4[m].plot(recs["sampTimes"], hh[n])
		ax5[m].plot(recs["sampTimes"], ii[n])

	ax1[m].set_ylim(recs["transCa"].min()*0.97, recs["transCa"].max()*1.01)
	ax2[m].set_ylim(recs["transCa"].min()*0.97, recs["transCa"].max()*1.01)
	ax3[m].set_ylim(-0.01, 1.01)
	ax4[m].set_ylim(-0.01, 0.1)
	ax5[m].set_ylim(-0.01, 1.01)
fig.subplots_adjust(left=0.05, right=0.98, top=0.98, bottom=0.05)
gr.ion(); gr.draw()

