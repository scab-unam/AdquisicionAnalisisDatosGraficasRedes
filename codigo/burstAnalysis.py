from spikeTrainAnalysis import *
gr.ion()
# Data from SN Reticulata
dDir="./datosTrenes/trenesBursts/"
prefs= ["Cnt", "Ant", "Les"]
# -------------------------------------
# Thresholds from ifr
# -------------------------------------
thrISI=dict(); thrIFR=dict(); spBr=dict()
# Thresholds Ale
thrISI["Cnt"]=sc.array([0.1, 0.07079, 0.07943, 0.07079, 0.3548, 0.7079, 0.2239, 0.1122, 0.1259, 0.07079])
thrISI["Ant"]=sc.array([1.0, 0.3548, 0.4467, 0.631, 0.7079, 0.2239, 0.7079, 0.631, 0.2239, 0.3981, 0.1995])
thrISI["Les"]=sc.array([0.2512, 0.3548, 0.5623, 0.5623, 0.4467, 0.1585, 0.2818, 0.7079])
# Thresholds Ale, IFR
thrIFR["Cnt"]=1/thrISI["Cnt"]
thrIFR["Ant"]=1/thrISI["Ant"]
thrIFR["Les"]=1/thrISI["Les"]
# Proposed changes
# Mean for the number of spikes in burst
spBr["Cnt"]=sc.array([3, 2, 2, 2, 2, 3, 2, 4, 3, 2])-1
spBr["Ant"]=sc.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])-1
spBr["Les"]=sc.array([2, 2, 3, 3, 2, 2, 2, 2])-1
#thrIFR["Cnt"][0:7]=19.0,10.0, 9.0, 9.0, 9.0, 3.0, 5.0
#thrIFR["Cnt"][-1]=2.0
#thrISI["Cnt"]=1/thrIFR["Cnt"]

# ------------------------------
# Data extraction and preliminary analyisis 
# ------------------------------
trains= dict()
for pref in prefs:
	trains[pref]= xtractSpikeTrains(dataDir=dDir,prefix=pref+"St",suff=".mat")
	#trains["Cnft"]= xtractTrains(dataDir=dDir,prefix=pref+"St",suff=".mat")
	nTrains=len(trains[pref])
	for cell in sc.arange(nTrains):
		st= trains[pref][cell]
		st.thrIFR = thrIFR[pref][cell]
		st.thrISI = thrISI[pref][cell]
		st.spBr= spBr[pref][cell]
		st.ifr=st.calcIFR(order=st.spBr)
		st.burstInds=st.burstIFR(st.thrIFR)
		#st.burstPercent= 100.0*len(sc.where(st.burstInds>st.thrIFR)[0])/st.nSpikes
		st.burstPercent= 100.0*len(st.burstInds)/st.nSpikes

# ------------------------------
# Regroup bursting measurements for analysis
# ------------------------------
percentiles=sc.arange(0,100.0000001)
bDurs=dict();bDursQuants=dict()
bIntervals=dict(); bIntervalsQuants=dict()
bDuty=dict(); bDutyQuants=dict()
bQuiesc=dict(); bQuiescQuants=dict()
fNs=dict()
for pref in prefs:
	nTrains=len(trains[pref])
	bDurs[pref]=list(); bDursQuants[pref]=sc.zeros((nTrains,len(percentiles)))
	bIntervals[pref]=list(); bIntervalsQuants[pref]=sc.zeros((nTrains,len(percentiles)))
	bQuiesc[pref]=list(); bQuiescQuants[pref]=sc.zeros((nTrains,len(percentiles)))
	bDuty[pref]=list(); bDutyQuants[pref]=sc.zeros((nTrains,len(percentiles)))
	fNs[pref]=list()
	for cell in sc.arange(nTrains):
		#print("Starting with cell %d"%cell)
		st=trains[pref][cell]
		fNs[pref].append(st.fName[:-4])
		if st.nBursts>0:
			bDurs[pref].append(st.bDurations)
			bIntervals[pref].append(st.bIntervals)
			bQuiesc[pref].append(st.bQuiesc)
			bDuty[pref].append(st.bDurations)
			if st.nBursts>1:
				bDursQuants[pref][cell]=stats.scoreatpercentile(a=st.bDurations,per=percentiles)
				bIntervalsQuants[pref][cell]=stats.scoreatpercentile(a=st.bIntervals,per=percentiles)
				bDutyQuants[pref][cell]=stats.scoreatpercentile(a=st.bDuty,per=percentiles)
				bQuiescQuants[pref][cell]=stats.scoreatpercentile(a=st.bQuiesc,per=percentiles)
			else:
				bDursQuants[pref][cell]=st.bDurations.sort()
				bIntervalsQuants[pref][cell]=st.bIntervals.sort()
				bDutyQuants[pref][cell]=st.bDuty.sort()
				bQuiescQuants[pref][cell]=st.bQuiesc.sort()
		else:
			bDurs[pref].append([0])
			bIntervals[pref].append([0])
			bQuiesc[pref].append([0])
			bDuty[pref].append([0])
		#print("Done with cell %d"%cell)

# ------------------------------
# Graphical comparison between burst measurements WITHIN GROUPS
# ------------------------------
compareBurstMeasures=1
if compareBurstMeasures:
	for m in sc.arange(len(prefs)):
		f= gr.figure(figsize=(11,11)); gr.ioff()
		colors= ["blue","orange","purple"]
		rows=4; cols=1
		axBInter=f.add_subplot(rows,cols,1)
		axBDur=f.add_subplot(rows,cols,2)
		axBDuty=f.add_subplot(rows,cols,3)
		axBQuiesc=f.add_subplot(rows,cols,4)
		pref=prefs[m]
		nTrains=len(bDurs[pref])
		for n in sc.arange(nTrains):
			x = bDurs[pref][n]
			fMean,fMeanInv= calcCDF(x)
			axBDur.plot(x, fMean(x),'.',alpha=0.3,lw=1)
		for n in sc.arange(nTrains):
			x = bIntervals[pref][n]
			fMean,fMeanInv= calcCDF(x)
			axBInter.plot(x, fMean(x),'.',alpha=0.3,lw=1)
		for n in sc.arange(nTrains):
			xDuty = bDuty[pref][n]
			fMean,fMeanInv= calcCDF(xDuty)
			axBDuty.plot(xDuty, fMean(xDuty),'.',alpha=0.3,lw=1)
		for n in sc.arange(nTrains):
			x = bQuiesc[pref][n]
			fMean,fMeanInv= calcCDF(x)
			axBQuiesc.plot(x, fMean(x),'.',alpha=0.3,lw=1)
		bi = bIntervalsQuants[pref].mean(0)
		bd = bDursQuants[pref].mean(0)
		bduty = bDutyQuants[pref].mean(0)
		bq = bQuiescQuants[pref].mean(0)
		fMean,fMeanInv= calcCDF(bi)
		axBInter.plot(bi, fMean(bi),colors[m],alpha=1,lw=3,label=pref)
		#axBInter.plot(bi, fMean(bi),colors[m]+'o',alpha=1,ms=5,label=pref)
		fMean,fMeanInv= calcCDF(bd)
		axBDur.plot(bd, fMean(bd),colors[m],alpha=1,lw=3,label=pref)
		#axBDur.plot(bd, fMean(bd),colors[m]+'o',alpha=1,ms=5,label=pref)
		fMean,fMeanInv= calcCDF(bduty)
		axBDuty.plot(100*bduty, fMean(bduty),colors[m],alpha=1,lw=3,label=pref)
		#axBDuty.plot(100*bduty, fMean(bduty),colors[m]+'o',alpha=1,ms=5,label=pref)
		fMean,fMeanInv= calcCDF(bq)
		axBQuiesc.plot(bq, fMean(bq),colors[m],alpha=1,lw=3,label=pref)
		#axBQuiesc.plot(bq, fMean(bq),colors[m]+'o',alpha=1,ms=5,label=pref)
		axBInter.set_xlabel('seconds'); axBInter.legend(); axBInter.set_title("Burst interval")
		axBDur.set_xlabel('seconds'); axBDur.legend(); axBDur.set_title("Burst duration")
		axBDuty.set_xlabel('percentage'); axBDuty.legend(); axBDuty.set_title("Burst duty percentage")
		axBQuiesc.set_xlabel('seconds'); axBQuiesc.legend(); axBQuiesc.set_title("Burst quiescence")
		f.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95,wspace=0.2, hspace=0.3)
		gr.ion(); gr.draw()
	
	#f.savefig("./ifrThresholding/burstMeasuresAll.png")


# ------------------------------
# Graphical comparison between burst measurements
# ------------------------------
compareBurstMeasures=1
if compareBurstMeasures:
	f= gr.figure(figsize=(11,13)); gr.ioff()
	colors= ["blue","orange","purple"]
	rows=4; cols=1
	axBInter=f.add_subplot(rows,cols,1)
	axBDur=f.add_subplot(rows,cols,2)
	axBDuty=f.add_subplot(rows,cols,3)
	axBQuiesc=f.add_subplot(rows,cols,4)
	for m in sc.arange(len(prefs)):
		pref=prefs[m]
		nTrains=len(bDurs[pref])
		#for n in sc.arange(nTrains):
		#	x = bDurs[pref][n]
		#	fMean,fMeanInv= calcCDF(x)
		#	gr.plot(x, fMean(x),colors[m]+'.',alpha=0.59,lw=1)
		bi = bIntervalsQuants[pref].mean(0)
		bd = bDursQuants[pref].mean(0)
		bduty = bDutyQuants[pref].mean(0)
		bq = bQuiescQuants[pref].mean(0)
		fMean,fMeanInv= calcCDF(bi)
		axBInter.plot(bi, fMean(bi),colors[m],alpha=1,lw=3,label=pref)
		#axBInter.plot(bi, fMean(bi),colors[m]+'o',alpha=1,ms=5,label=pref)
		fMean,fMeanInv= calcCDF(bd)
		axBDur.plot(bd, fMean(bd),colors[m],alpha=1,lw=3,label=pref)
		#axBDur.plot(bd, fMean(bd),colors[m]+'o',alpha=1,ms=5,label=pref)
		fMean,fMeanInv= calcCDF(bduty)
		axBDuty.plot(100*bduty, fMean(bduty),colors[m],alpha=1,lw=3,label=pref)
		#axBDuty.plot(100*bduty, fMean(bduty),colors[m]+'o',alpha=1,ms=5,label=pref)
		fMean,fMeanInv= calcCDF(bq)
		axBQuiesc.plot(bq, fMean(bq),colors[m],alpha=1,lw=3,label=pref)
		#axBQuiesc.plot(bq, fMean(bq),colors[m]+'o',alpha=1,ms=5,label=pref)
	axBInter.set_xlabel('seconds'); axBInter.legend(); axBInter.set_title("Burst interval")
	axBDur.set_xlabel('seconds'); axBDur.legend(); axBDur.set_title("Burst duration")
	axBDuty.set_xlabel('percentage'); axBDuty.legend(); axBDuty.set_title("Burst duty percentage")
	axBQuiesc.set_xlabel('seconds'); axBQuiesc.legend(); axBQuiesc.set_title("Burst quiescence")
	f.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95,wspace=0.2, hspace=0.3)
	gr.ion(); gr.draw()
	f.savefig("./figsIFRThresholding/burstMeasuresAll.png")





percentiles=[0,.25,.5,.75,1.0]
if 0:
	strBDur="Data &  & %g & %g & %g & %g \\\\"%(st.fName[:-4], q0,q25,q50,q75,q100)
	print(strBDur+"\n")
	for pref in prefs:
		nTrains=len(trains[pref])
		for cell in sc.arange(nTrains):
			st=trains[pref][cell]
			q0,q25,q50,q75,q100= stats.scoreatpercentile(a=st.bDurations,per=percentiles)
			strBDur="%s %g & %g & %g & %g & %g \\\\"%(st.fName[:-4], q0,q25,q50,q75,q100)
			print(strBDur+"\n")

showIFRThresholding=1
ifrMax=dict()
ifrMax["Cnt"]=[40,40,50,50,10,10,21,21,40,40]
ifrMax["Les"]=[84,24,13,13,15,45,25,13]
ifrMax["Ant"]=[20,20,20,20,20,20,20,20,20,20,20]
if showIFRThresholding:
	for pref in prefs:
		nTrains=len(trains[pref])
		f=gr.figure(figsize=(19,15)); gr.ioff()
		f.canvas.set_window_title(pref)
		cols=2; rows = sc.ceil(nTrains/cols)
		ax=list()
		for cell in sc.arange(nTrains):
			#pref="Cnt"; cell=6
			st=trains[pref][cell]
			ax.append(f.add_subplot(rows,cols,cell+1))
			str0="%s %g percent burst"%(st.fName[:-4],st.burstPercent)
			str1=r"ISI thresh order %d, %g secs"%(st.spBr, st.thrISI)
			gr.plot(st.spikeTimes,st.thrIFR*sc.ones(st.nSpikes),"b",label=str0+"\n"+str1); 
			ifr=st.calcIFR(1)
			gr.plot(st.spikeTimes,st.calcIFR(1),"wo", label=r"IFR order %d"%1); 
			if st.spBr>1:
				gr.plot(st.spikeTimes,st.calcIFR(st.spBr),"b.",label=r"IFR order %d"%st.spBr); 
			gr.xlim(0,st.spikeTimes.max())
			gr.ylim(0,1.1*ifrMax[pref][cell])
			gr.legend(ncol=1,fontsize=10); 
		gr.ion(); gr.draw()
		f.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.98,wspace=0.1, hspace=0.175)
		f.savefig("./figsIFRThresholding/"+pref+"IFRThresholding.png")
		gr.close("all")
#
showBurstDurations=1
if showBurstDurations:
	f=gr.figure(figsize=(17,13)); gr.ioff()
	f.canvas.set_window_title("Burst durations")
	f.suptitle("Burst durations")
	cols=1; rows=len(prefs)
	ax=list()
	p=1
	for pref in prefs:
		ax.append(f.add_subplot(rows,cols,p))
		nTrains=len(trains[pref])
		bDs=list(); fNs=list()
		for cell in sc.arange(nTrains):
			st=trains[pref][cell]
			fNs.append(st.fName[:-4])
			if len(st.bDurations)>0:
				bDs.append(st.bDurations)
			else:
				bDs.append([0])
		#gr.figure(figsize=(13,13))
		#gr.boxplot(st.bDurations, positions=[cell+1,]) #, labels=[st.fName[:-4]])
		ax[p-1].boxplot(bDs, 0, 'ko', labels=fNs) #, labels=[st.fName[:-4]])
		ax[p-1].set_ylim(0,12)
		p=p+1
	f.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95, wspace=0.1, hspace=0.2)
	gr.ion(); gr.draw()
	f.savefig("./figsIFRThresholding/burstDurationBoxes.png")

# Quiescence intervals
showQuiesInterval=1
if showQuiesInterval:
	f=gr.figure(figsize=(17,13)); gr.ioff()
	f.canvas.set_window_title("Burst quiescence intervals")
	f.suptitle("Burst quiescence")
	cols=1; rows=len(prefs)
	ax=list()
	p=1
	for pref in prefs:
		ax.append(f.add_subplot(rows,cols,p))
		nTrains=len(trains[pref])
		bQs=list(); fNs=list()
		for cell in sc.arange(nTrains):
			st=trains[pref][cell]
			fNs.append(st.fName[:-4])
			if len(st.bQuiesc)>0:
				bQs.append(st.bQuiesc)
			else:
				bQs.append([0])
		#gr.figure(figsize=(13,13))
		#gr.boxplot(st.bDurations, positions=[cell+1,]) #, labels=[st.fName[:-4]])
		ax[p-1].boxplot(bQs, 0, 'ko', labels=fNs) #, labels=[st.fName[:-4]])
		#ax[p-1].set_ylim(0,12)
		p=p+1
	f.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95, wspace=0.1, hspace=0.2)
	gr.ion(); gr.draw()
	f.savefig("./figsIFRThresholding/burstQuiescenceBoxes.png")

# Duty cycles
showDuty=1
if showDuty:
	f=gr.figure(figsize=(17,13)); gr.ioff()
	f.canvas.set_window_title("Burst duty percentage")
	f.suptitle("Burst duty percentage (burst duration/burst interval)")
	cols=1; rows=len(prefs)
	ax=list()
	p=1
	for pref in prefs:
		ax.append(f.add_subplot(rows,cols,p))
		nTrains=len(trains[pref])
		bDuties=list(); fNs=list()
		for cell in sc.arange(nTrains):
			st=trains[pref][cell]
			fNs.append(st.fName[:-4])
			if len(st.bDuty)>0:
				bDuties.append(100*st.bDuty)
			else:
				bDuties.append([0])
		#gr.figure(figsize=(13,13))
		#gr.boxplot(st.bDurations, positions=[cell+1,]) #, labels=[st.fName[:-4]])
		ax[p-1].boxplot(bDuties, 0, 'ko', labels=fNs) #, labels=[st.fName[:-4]])
		#ax[p-1].set_ylim(0,12)
		p=p+1
	f.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95, wspace=0.1, hspace=0.2)
	gr.ion(); gr.draw()
	f.savefig("./ifrThresholding/burstDutyBoxes.png")

# Interburst intervals
showBIntervals=1
if showBIntervals:
	f=gr.figure(figsize=(17,13)); gr.ioff()
	f.canvas.set_window_title("Interburst intervals")
	f.suptitle("Interburst intervals")
	cols=1; rows=len(prefs)
	ax=list()
	p=1
	for pref in prefs:
		ax.append(f.add_subplot(rows,cols,p))
		nTrains=len(trains[pref])
		bInt=list(); fNs=list()
		for cell in sc.arange(nTrains):
			st=trains[pref][cell]
			fNs.append(st.fName[:-4])
			if len(st.bIntervals)>0:
				bInt.append(st.bIntervals)
			else:
				bInt.append([0])
		#gr.figure(figsize=(13,13))
		#gr.boxplot(st.bDurations, positions=[cell+1,]) #, labels=[st.fName[:-4]])
		ax[p-1].boxplot(bInt, 0, 'ko', labels=fNs) #, labels=[st.fName[:-4]])
		ax[p-1].set_ylim(0,100)
		p=p+1
	f.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95, wspace=0.1, hspace=0.2)
	gr.ion(); gr.draw()
	f.savefig("./figsIFRThresholding/burstIntervalBoxes.png")
#


