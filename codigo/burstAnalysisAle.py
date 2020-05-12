from spikeTrainAnalysis import *
# Data from SN Reticulata
dDir="./datosTrenes/trenesBursts/"
prefs= ["Cnt", "Ant", "Les"]
# -------------------------------------
# Thresholds from ifr
# -------------------------------------
thrISI=dict(); thrIFR=dict()
# Thresholds Ale
thrISI["Cnt"]=sc.array([0.1, 0.07079, 0.07943, 0.07079, 0.3548, 0.7079, 0.2239, 0.1122, 0.1259, 0.07079])
thrISI["Ant"]=sc.array([1.0, 0.3548, 0.4467, 0.631, 0.7079, 0.2239, 0.7079, 0.631, 0.2239, 0.3981, 0.1995])
thrISI["Les"]=sc.array([0.2512, 0.3548, 0.5623, 0.5623, 0.4467, 0.1585, 0.2818, 0.7079])
#min ob spikes in burst
minSpBr["Cnt"]=sc.array([3, 2, 2, 2, 2, 3, 2, 4, 3, 2])
minSpBr["Ant"]=sc.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
minSpBr["Les"]=sc.array([2, 2, 3, 3, 2, 2, 2, 2])
# Thresholds Ale, IFR
thrIFR["Cnt"]=1/thrISI["Cnt"]
thrIFR["Ant"]=1/thrISI["Ant"]
thrIFR["Les"]=1/thrISI["Les"]
# Proposed changes
#thrIFR["Cnt"][0:7]=19.0,10.0, 9.0, 9.0, 9.0, 3.0, 5.0
#thrIFR["Cnt"][-1]=2.0

trains= dict()
for pref in prefs:
	trains[pref]= xtractSpikeTrains(dataDir=dDir,prefix=pref+"St",suff=".mat")
	#trains["Cnft"]= xtractTrains(dataDir=dDir,prefix=pref+"St",suff=".mat")
	nTrains=len(trains[pref])
	for cell in sc.arange(nTrains):
		st= trains[pref][cell]
		st.thrIFR= thrIFR[pref][cell]
		st.burstInds=st.burstIFR(st.thrIFR)

#
for pref in prefs:
	nTrains=len(trains[pref])
	f=gr.figure(figsize=(17,13)); gr.ioff()
	f.canvas.set_window_title(pref)
	cols=2; rows = sc.ceil(nTrains/cols)
	ax=list()
	for cell in sc.arange(nTrains):
		st= trains[pref][cell]
		ax.append(f.add_subplot(rows,cols,cell+1))
		aa=sc.array([0,st.spikeTimes.max()])
		ax[cell].plot(st.spikeTimes,st.ifr,'b.',label=st.fName[:-4]+' IFR (Hz)')
		thrStr= r"Threshold: (ISI,IFR)=(%g, %g)"%(thrISI[pref][cell],thrIFR[pref][cell])
		ax[cell].plot(aa, st.thrIFR*sc.array([1,1]),'r',label=thrStr)
		ax[cell].plot(st.bTimes, 30*sc.ones(len(st.bTimes)), 'k|',ms=10)
		ax[cell].set_xlabel("time (secs)")
		ax[cell].set_xlim(0,st.spikeTimes.max())
		ax[cell].set_ylim(0,sc.minimum(1.2*st.ifr.max(),200))
		#ax[cell].set_title(st.fName,fontsize=12)
		ax[cell].legend(fontsize=10)
		if len(st.bInds)>0:
			str0="%2.2f percent burst"%(100*len(st.bInds)/sc.float64(len(st.spikeTimes)))
			str1="burst percentage %g"%(100*len(st.bInds)/sc.float64(len(st.spikeTimes)))
			ax[cell].text(1,st.ifr.max(),str0,fontsize=10)
	f.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.98,wspace=0.1, hspace=0.175)
	f.savefig("./burstingFigs/bursting"+pref+".png")
	f.clf()
	gr.ion(); gr.draw()
	gr.close("all")

