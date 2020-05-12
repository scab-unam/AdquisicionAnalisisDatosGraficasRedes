from transientPreAnalysis import *

if 0: 
    # Dictionary containing the names of recorded variables during the experiment
    dataNames= st.split(allData.children[0].userstr["S_waveNames"], ";")
    extractedData1=list()
    extractedData2=list()

    for nam in dataNames:
        if len(nam)>3:
            myInd=sc.int32(nam[4:]) -1
            str1= "w1=allData." + nam + "_%dF_min.data"%(myInd)
            str2= "w2=allData." + nam + ".data"%(myInd)
            #print("%s, %s"% (nam[4:],str1))
            exec(str2)
            #print("%s, %s"% (nam[4:],str2))
            exec(str1)
            extractedData1.append(w1)
            extractedData2.append(w2)
        else:
            print("Found empty string")

    dBleachData=sc.array(extractedData1)
    rawData=sc.array(extractedData2)
    nCells,nFrames= rawData.shape

    # Grafica exploratoria de muestras de trazas    
    sampleSize=10
    sampleInds = sc.random.random_integers(0,nCells,sampleSize)
    f=gr.figure(figsize=(13,11))
    ax1=f.add_subplot(311)
    ax2=f.add_subplot(312)
    ax3=f.add_subplot(313)
    gr.ioff();
    for n in sampleInds:
        ax1.plot(rawData[n])
        ax2.plot(dBleachData[n])
        ax3.plot(fMinData[n])
    ax1.plot(rawData.max(0),'k',lw=4,alpha=0.4)
    ax1.plot(rawData.min(0),'k',lw=4,alpha=0.4)

    gr.ion(); gr.draw()




