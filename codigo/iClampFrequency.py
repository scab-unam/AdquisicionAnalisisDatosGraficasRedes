from spikeTrainAnalysis import *

srcDir="./freqFabiola/"
file1= "datosFrecD1L.csv"
file2= "datosFrecD2L.csv"

# Read generic csv file
def readMultipleColumnsCSV(delimiter=',', fileName ='.csv',nHeaderLines=1):
    data=dict()
    data['name'] = fileName
    #with open(fileName, newline='', encoding='utf-8','rb') as f:
    with open(fileName, newline='', encoding='utf-8') as f:
        #reader=csv.reader(f, delimiter=delimiter)
        reader=csv.reader(f)
        dataList=list()
        nRows=0
        nRead=0
        # Skip the first nHeaderLines
        nCols=0
        for row in reader:
            nRead=nRead+1
            if nRead== nHeaderLines:
                break

        for row in reader:
            #aa=sc.float64(row)
            dataList.append(row)
            print(row)
            nCols = sc.maximum(nCols,len(row))
            nRows=nRows+1

        dataList = sc.array(dataList).transpose()
        dd = list()
        for c in sc.arange(nCols):
            x= dataList[c]
            print(x)
            rr= 0
            for r in sc.arange(nRows):
                xx=dataList[c][r]
                if len(xx)>0:
                    rr=rr+1
                else:
                    print(xx)
            dd.append(sc.zeros(rr))
            for r in sc.arange(rr):
                xx=dataList[c][r]
                if len(xx)>0:
                    print(sc.float64(xx))
                    dd[c][r]=sc.float64(xx)
                    rr=rr+1
                else:
                    print(xx)

        data['values']= dd
        data['nRows']= nRows
        data['nCols']= nCols
        f.close()
    return data


d1= readMultipleColumnsCSV(delimiter=',', fileName =srcDir+file1,nHeaderLines=1)
d2= readMultipleColumnsCSV(delimiter=',', fileName =srcDir+file2,nHeaderLines=1)





