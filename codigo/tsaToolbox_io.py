# analysisToolbox_io.py
import csv
import numpy as np

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

# -----------------------------------------------------------------
# read Matlab files from a list
# -----------------------------------------------------------------
def readMatlabFileList(dataDir,files,varName):
    nFiles=len(files)
    data= list()
    for n in sc.arange(nFiles):
        fName = files[n]
        ppp = sc.squeeze(io.loadmat(dataDir + fName)[varName][0][0])
        data.append(ppp)
    return data


# Read generic csv nFiles (require csv module)

def readColumnCSV(dataDir='./',  fileName ='Neuron.csv', delimiter=',', nHeaderRows=0):
    data=dict()
    data['name'] = fileName
    with open(dataDir+fileName,'rb') as f:
        reader=csv.reader(f, delimiter=delimiter)
        dataList=list()
        nRows=0
        for row in reader:
            if nRows>nHeaderRows-1:
                aa=np.float32(row)
                dataList.append(aa)
            nRows=nRows+1

        #for row in reader:
        #    aa=sc.float32(row)
        #    dataList.append(aa)
        #    nRows=nRows+1

        data['values']= np.array(dataList).squeeze()
        data['nRows']= nRows
        f.close()
    return data


def readMultipleColumnsCSV(delimiter=',', fileName ='.csv',nHeaderLines=1,nLabelCols=1):
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
        # Put rows into the list
        for row in reader:
            dataList.append(row)
            print(row)
            nCols = np.maximum(nCols,len(row))
            nRows=nRows+1

        dataList = np.array(dataList).transpose()
        dd = list()
        for c in range(nLabelCols,nCols):
            x= dataList[c]
            print(x)
            rr= 0
            for r in range(nRows):
                xx=dataList[c][r]
                if len(xx)>0:
                    rr=rr+1
                else:
                    print(xx)
            dd.append(np.zeros(rr))
            for r in range(rr):
                xx=dataList[c][r]
                if len(xx)>0:
                    print(np.float64(xx))
                    dd[c][r]=np.float64(xx)
                    rr=rr+1
                else:
                    print(xx)

        data['values']= dd
        data['nRows']= nRows
        data['nCols']= nCols
        f.close()
    return data
