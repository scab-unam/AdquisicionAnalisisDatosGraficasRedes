
# coding: utf-8

# # Visibility graphs from spike trains
# ### Verónica Alejandra Cáceres$^1$, José Bargas $^1$ & Marco Arieli Herrera-Valdez$^{2}$
# #### $^1$Instituto de Fisiología Celular & $^2$ Facultad de Ciencias, Universidad Nacional Autónoma de México
# 
# Created: 20151131 MAHV
# 
# **Description** 
# 
# Visibility graphs from the instantaneous firing rate of a bursting neron are constructed to explore a possible classification scheme for spike trains. The analysis is implemented in Python 3 (www.python.org). 

# ## Exploration: Visibility graph from a single spike train
# 
# Import the necessary modules. Many of the functions created to perform burst analyisis for these data sets are contained in the python module burstDetection.py

# In[1]:

#get_ipython().magic('matplotlib inline')
from burstDetection import *


# Open one data file containing the spike times and preprocess the data to obtain instantaneous firing rates
# In[2]:

dataDir='./trenesBursts/aleSpikeTimes/'
fName = "CntSt1.mat"
pigas = sc.squeeze(io.loadmat(dataDir + fName)["savevar"][0][0])
print("Found %d spikes"%(len(pigas)))

# Now obtain interspike intervals and instantáneous firing rates
# In[3]:
isi = sc.zeros(len(pigas)); ifr = sc.zeros(len(pigas)) #Note time stamps are in milliseconds
isi[1:]= pigas[1:]-pigas[:-1]
ifr[1:]= 1/isi[1:]

# Calculate the visility graph obtained from different segments of the train. 
# The resulting graphs will be put into a list called $g$.
# Then, plot the visibility of the spike train ifr trace in segments

# In[6]:

g=list()
nSeg=8; nPtsSeg= len(pigas)/nSeg
print("Dividing the spike train into %d segments with %d points each"%(nSeg,nPtsSeg))
for n in sc.arange(nSeg):
    a=n*nPtsSeg; 
    if n==nSeg-1:
        b=-1
    else:
        b=a+nPtsSeg
    print("Segment 1 between %g and %g secs"%(pigas[a],pigas[b]))
    g.append(visibilityGraph(pigas[a:b],ifr[a:b]))
#
nodeSize=10; 
nodeColor1='blue'; nodeColor2='yellow'
#
f=gr.figure(figsize=(17,11))
ax=list(); 
if nSeg>3:
    rows=2; cols=sc.ceil(nSeg/2); 
else:
    rows=1; cols=nSeg
gr.ioff();
for n in sc.arange(nSeg):
    a=n*nPtsSeg; 
    if n==nSeg-1:
        b=-1
    else:
        b=a+nPtsSeg
    ax.append(f.add_subplot(rows,cols,n+1))
    #nx.draw(g[n],ax=ax[n],pos=nx.shell_layout(g[n]), node_size=nodeSize,node_color=nodeColor1,alpha=0.5)
    #nx.draw(g[n],ax=ax[n],pos=nx.graphviz_layout(g[n]), node_size=nodeSize,node_color=nodeColor1,alpha=0.5)
    #nx.draw(g[n],ax=ax[n],pos=nx.spectral_layout(g[n]),node_size=nodeSize,node_color=nodeColor1,alpha=0.5)
    #nx.draw(g[n],ax=ax[n],pos=nx.random_layout(g[n]), node_size=nodeSize,node_color=nodeColor1,alpha=0.5)
    nx.draw(g[n],ax=ax[n],pos=nx.spring_layout(g[n]), node_size=nodeSize,node_color=nodeColor1,alpha=0.5)
    #nx.draw(g[n],ax=ax[n],pos=nx.circular_layout(g[n]), node_size=nodeSize,node_color=nodeColor1,alpha=0.5)
    #nx.draw_networkx_labels(g[n],ax=ax[n],pos=nx.circular_layout(g[n]))
    ax[n].set_title(fName+" %g-%g secs"%(pigas[a],pigas[b]) )
gr.ion(); gr.draw()


# ## Comparison between different spike trains recorded in control and lesioned animals 

# In[8]:

ctrFiles=getFileList(dataDir, prefix="CntSt", suffix=".mat")
lesFiles=getFileList(dataDir, prefix="LesSt", suffix=".mat")


# Use the script to pre-process the data and obtain lists that can be used to produce several visibility graphs

# In[34]:

ctrData= preprocessTrainsFromMatFiles(ctrFiles)
lesData= preprocessTrainsFromMatFiles(lesFiles)


# Create a function that calculates visibility graphs from the spike trains in a list as pre-processed by the function _preprocessTrainsFromMatFiles_.

# In[62]:




# Now get the visibility graphs from the spike trains 

# In[58]:

gCtr=visibilityFromSpikeTrains(ctrData["spikeTrains"], ctrData["IFRs"],a=0, b=6000)
gLes=visibilityFromSpikeTrains(lesData["spikeTrains"], lesData["IFRs"],a=0, b=6000)


# The visibility graphs from the spike trains recorded from control animals are

# In[59]:

fCtr=showGraphList(gCtr, ctrFiles)


# The visibility graphs from the spike trains recorded from lesioned animals are

# In[60]:

fLes=showGraphList(gLes, lesFiles)


# ### Analysis of the topological properties of the graphs

# In[ ]:



