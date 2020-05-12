using MAT
using PyPlot

dataDir = "microcircuitsNetworks/"
dataFile = "duhne20140816-1.mat"
fName = dataDir* dataFile
data= matread(fName)
spikes=data["Spikes"]
nTPts=size(spikes,1)
nCells=size(spikes,2)

figure()
for n in 1:nCells
    plot([1:nTPts], n*spikes[:,n],"|")
end



