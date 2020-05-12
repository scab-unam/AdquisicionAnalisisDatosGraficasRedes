# Visibility graph analysis on bursting time series
using PyCall
@pyimport matplotlib.pyplot as gr
using MAT

dataDir= "trenesBursts/aleSpikeTimes/"
fName = "CntSt1.mat"
file= matopen(dataDir * fName)
pigas= read(file)["savevar"][1]
ifr = zeros(size(pigas))
ifr[2:end]= 1./(pigas[2:end]-pigas[1:end-1])
gr.figure()
gr.ioff()
gr.plot(pigas,ifr, linestyle="--", color="black", linewidth=3,alpha=0.7)
gr.ion()
gr.draw()

function pVisibleAB(A,B)
	Int32(p[2] < (p[1] - A[1]) * (B[2]-A[2]) / (B[1]-A[1]))
end



