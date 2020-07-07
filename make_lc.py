import numpy as np
import matplotlib as mp
import glob
import matplotlib.pyplot as plt


#~~~~~~~~~~~~~~~~~~~~~~for plotting purposes~~~~~~~~#
mp.rcParams['font.family']='serif'
mp.rcParams['xtick.major.size']=10
mp.rcParams['xtick.major.width']=2
mp.rcParams['xtick.minor.size']=7
mp.rcParams['xtick.minor.width']=2
mp.rcParams['ytick.major.size']=10
mp.rcParams['ytick.major.width']=2
mp.rcParams['ytick.minor.size']=7
mp.rcParams['ytick.minor.width']=2
mp.rcParams['axes.linewidth']=1.5
mp.rcParams['xtick.labelsize']=36
mp.rcParams['ytick.labelsize']=36
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

files=sorted(glob.glob('*dat'))
agn=[]
s1=[]
s2=[]
s3=[]
s4=[]
s5=[]
s6=[]
julian=[]
agn_err=[]
s1_err=[]
s2_err=[]
for i in range(0,len(files)):
    jd=files[i]
    julian.append(float(jd[0:12])-2458850)
    bkg=np.genfromtxt(files[i],unpack=True,usecols=0)
    sources,source_err=np.genfromtxt(files[i],unpack=True,usecols=(3,4))
    agn.append(sources[0]-sources[6]+1)
    agn_err.append(source_err[0])
    s1.append(sources[6]-sources[5])
    s1_err.append(source_err[1])
    s2.append(sources[0]-sources[5]+2)
    s2_err.append(source_err[2])
    s3.append(sources[3])
    s4.append(sources[4])
    s5.append(sources[5])
    s6.append(sources[6])


fig,ax=plt.subplots()
plt.plot(julian,s1,'ko',markersize=15,label='star differential')
plt.plot(julian,agn,'bo',markersize=15,label='AGN differential 1')
plt.plot(julian,s2,'go',markersize=15,label='AGN differential 2')
plt.errorbar(julian,s1,yerr=s1_err,capsize=5,fmt=' ')
plt.errorbar(julian,agn,yerr=s1_err,capsize=5,fmt=' ')
plt.errorbar(julian,s2,yerr=s2_err,capsize=5,fmt=' ')
plt.ylabel("Magnitudes (instrumental)",fontsize=36)
plt.xlabel("Julian Dates (2458850+)",fontsize=36)
#plt.title("Deviation in lag values with respect to the variability",fontsize=36)
ax.tick_params(axis="both",which='minor',direction="in")
ax.tick_params(axis="both",which='major',direction="in")
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.minorticks_on()
ax.legend(fontsize=30)
plt.show()
