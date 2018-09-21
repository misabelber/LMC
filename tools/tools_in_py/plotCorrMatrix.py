import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm as cm

dat = pd.read_table("/afs/ciemat.es/user/b/bernardos/GitHub/LMC/results/CMatrix_rebin_0.1x100_Pointin5deg_OnlyExtended.dat",sep=" ")

fig = plt.figure()
ax1 = fig.add_subplot(111)
cmap = cm.get_cmap('jet',30)
cax = ax1.imshow(dat,cmap=cmap)
plt.title('Correlation Matrix of Extended Baryonic Sources')
labels=dat.columns
plt.locator_params(axis='x', nbins=6)
plt.locator_params(axis='y', nbins=6)
ax1.set_xticklabels(labels,fontsize=10,rotation='45',ha='right')
ax1.set_yticklabels(labels,fontsize=10)
fig.colorbar(cax)
locs = ([0,1,2,3,4,5])
plt.xticks(locs,labels)
plt.yticks(locs,labels)
plt.xlim(-0.5,5.5)
plt.ylim(5.5,-0.5)
plt.show()
