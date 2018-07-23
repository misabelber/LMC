import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm as cm

dat = pd.read_table("/afs/ciemat.es/user/b/bernardos/GitHub/LMC/results/Correlations1TeVDM.dat",sep=" ")

fig = plt.figure()
ax1 = fig.add_subplot(111)
cmap = cm.get_cmap('jet',30)
cax = ax1.imshow(dat,cmap=cmap)
plt.title('Correlation Matrix of Baryonic Sources')
labels=dat.columns
plt.locator_params(axis='x', nbins=17)
plt.locator_params(axis='y', nbins=17)
ax1.set_xticklabels(labels,fontsize=7,rotation='45')
ax1.set_yticklabels(labels,fontsize=7)
fig.colorbar(cax)
locs = ([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
plt.xticks(locs,labels)
plt.yticks(locs,labels)
plt.xlim(-0.5,16.5)
plt.ylim(16.5,-0.5)
plt.show()
