import numpy as np
import matplotlib.pylab as plt

myidx = 1

mine = np.load('/users/hunterbrooks/astronomy/results/1.2c/predictions.dat.npy')[:,myidx]
myerr = np.load('/users/hunterbrooks/astronomy/results/1.2c/errors.dat.npy')



his = np.load('/users/hunterbrooks/astronomy/results/1.2c/predictions.dat.npy')[:,4]
hiserr = np.load('/users/hunterbrooks/astronomy/results/1.2c/errors.dat.npy')[:,4]



a = []
b = []
c = []
d = []
e = []

for i in range(0, np.shape(mine)[0]):
    if mine[i] != 0:
        a.append(his[i])
        b.append(hiserr[i])
        c.append(mine[i])
        d.append(myerr[i,myidx,0])
        e.append(myerr[i,myidx,1])
            
# since errors are just percentile, need difference from median
d = np.asarray(c)-np.asarray(d)
e = np.asarray(e)-np.asarray(c)

plt.errorbar(c, a,  xerr=[d,e], fmt='o', #xerr=b,
             color='b', capsize=0, alpha=0.50)
             
plt.xlabel(r'$\delta$')
plt.ylabel('log10(Stellar Mass) [M_sol]')
plt.title(r'$\delta$ vs log10(Stellar Mass)')
plt.show()