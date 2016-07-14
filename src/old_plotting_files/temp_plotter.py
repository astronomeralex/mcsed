import numpy as np
import matplotlib.pylab as plt

myidx = 2

mine = np.load('/users/hunterbrooks/astronomy/results/1.2c/predictions.dat.npy')[:,myidx]
myerr = np.load('/users/hunterbrooks/astronomy/results/1.2c/errors.dat.npy')



his = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/data_house/alex.dat',
                 dtype='str')[:,18]
his_temp = np.zeros(np.shape(his))
for i in range(0, len(his)):
    try:
        his_temp[i] = float(his[i])
    except ValueError:
        his_temp[i] = 0
his = np.array(his_temp,copy=True)


        
hiserr = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/data_house/alex.dat',
                    dtype='str')[:,19]
hiserr_temp = np.zeros(np.shape(hiserr))
for i in range(0, len(hiserr)):
    try:
        hiserr_temp[i] = float(hiserr[i])
    except ValueError:
        hiserr_temp[i] = 0
hiserr = np.array(hiserr_temp,copy=True)





hisfield = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/data_house/alex.dat',
                      dtype='str')[:,7]
hisfield_temp = np.zeros(np.shape(hisfield))
for i in range(0, len(hisfield)):
    try:
        hisfield_temp[i] = float(hisfield[i])
    except ValueError:
        hisfield_temp[i] = 0
hisfield = np.array(hisfield_temp,copy=True)




hisid = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/data_house/alex.dat',
                   dtype='str')[:,8]
hisid_temp = np.zeros(np.shape(hisid))
for i in range(0, len(hisid)):
    try:
        hisid_temp[i] = float(hisid[i])
    except ValueError:
        hisid_temp[i] = 0
hisid = np.array(hisid_temp,copy=True)



myid = np.loadtxt('/users/hunterbrooks/astronomy/results/1.2c/superInfo.dat')[:,4] 
myfield = np.loadtxt('/users/hunterbrooks/astronomy/results/1.2c/superInfo.dat')[:,3] 
errors = np.load('/users/hunterbrooks/astronomy/results/1.2c/errors.dat.npy')

'''
EB vs DELTA VS:
     surface density(mass/pi*r**2),
     stellar mass,
     ellipticity(from b/a),
     effective radius (=size) (18)
     O/H. =12
'''

a = []
b = []
c = []
d = []
e = []

for i in range(0, np.shape(mine)[0]):
    if mine[i] != 0:
        print '----'
        print myfield[i], myid[i]
        
        pos_id_idx = np.where(hisid == myid[i])[0]
        
        if len(pos_id_idx) == 1:
            his_id_idx = int(pos_id_idx)
            
        if len(pos_id_idx) == 0:
            continue
            
        else:
            pos_field_idx = np.where(hisfield == myfield[i])[0]
            
            for ii in range(0, len(pos_field_idx)):
                for iii in range(0, len(pos_id_idx)):
                    if pos_field_idx[ii] == pos_id_idx[iii]:
                        his_id_idx = int(pos_id_idx[iii])
                        
        print hisfield[his_id_idx], hisid[his_id_idx]
            
        assert myfield[i] == hisfield[his_id_idx]
        assert myid[i] == hisid[his_id_idx]
        
        if his[his_id_idx] != 0:
            a.append(his[his_id_idx])
            #b.append(hiserr[his_id_idx])
            c.append(mine[i])
            d.append(myerr[i,myidx,0])
            e.append(myerr[i,myidx,1])
            
# since errors are just percentile, need difference from median
d = np.asarray(c)-np.asarray(d)
e = np.asarray(e)-np.asarray(c)

plt.errorbar(c, a,  xerr=[d,e], fmt='o', #xerr=b,
             color='b', capsize=0, alpha=0.50)
             
plt.xlabel(r'$E_b$')
plt.ylabel('Half Light Radius [kpc]')
plt.title(r'$E_b$ vs Half Light Radius for 34 objects')
plt.show()