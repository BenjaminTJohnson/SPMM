import numpy as np
import glob
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os, time

NUM_PROCESSES = 4
pcnt = 0
infiles = sorted(glob.glob("domain*.out"))
#infiles = sorted(glob.glob("spmelt3c_f0100000.out"))
print(infiles)
ld = len(infiles)-1
skip = 1 #arbitrary
#ld = 1 #temporary for testing

minr = [0,int(ld*1/4),int(ld*2/4),int(ld*3/4)]
maxr = [int(ld*1/4)+1,int(ld*2/4)+1,int(ld*3/4)+1,int(ld)+1]
#mx = 30
ms = 5.0

#get shapes
print('h1')
infile = infiles[0]
print('h1a', infile)
t = np.loadtxt(infile)
print('h2')
x = t[:,1]
y = t[:,2]
z = t[:,3]
mx = max(max(x),max(y),max(z))
#ms = round(mx**(1.0/3.0)/2.0)
print('h3', ms )


def fileloop(minr,maxr,infiles):
    fcnt = 0

    for ifile in range(minr,maxr,skip):
        fcnt = fcnt + 1
        infile = infiles[ifile]
        outfile = infile[:-4]+'.png'

        t = np.loadtxt(infile)
        
        x = t[:,1]
        y = t[:,2]
        z = t[:,3]

        ctype = t[:,4].astype(np.int)
        qqw = np.nonzero(ctype == -1)
        qqi = np.nonzero(ctype == 1)
        qqo = np.nonzero(ctype == 0)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111,projection='3d')
        m = '.'
        lenqqo=np.shape(qqo)[1]
        lenqqi=np.shape(qqi)[1]
        lenqqw=np.shape(qqw)[1]


        ax.scatter(x[qqw],y[qqw],z[qqw],c='r',marker=m,lw=0,s=ms)
        ax.scatter(x[qqi],y[qqi],z[qqi],c='b',marker=m,lw=0,s=ms)
        ax.set_xlim(0,mx)
        ax.set_ylim(0,mx)
        ax.set_zlim(0,mx)
        ax.view_init(-30, -40)

        plt.gca().set_aspect('auto', adjustable='box')


        fig.savefig(outfile)
        print('output filename: ', outfile)
        plt.close(fig)  
     

children = []

start_time = time.time()
print('start:',os.getpid())
cnt = 0
for process in range(NUM_PROCESSES):
    pid = os.fork()
    print('basepid:',os.getpid())

    if pid:
        cnt = cnt + 1
        children.append(pid)
    else:
        print('childpid:',os.getpid(), cnt)
        fileloop(minr[cnt],maxr[cnt],infiles)
        os._exit(0)
        
for i, child in enumerate(children):
    os.waitpid(child, 0)
    

