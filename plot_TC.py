import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
#plt.style.use('presentation')
plt.rcParams['legend.numpoints'] = 1
params = {
   'axes.labelsize': 12,
#   'text.fontsize': 8,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'lines.markeredgewidth' : 0,
   'lines.markersize' : 9,
   'figure.figsize': [4.2, 3.6]
}
plt.rcParams.update(params)
#plt.rcParams['text.usetex'] = True
kB = 1.38e-23
ps2s = 1e-12
A2m = 1e-10
e2J = 1.6e-19

def main():
    nstep = 8000
    deltaT = 0.005  # ps
    tmax = nstep * deltaT

    aveHCACF = np.zeros(nstep)
    aveTC = np.zeros(nstep)
    count = 0
    t = np.arange(nstep)*deltaT*1 # 100 ps
    print("Total time: {} ps".format(tmax))

    V = 42.432**3
    T = 40
    convert = e2J*e2J/ps2s/A2m
    ratio = convert/kB/T/T/V*1*deltaT
    for iseed in np.arange(0,10,1):
        print('seed: ',iseed+1)
        count +=1
        dirName = 'seed'+str(iseed+1)
        os.chdir(dirName)
        os.system("tail -n 8000 J0Jt.dat > data.txt")
        data = np.loadtxt('data.txt')
        HCACF = data[:,3:6]
        aveHCACF +=np.sum(HCACF,axis = 1)
        plt.figure(211)
        plt.plot(t,HCACF[:,0],'-',lw = 1,color='0.5',label='seed%d' %(iseed+1))
        plt.plot(t,HCACF[:,1],'-',lw = 1,color='0.5')
        plt.plot(t,HCACF[:,2],'-',lw = 1,color = '0.5')
        plt.xlim(0,300)
        plt.figure(212)
        TC = np.cumsum(HCACF, axis = 0) *ratio
        plt.plot(t,TC[:,0],'-',lw = 1, color='0.8')
        plt.plot(t,TC[:,1],'-',lw = 1, color='0.8')
        plt.plot(t,TC[:,2],'-',lw = 1, color = '0.8')
        aveTC +=np.sum(TC,axis=1)
        os.chdir('..')
    aveHCACF /=(count*3.0)
    aveTC /=(count*3.0)

    # save heat current correlation function
    plt.figure(211)
    plt.xlim(-0,tmax)
    plt.plot(t,aveHCACF, 'b-',lw=1.5)
    plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    plt.xlabel('Time (ps)')
    plt.ylabel('HCACF (a.u.)')
    plt.subplots_adjust(left=0.25,bottom=0.15, right = 0.95, top = 0.9)
    #plt.legend()
    plt.savefig('PbS_Ligand_HCACF_300K.png',dpi=300)

    # save the Tc
    plt.figure(212)
    plt.plot(t,aveTC, 'b-',lw=1.5)
    N = 10
    #temp=np.convolve(aveTC, np.ones((N))/N, mode='same')
    #plt.plot(t,temp, 'r-',lw=1.5)
    middle = int(len(aveTC)/2 )
    meanTC = np.mean(aveTC[middle:])
    plt.plot([t[middle],t[-1]],[meanTC,meanTC],'r-',lw=1.5)
    plt.text(2,1,'TC=%.3f W/mK'%meanTC)
    print('thermal conductivity:', meanTC)
    plt.xlim(-0,tmax)
    plt.ylim(-0.4,1.5)
    #plt.ylim(0,8)
    plt.xlabel('Time (ps)')
    plt.ylabel('Thermal conductivity (W/mK)')
    plt.subplots_adjust(left=0.25,bottom=0.15, right = 0.95, top = 0.9)
    plt.savefig('PbS_Ligand_TC_300K.png',dpi=300)
    out=np.column_stack((t,aveHCACF,aveTC))
    np.savetxt('PbS_GK_300K.dat',out)
    plt.show()

if __name__=="__main__":
    main()
