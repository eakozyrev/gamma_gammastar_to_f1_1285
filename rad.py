from __future__ import print_function
from math import sin
from array import array
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    rad = []
    en = []
    drad = array( 'd' )
    den = array( 'd' )
    i = 0
    file_of = open("results/rad.dat",'w')
    with open('testgg/rad.dat') as fp:
        en0,cr,dcr,cr0,dcr0 = np.loadtxt(fp,usecols=(0,1,2,3,4),unpack=True)
        for el in en0:
            den.append(0)
            if i < 6:
                cr[i] =  cr[i]-cr[i+1]
                cr0[i] =  cr0[i]-cr0[i+1]
                file_of.write(f'{en0[i]} {en0[i]} {float(cr[i])/float(cr0[i])}\n');
            rad.append(float(cr[i])/float(cr0[i]))
            drad.append(float(dcr[i])/float(cr0[i]))
            if i < 6: en.append(en0[i]/2.+en0[i+1]/2.)
            else:  en.append(100)
            i+=1
    print(en)
    plt.errorbar(en,rad,drad,marker='s', mfc='red',mec='blue', ms=4, mew=4)
    plt.xlabel("Q^{2}, GeV^2",fontsize=20)
    plt.ylabel('rad corr.',fontsize=20)
    plt.show()
        
    #en = np.array(en)
    #f = TGraphErrors(en,den,rad,drad)
    #f.Draw("AP")
        # plt.errorbar(en,rad)
