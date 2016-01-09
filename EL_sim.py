# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 18:08:39 2015

@author: patrick
"""

###-------------------------------------------------------------------------------
# Name:        EL_sim
# Purpose:     Simulate Early-Late sampling regeimes to observe the probability that an allele with a consistent average
# effect in both years is found to be significant in both years.  
#       
#-------------------------------------------------------------------------------

from random import *
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats as sp
import csv

OUTPUT_FILE="EL_sim.csv"
out1=open(OUTPUT_FILE,'wb')
out1x=csv.writer(out1,delimiter=",",dialect="excel")
out1x.writerow(["rep","A","rd","LL","dp","pval","psig"])

Reps=100

p_cutoff = 4.62448137005e-06

p=0.5 #allele frequency
a=9 #average effect
Ve=10
PopMean=50

N=10000 #Pop Size
B=100 #Bulk Size

DD=80 #Dry down
LS=30 #Life Span

LvE=0.01
LvL=0.01 #Library Prep vars

RdE=20
RdL=20


for A in range(a-2,a+3):
    for rd in range(RdE-10,RdE+30,10):
        for rep in range(0,Reps):
            
            vE=LvE+(1.0/float(rd))
            vL=LvL+(1.0/float(rd))            
            wE=1/vE
            wL=1/vL
            
            Phenos=[]
            Genos=[]
            for i in range(N):
                rx1=random()
                rx2=random()
                ac=0
                if rx1<p:
                    ac+=1
                    if rx2<p:
                        ac+=1        
                else:
                    if rx2<p:
                        ac+=1
                        
                z=ac*A+normalvariate(PopMean, Ve)
            
                Phenos.append(z)
                Genos.append(ac)
            
            E=[] 
            for j, ind in enumerate(Phenos):
                if ind < np.percentile(Phenos,10):
                    E.append(Genos[j])
            
            E=sample(E,B)
            
            L=[]
            for j, ind in enumerate(Phenos):
                if ind > np.percentile(Phenos,90):
                    L.append(Genos[j])
                    
            L=sample(L,B)
            
            pE=float(sum(E))/float((2*len(E)))
            pL=float(sum(L))/float((2*len(L)))
            
            dp=pE-pL
            
            zE=2*math.asin(pE**0.5)
            zL=2*math.asin(pL**0.5)
            
            ze=normalvariate(zE,vE**0.5)
            zl=normalvariate(zL,vL**0.5)
            
            print ze,zl
            print zE,zL
            
            u=(ze*wE+zl*wL)/(wE+wL)
            
            LL=(-((ze-u)**2)/(2*vE))+(-((zl-u)**2)/(2*vL))
            
            pval=1-sp.chi2.cdf(-(2*LL),1)
            if pval<p_cutoff:
                psig=1
            else:
                psig=0
            print u, LL
            print pL, pE
            print pval
            print dp
            print rep,A,rd
            out1x.writerow([rep,A,rd,-LL,dp,pval,psig])


#
#plt.hist(Phenos)
#plt.show()



#mu, sigma = 100, 15
#x = mu + sigma * np.random.randn(10000)
## the histogram of the data
#n, bins, patches = plt.hist(x, 50, normed=1, facecolor='g', alpha=0.75)
#
#
#plt.xlabel('Smarts')
#plt.ylabel('Probability')
#plt.title('Histogram of IQ')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
#plt.grid(True)
#plt.show()
    


