# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 09:51:44 2015

@author: patrick
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 16:26:57 2015

@author: patrick
"""

import csv

#INDIR="/Volumes/avery/Research/EarlyLate/"
#OUTDIR="/Volumes/avery/Research/EarlyLate/"
INDIR="/Users/monnahap/Documents/Research/EarlyLate/"
OUTDIR="/Users/monnahap/Documents/Research/EarlyLate/"

INPUT_FILE="EL_Likelihoods_IndTests_20170406-1007.csv"

p_B_im13_cutoff = 1.55238590731e-05
p_B_br13_cutoff = 3.69770185093e-07
p_B_q14_cutoff = 4.6588917357e-06
p_B_im14_cutoff = 5.43576299401e-06
p_B_q13_cutoff = 0.000256908772079


sigfile=open(OUTDIR+INPUT_FILE[:-4]+"_sig.csv","wb")
sigwriter=csv.writer(sigfile,delimiter=",",dialect='excel')    

nonsigfile=open(OUTDIR+INPUT_FILE[:-4]+"_nonsig.csv","wb")
nonsigwriter=csv.writer(nonsigfile,delimiter=",",dialect='excel')    

sites=0
numsigB_im13=0
numsigB_q13=0
numsigB_q14=0
numsigB_im14=0
numsigB_br13=0

with open(INDIR+INPUT_FILE,"rb") as sites_file:
    scaffs=[]
    snps=[]
    counts=[]
    for i,site in enumerate(sites_file):
        sites+=1
        if sites%10000==0:
            print sites
        sigB=0
        site=site.strip("\n")
        site=site.strip("\r")
        site=site.split(",")
        if i == 0:      
            
            for k,j in enumerate(site):
                print k,j
            sigwriter.writerow(site+["sigBim13","sigBq13","sigBbr13","sigBim14","sigBq14"])
            nonsigwriter.writerow(site+["sigBim13","sigBq13","sigBbr13","sigBim14","sigBq14"])
        
        else:
            scaff=site[0]
            try:
                pos=float(site[1])
            except IndexError:
                pass
            try:
                p_B=float(site[3])               
            except (ValueError, IndexError):
                p_B=99
                sigB="-"
            try:
                pop=site[4]
            except (ValueError, IndexError):
                pass
            
            
            
            if pop=="IM13":
                if p_B <= p_B_im13_cutoff:
                    sigB=1
                    numsigB_im13+=1
            if pop=="IM14":
                if p_B <= p_B_im14_cutoff:
                    sigB=1
                    numsigB_im14+=1
            if pop=="Q13":
                if p_B <= p_B_q13_cutoff:
                    sigB=1
                    numsigB_q13+=1
            if pop=="Q14":
                if p_B <= p_B_q14_cutoff:
                    sigB=1
                    numsigB_q14+=1
            if pop=="BR13":
                if p_B <= p_B_br13_cutoff:
                    sigB=1
                    numsigB_br13+=1

            if sigB==1:
                sigwriter.writerow(site+[sigB])                
            else:
                nonsigwriter.writerow(site+[sigB]) 
                
print "Number of B sig IM13 = ", numsigB_im13
print "Number of B sig IM14 = ", numsigB_im14
print "Number of B sig Q13 = ", numsigB_q13
print "Number of B sig Q14 = ", numsigB_q14
print "Number of B sig BR = ", numsigB_br13

                
    