# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 16:26:57 2015

@author: patrick
"""

import csv

INDIR="/Users/patrick/Documents/Research/EarlyLate/"
OUTDIR="/Volumes/TOSHIBA EXT/EarlyLate/"
INPUT_FILE="EL_Likelihoods_both20160209-2211.csv"

p_Y_cutoff = 0.000112266628638
p_B_cutoff = 0.00015940722695
p_P_cutoff = 0.060793629058


sigfile=open(OUTDIR+INPUT_FILE[:-4]+"_sig.csv","wb")
sigwriter=csv.writer(sigfile,delimiter=",",dialect='excel')    

nonsigfile=open(OUTDIR+INPUT_FILE[:-4]+"_nonsig.csv","wb")
nonsigwriter=csv.writer(nonsigfile,delimiter=",",dialect='excel')    
 
sites=0
numsigY=0
numsigB=0
numsigP=0
numsig=0
with open(OUTDIR+INPUT_FILE,"rb") as sites_file:
    scaffs=[]
    snps=[]
    counts=[]
    for i,site in enumerate(sites_file):
        sites+=1
        if sites%10000==0:
            print sites
        sigY=0
        sigP=0
        sigB=0
        site=site.strip("\n")
        site=site.strip("\r")
        site=site.split(",")
        if i == 0:            
            for k,j in enumerate(site):
                print k,j
            sigwriter.writerow(site+["sigY","sigB","sigP"])
            nonsigwriter.writerow(site+["sigY","sigB","sigP"])
        else:
            try:
                p_Y=site[35]
                p_P=site[36]
                p_B=site[37]
            except IndexError:
                print site #stopped here
            scaff=site[0]
            pos=float(site[1])
            try:
                p_Y=float(site[35])                
            except (ValueError, IndexError):
                p_Y=99
                sigY="-"                
            try:
                p_B=float(site[37])               
            except (ValueError, IndexError):
                p_B=99
                sigB="-"
            try:
                p_P=float(site[36])              
            except (ValueError, IndexError):
                p_P = 99
                sigP="-"
            if (p_Y < p_Y_cutoff or p_B < p_B_cutoff or p_P < p_P_cutoff):
                numsig+=1
                if p_Y < p_Y_cutoff:
                    sigY=1
                    numsigY+=1
                if p_B < p_B_cutoff:
                    sigB=1
                    numsigB+=1
                if p_P < p_P_cutoff:
                    sigP=1
                    numsigP+=1
#                if sigY == 1:
#                    print "here"
                sigwriter.writerow(site+[sigY,sigB,sigP])
                
            else:
                nonsigwriter.writerow(site+[sigY,sigB,sigP]) 
                
print "Number of Y sig = ", numsigY
print "Number of B sig = ", numsigB
print "Number of P sig = ", numsigP
print "Number of total hits = ", numsig  
                
    