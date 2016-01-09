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

INDIR="/Volumes/TOSHIBA EXT/EarlyLate/"
OUTDIR="/Volumes/TOSHIBA EXT/EarlyLate/"
INPUT_FILE="EL_Likelihoods_PerPop_201320151118-1514.csv"

p_Y_im_cutoff = 0
p_YCB_im_cutoff = 0
p_B_im_cutoff = 4.62448137005e-06
p_Y_br_cutoff = 0
p_B_br_cutoff = 1.60241798643e-06
p_Y_q_cutoff = 0
p_B_q_cutoff = 0.000138586071041
p_YCB_q_cutoff = 0

sigfile=open(OUTDIR+INPUT_FILE[:-4]+"_sig.csv","wb")
sigwriter=csv.writer(sigfile,delimiter=",",dialect='excel')    

nonsigfile=open(OUTDIR+INPUT_FILE[:-4]+"_nonsig.csv","wb")
nonsigwriter=csv.writer(nonsigfile,delimiter=",",dialect='excel')    

sites=0
numsigY_im=0
numsigB_im=0
numsigYCB_im=0
numsig_im=0
numsigY_q=0
numsigB_q=0
numsigYCB_q=0
numsig_q=0
numsigY_br=0
numsigB_br=0
numsig_br=0
with open(INDIR+INPUT_FILE,"rb") as sites_file:
    scaffs=[]
    snps=[]
    counts=[]
    for i,site in enumerate(sites_file):
        sites+=1
        if sites%10000==0:
            print sites
        sigY_im=0
        sigYCB_im=0
        sigB_im=0
        sigY_q=0
        sigYCB_q=0
        sigB_q=0
        sigY_br=0
        sigB_br=0
        site=site.strip("\n")
        site=site.strip("\r")
        site=site.split(",")
        if i == 0:      
            
            for k,j in enumerate(site):
                print k,j
            sigwriter.writerow(site+["sigYim","sigBim","sigYCBim","sigYq","sigBq","sigYCBq","sigYbr","sigBbr"])
            nonsigwriter.writerow(site+["sigYim","sigBim","sigYCBim","sigYq","sigBq","sigYCBq","sigYbr","sigBbr"])
        
        else:
            scaff=site[0]
            pos=float(site[1])
            try:
                p_Y_im=float(site[32])                
            except (ValueError, IndexError):
                p_Y_im=99                
            try:
                p_B_im=float(site[35])               
            except (ValueError, IndexError):
                p_B_im=99
            try:
                p_YCB_im=float(site[44])              
            except (ValueError, IndexError):
                p_YCB_im = 99
            try:
                p_Y_q=float(site[38])                
            except (ValueError, IndexError):
                p_Y_q=99                
            try:
                p_B_q=float(site[41])               
            except (ValueError, IndexError):
                p_B_q=99
            try:
                p_YCB_q=float(site[47])              
            except (ValueError, IndexError):
                p_YCB_q = 99
            try:
                p_Y_br=float(site[26])                
            except (ValueError, IndexError):
                p_Y_br=99                
            try:
                p_B_br=float(site[29])               
            except (ValueError, IndexError):
                p_B_br=99
            
            
            if p_Y_im < p_Y_im_cutoff:
                sigY_im=1
                numsigY_im+=1
            if p_B_im < p_B_im_cutoff:
                sigB_im=1
                numsigB_im+=1
            if p_YCB_im < p_YCB_im_cutoff:
                sigYCB_im=1
                numsigYCB_im+=1
            if p_Y_q < p_Y_q_cutoff:
                sigY_q=1
                numsigY_q+=1
            if p_B_q < p_B_q_cutoff:
                sigB_q=1
                numsigB_q+=1
            if p_YCB_q < p_YCB_q_cutoff:
                sigYCB_q=1
                numsigYCB_q+=1
            if p_Y_br < p_Y_br_cutoff:
                sigY_br=1
                numsigY_br+=1
            if p_B_br < p_B_br_cutoff:
                sigB_br=1
                numsigB_br+=1

            if any(k == 1 for k in [sigY_im,sigB_im,sigYCB_im,sigY_q,sigB_q,sigYCB_q,sigY_br,sigB_br]):
                sigwriter.writerow(site+[sigY_im,sigB_im,sigYCB_im,sigY_q,sigB_q,sigYCB_q,sigY_br,sigB_br])                
            else:
                nonsigwriter.writerow(site+[sigY_im,sigB_im,sigYCB_im,sigY_q,sigB_q,sigYCB_q,sigY_br,sigB_br]) 
                
print "Number of Y sig IM = ", numsigY_im
print "Number of B sig IM = ", numsigB_im
print "Number of YIB sig IM = ", numsigYCB_im
print "Number of Y sig Q = ", numsigY_q
print "Number of B sig Q = ", numsigB_q
print "Number of YIB sig Q = ", numsigYCB_q
print "Number of Y sig BR = ", numsigY_br
print "Number of B sig BR = ", numsigB_br

                
    