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
INPUT_FILE="EL_Likelihoods_20160404-1106.csv"

p_Y_im_cutoff = 2.84717126029e-05
p_YCB_im_cutoff = 4.41809928877e-06
p_B_im13_cutoff = 1.8735352155e-05 
p_B_br13_cutoff = 3.46567027194e-07
p_B_q14_cutoff = 1.8735352155e-05 
p_B_im14_cutoff = 3.46567027194e-07
p_Y_q_cutoff = 9.4924509668e-05
p_B_q13_cutoff = 0.000222275904042
p_YCB_q_cutoff = 6.75887478276e-06
p_PCB_13_cutoff=2.88344135854e-05
p_PCB_14_cutoff=0.0

sigfile=open(OUTDIR+INPUT_FILE[:-4]+"_sig.csv","wb")
sigwriter=csv.writer(sigfile,delimiter=",",dialect='excel')    

nonsigfile=open(OUTDIR+INPUT_FILE[:-4]+"_nonsig.csv","wb")
nonsigwriter=csv.writer(nonsigfile,delimiter=",",dialect='excel')    

sites=0
numsigY_im=0
numsigB_im13=0
numsigYCB_im=0
numsigY_q=0
numsigB_q13=0
numsigB_q14=0
numsigB_im14=0
numsigYCB_q=0
numsigY_br=0
numsigB_br13=0
numsigPCB_13=0
numsigPCB_14=0

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
        sigB_im13=0
        sigB_im14=0
        sigY_q=0
        sigYCB_q=0
        sigB_q13=0
        sigB_q14=0
        sigY_br=0
        sigB_br13=0
        sigPCB_13=0
        sigPCB_14=0
        site=site.strip("\n")
        site=site.strip("\r")
        site=site.split(",")
        if i == 0:      
            
            for k,j in enumerate(site):
                print k,j
            sigwriter.writerow(site+["sigYim","sigBim13","sigYCBim","sigYq","sigBq13","sigYCBq","sigYbr","sigBbr13","sigBim14","sigBq14","sigPCB_13","sigPCB_14"])
            nonsigwriter.writerow(site+["sigYim","sigBim13","sigYCBim","sigYq","sigBq13","sigYCBq","sigYbr","sigBbr13","sigBbr13","sigBim14","sigBq14","sigPCB_13","sigPCB_14"])
        
        else:
            scaff=site[0]
            try:
                pos=float(site[1])
            except IndexError:
                pass
            try:
                p_Y_im=float(site[60])                
            except (ValueError, IndexError):
                p_Y_im=99
                sigY_im="-" 
            try:
                p_B_im13=float(site[27])               
            except (ValueError, IndexError):
                p_B_im13=99
                sigB_im13="-"
            try:
                p_B_im14=float(site[31])               
            except (ValueError, IndexError):
                p_B_im14=99
                sigB_im14="-"
            try:
                p_YCB_im=float(site[35])              
            except (ValueError, IndexError):
                p_YCB_im = 99
                sigYCB_im="-"
            try:
                p_Y_q=float(site[60])                
            except (ValueError, IndexError):
                p_Y_q=99     
                sigY_q="-"
            try:
                p_B_q13=float(site[29])               
            except (ValueError, IndexError):
                p_B_q13=99
                sigB_q13="-"
            try:
                p_B_q14=float(site[33])               
            except (ValueError, IndexError):
                p_B_q14=99
                sigB_q14="-"
            try:
                p_YCB_q=float(site[37])              
            except (ValueError, IndexError):
                p_YCB_q = 99
                sigYCB_q="-"
                              
            try:
                p_B_br13=float(site[25])               
            except (ValueError, IndexError):
                p_B_br13=99
                sigB_br13="-"
            try:
                p_PCB_13=float(site[39])              
            except (ValueError, IndexError):
                p_PCB_13 = 99
                sigPCB_13="-"
            try:
                p_PCB_14=float(site[41])              
            except (ValueError, IndexError):
                p_PCB_14 = 99
                sigPCB_14="-"
            
            
            if p_Y_im <= p_Y_im_cutoff:
                sigY_im=1
                numsigY_im+=1
            if p_B_im13 <= p_B_im13_cutoff:
                sigB_im13=1
                numsigB_im13+=1
            if p_B_im14 <= p_B_im14_cutoff:
                sigB_im14=1
                numsigB_im14+=1
            if p_YCB_im <= p_YCB_im_cutoff:
                sigYCB_im=1
                numsigYCB_im+=1
            if p_Y_q <= p_Y_q_cutoff:
                sigY_q=1
                numsigY_q+=1
            if p_B_q13 <= p_B_q13_cutoff:
                sigB_q13=1
                numsigB_q13+=1
            if p_B_q14 <= p_B_q14_cutoff:
                sigB_q14=1
                numsigB_q14+=1
            if p_YCB_q <= p_YCB_q_cutoff:
                sigYCB_q=1
                numsigYCB_q+=1
            if p_B_br13 <= p_B_br13_cutoff:
                sigB_br13=1
                numsigB_br13+=1
            if p_PCB_13 <= p_PCB_13_cutoff:
                sigPCB_13=1
                numsigPCB_13+=1
            if p_PCB_14 <= p_PCB_14_cutoff:
                sigPCB_14=1
                numsigPCB_14+=1

            if any(k == 1 for k in [sigY_im,sigB_im13,sigB_im14,sigYCB_im,sigY_q,sigB_q13,sigB_q14,sigYCB_q,sigB_br13,sigPCB_13,sigPCB_14]):
                sigwriter.writerow(site+[sigY_im,sigB_im13,sigYCB_im,sigY_q,sigB_q13,sigYCB_q,sigY_br,sigB_br13,sigB_im14,sigB_q14,sigPCB_13,sigPCB_14])                
            else:
                nonsigwriter.writerow(site+[sigY_im,sigB_im13,sigYCB_im,sigY_q,sigB_q13,sigYCB_q,sigY_br,sigB_br13,sigB_im14,sigB_q14,sigPCB_13,sigPCB_14]) 
                
print "Number of Y sig IM = ", numsigY_im
print "Number of B sig IM13 = ", numsigB_im13
print "Number of B sig IM14 = ", numsigB_im14
print "Number of YIB sig IM = ", numsigYCB_im
print "Number of Y sig Q = ", numsigY_q
print "Number of B sig Q13 = ", numsigB_q13
print "Number of B sig Q14 = ", numsigB_q14
print "Number of YIB sig Q = ", numsigYCB_q
print "Number of PIB sig 13 = ",numsigPCB_13
print "Number of PIB sig 14 = ",numsigPCB_14
print "Number of B sig BR = ", numsigB_br13

                
    