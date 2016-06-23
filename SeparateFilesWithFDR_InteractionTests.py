# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 09:46:56 2016

@author: patrick
"""

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
INDIR="/Users/patrick/Documents/Research/EarlyLate/"
OUTDIR="/Users/patrick/Documents/Research/EarlyLate/"

INPUT_FILE="EL_Likelihoods_Interactions_20160602-0925.csv"

p_B_im_M12_cutoff = 4.37937548774e-06
p_B_im_M13_cutoff = 2.03748819271e-05
p_B_im_M23_cutoff = 4.41809928877e-06
p_B_q_M12_cutoff = 8.9285116539e-05
p_B_q_M13_cutoff = 0.000124947864633
p_B_q_M23_cutoff = 6.75887478276e-06
p_B_13_M12_cutoff = 3.59405881014e-05
p_B_13_M13_cutoff = 0.000127648970888
p_B_13_M23_cutoff = 2.88344135854e-05
p_B_14_M12_cutoff = 9.2059703091e-06
p_B_14_M13_cutoff = 1.22117519419e-05
p_B_14_M23_cutoff = 2.99643132326e-07

sigfile=open(OUTDIR+INPUT_FILE[:-4]+"_sig.csv","wb")
sigwriter=csv.writer(sigfile,delimiter=",",dialect='excel')    

nonsigfile=open(OUTDIR+INPUT_FILE[:-4]+"_nonsig.csv","wb")
nonsigwriter=csv.writer(nonsigfile,delimiter=",",dialect='excel')    

sites=0
numsigB_im_M12=0
numsigB_im_M13=0
numsigB_im_M23=0
numsigB_q_M12=0
numsigB_q_M13=0
numsigB_q_M23=0
numsigB_13_M12=0
numsigB_13_M13=0
numsigB_13_M23=0
numsigB_14_M12=0
numsigB_14_M13=0
numsigB_14_M23=0


with open(INDIR+INPUT_FILE,"rb") as sites_file:
    scaffs=[]
    snps=[]
    counts=[]
    for i,site in enumerate(sites_file):
        sites+=1
        if sites%10000==0:
            print sites
        sigB_M12=0
        sigB_M13=0
        sigB_M23=0
        
        site=site.strip("\n")
        site=site.strip("\r")
        site=site.split(",")
        if i == 0:      
            
            for k,j in enumerate(site):
                print k,j
            sigwriter.writerow(site+["sigB_M12","sigB_M13","sigB_M23"])
            nonsigwriter.writerow(site+["sigB_M12","sigB_M13","sigB_M23"])
        
        else:
            scaff=site[0]
            try:
                pos=float(site[1])
            except IndexError:
                pass
            try:
                p_B_M12=float(site[3])               
            except (ValueError, IndexError):
                p_B_M12=99
                sigB_M12="-"
            try:
                p_B_M13=float(site[5])               
            except (ValueError, IndexError):
                p_B_M13=99
                sigB_M13="-"
            try:
                p_B_M23=float(site[7])               
            except (ValueError, IndexError):
                p_B_M23=99
                sigB_M23="-"
            try:
                pop=site[8]             
            except (ValueError, IndexError):
                pass                           
            
            
            if pop=="IM":
                if p_B_M12 <= p_B_im_M12_cutoff:
                    sigB_M12=1
                    numsigB_im_M12+=1
                if p_B_M13 <= p_B_im_M13_cutoff:
                    sigB_M13=1
                    numsigB_im_M13+=1
                if p_B_M23 <= p_B_im_M23_cutoff:
                    sigB_M23=1
                    numsigB_im_M23+=1
            if pop=="Q":
                if p_B_M12 <= p_B_q_M12_cutoff:
                    sigB_M12=1
                    numsigB_q_M12+=1
                if p_B_M13 <= p_B_q_M13_cutoff:
                    sigB_M13=1
                    numsigB_q_M13+=1
                if p_B_M23 <= p_B_q_M23_cutoff:
                    sigB_M23=1
                    numsigB_q_M23+=1
            if pop=="2013":
                if p_B_M12 <= p_B_13_M12_cutoff:
                    sigB_M12=1
                    numsigB_13_M12+=1
                if p_B_M13 <= p_B_13_M13_cutoff:
                    sigB_M13=1
                    numsigB_13_M13+=1
                if p_B_M23 <= p_B_13_M23_cutoff:
                    sigB_M23=1
                    numsigB_13_M23+=1
            if pop=="2014":
                if p_B_M12 <= p_B_14_M12_cutoff:
                    sigB_M12=1
                    numsigB_14_M12+=1
                if p_B_M13 <= p_B_14_M13_cutoff:
                    sigB_M13=1
                    numsigB_14_M13+=1
                if p_B_M23 <= p_B_14_M23_cutoff:
                    sigB_M23=1
                    numsigB_14_M23+=1

            if any(k == 1 for k in [sigB_M12,sigB_M13,sigB_M23]):
                sigwriter.writerow(site+[sigB_M12,sigB_M13,sigB_M23])                
            else:
                nonsigwriter.writerow(site+[sigB_M12,sigB_M13,sigB_M23]) 
                
print "Number of B sig IM M12,M13,M23 = ", numsigB_im_M12,numsigB_im_M13,numsigB_im_M23
print "Number of B sig Q M12,M13,M23 = ", numsigB_q_M12,numsigB_q_M13,numsigB_q_M23
print "Number of B sig 2013 M12,M13,M23 = ", numsigB_13_M12,numsigB_13_M13,numsigB_13_M23
print "Number of B sig 2014 M12,M13,M23 = ", numsigB_14_M12,numsigB_14_M13,numsigB_14_M23

