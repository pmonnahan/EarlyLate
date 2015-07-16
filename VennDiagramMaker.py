# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 20:03:26 2014

@author: patrick
"""
from scipy.stats import *
import csv
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles

popT=20.0 #Number of reads required to analye a site for a particular bulk
bulkT=10.0

with open("/Users/patrick/Google Drive/Research/EarlyLate/2013 All sequence data/counts_EL_all.txt","rb") as sites_file:  
    BR_fixed_ref=0
    BR_fixed_alt=0
    BR_only_poly=0
    IM_fixed_ref=0
    IM_fixed_alt=0
    IM_only_poly=0
    Q_fixed_ref=0
    Q_fixed_alt=0
    Q_only_poly=0
    All_poly=0
    IMQ_poly=0
    IMBR_poly=0
    QBR_poly=0
            
    BR_only_count=0
    IM_only_count=0
    Q_only_count=0
    IMBR_count=0
    QBR_count=0
    IMQ_count=0
    All3_count=0 
        
    for i,site in enumerate(sites_file):
        site=site.strip("\n")
        site=site.split("\t")
        site=site[0:16]
        scaff,pos,ref,alt,bre_ref,bre_alt,brl_ref,brl_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe_ref,qe_alt,ql_ref,ql_alt=site  
        
        bre_ref=float(bre_ref)
        bre_alt=float(bre_alt)
        brl_ref=float(brl_ref)
        brl_alt=float(brl_alt)
        ime_ref=float(ime_ref)
        ime_alt=float(ime_alt)
        iml_ref=float(iml_ref)
        iml_alt=float(iml_alt)
        qe_ref=float(qe_ref)
        qe_alt=float(qe_alt)
        ql_ref=float(ql_ref)
        ql_alt=float(ql_alt)
        
        #Calculate read counts for bulks and populations
        bre_count=bre_ref+bre_alt
        brl_count=brl_ref+brl_alt
        ime_count=ime_ref+ime_alt
        iml_count=iml_ref+iml_alt
        qe_count=qe_ref+qe_alt
        ql_count=ql_ref+ql_alt
        br_count=bre_count+brl_count
        im_count=ime_count+iml_count
        q_count=qe_count+ql_count
        
        #Calculate allele frequencies
        bre_rf=None 
        brl_rf=None
        br_rfd=None
        im_rf=None
        ime_rf=None 
        iml_rf=None
        im_rfd=None
        q_rf=None
        qe_rf=None 
        ql_rf=None
        q_rfd=None
        br_rf=None
        im_rf=None
        q_rf=None
        if bre_count>0:
            try:
                bre_rf=bre_ref/bre_count
            except ZeroDivisionError:
                bre_rf=0.0
        if brl_count >0:
            try:
                brl_rf=brl_ref/brl_count
            except ZeroDivisionError:
                brl_rf=0.0
        if ime_count>0:
            try:
                ime_rf=ime_ref/ime_count
            except ZeroDivisionError:
                ime_rf=0.0
        if iml_count>0:
            try:
                iml_rf=iml_ref/iml_count
            except ZeroDivisionError:
                iml_rf=0.0
        if qe_count>0:
            try:
                qe_rf=qe_ref/qe_count
            except ZeroDivisionError:
                qe_rf=0.0
        if ql_count>0:
            try:
                ql_rf=ql_ref/ql_count
            except ZeroDivisionError:
                ql_rf=0.0
        if br_count>0:
            try:
                br_rf=(bre_ref+brl_ref)/br_count
            except ZeroDivisionError:
                br_rf=0.0
        if im_count>0:
            try:
                im_rf=(ime_ref+iml_ref)/im_count
            except ZeroDivisionError:
                im_rf=0.0
        if q_count>0:
            try:
                q_rf=(qe_ref+ql_ref)/q_count
            except ZeroDivisionError:
                q_rf=0.0
        
        #Determine if population passes filter
        if bre_count >= bulkT and brl_count >= bulkT:
            BR_pass=True
        else:
            BR_pass=False
        if ime_count >= bulkT and iml_count >= bulkT:
            IM_pass=True
        else:
            IM_pass=False
        if qe_count >= bulkT and ql_count >= bulkT:
            Q_pass=True
        else:
            Q_pass=False
            
        #Determine if population is polymorphic for the site
        if BR_pass==True and br_rf > 0.0 and br_rf < 1.0:
            BR_poly=True
        else:
            BR_poly=False
        if IM_pass==True and im_rf > 0.0 and im_rf < 1.0:
            IM_poly=True
        else:
            IM_poly=False
        if Q_pass==True and q_rf > 0.0 and q_rf < 1.0:
            Q_poly=True
        else:
            Q_poly=False
            
        #Calc venn_diagram stats.            
        if BR_pass==True and IM_pass==True and Q_pass==True:
            All3_count+=1
            if BR_poly==True and IM_poly==True and Q_poly==True:
                All_poly+=1
#                print im_rf,br_rf,q_rf,bre_ref+brl_ref,br_count,ime_ref+iml_ref,im_count,qe_ref+ql_ref,q_count
            elif BR_poly==False and IM_poly==True and Q_poly==True:
                IMQ_poly+=1
                if br_rf==0.0:
                    BR_fixed_alt+=1
                else:
                    BR_fixed_ref+=1
            elif BR_poly==True and IM_poly==False and Q_poly==True:
                QBR_poly+=1
                if im_rf==0.0:
                    IM_fixed_alt+=1
                else:
                    IM_fixed_ref+=1
            elif BR_poly==True and IM_poly==True and Q_poly==False:
                IMBR_poly+=1
                if q_rf==0.0:
                    Q_fixed_alt+=1
                else:
                    Q_fixed_ref+=1
            elif BR_poly==False and IM_poly==False and Q_poly==True:
                Q_only_poly+=1
                if im_rf==0.0:
                    IM_fixed_alt+=1
                else:
                    IM_fixed_ref+=1
                if br_rf==0.0:
                    BR_fixed_alt+=1
                else:
                    BR_fixed_ref+=1
            elif BR_poly==False and IM_poly==True and Q_poly==False:
                IM_only_poly+=1
                if br_rf==0.0:
                    BR_fixed_alt+=1
                else:
                    BR_fixed_ref+=1
                if q_rf==0.0:
                    Q_fixed_alt+=1
                else:
                    Q_fixed_ref+=1
            elif BR_poly==True and IM_poly==False and Q_poly==False:
                BR_only_poly+=1
                if im_rf==0.0:
                    IM_fixed_alt+=1
                else:
                    IM_fixed_ref+=1
                if q_rf==0.0:
                    Q_fixed_alt+=1
                else:
                    Q_fixed_ref+=1
                
        elif BR_pass==True and IM_pass==True and Q_pass==False:
            IMBR_count+=1
        elif BR_pass==True and IM_pass==False and Q_pass==True:
            QBR_count+=1
        elif BR_pass==False and IM_pass==True and Q_pass==True:
            IMQ_count+=1
        elif BR_pass==False and IM_pass==False and Q_pass==True:
            Q_only_count+=1
        elif BR_pass==False and IM_pass==True and Q_pass==False:
            IM_only_count+=1            
        elif BR_pass==True and IM_pass==False and Q_pass==False:
            BR_only_count+=1
        
        if im_rf==0.0 or br_rf==0.0 or q_rf==0.0:
            print im_rf,br_rf, q_rf
        
                    
print (BR_fixed_alt,BR_fixed_ref,IM_fixed_alt,IM_fixed_ref,Q_fixed_alt,Q_fixed_ref)

plt.figure(figsize=(4,4))
v = venn3(subsets=(BR_only_count,IM_only_count,IMBR_count,Q_only_count,QBR_count,IMQ_count,All3_count), set_labels = ('BR', 'IM', 'Q'))
v.get_patch_by_id('100').set_alpha(1.0)
plt.title("Sites with data")
plt.show()

plt.figure(figsize=(4,4))
v = venn3(subsets=(BR_only_poly,IM_only_poly,IMBR_poly,Q_only_poly,QBR_poly,IMQ_poly,All_poly), set_labels = ('BR', 'IM', 'Q'))
v.get_patch_by_id('100').set_alpha(1.0)
plt.title("Polymorphism for sites")
plt.show()

#BR_only=list(csv.reader(open("/Users/patrick/Google Drive/Research/EarlyLate/2013 All sequence data/EL_BR_only.csv")))
#IM_only=list(csv.reader(open("/Users/patrick/Google Drive/Research/EarlyLate/2013 All sequence data/EL_IM_only.csv")))
#Q_only=list(csv.reader(open("/Users/patrick/Google Drive/Research/EarlyLate/2013 All sequence data/EL_Q_only.csv")))
#IM_BR=list(csv.reader(open("/Users/patrick/Google Drive/Research/EarlyLate/2013 All sequence data/EL_IMBR.csv")))
#Q_BR=list(csv.reader(open("/Users/patrick/Google Drive/Research/EarlyLate/2013 All sequence data/EL_QBR.csv")))
#IM_Q=list(csv.reader(open("/Users/patrick/Google Drive/Research/EarlyLate/2013 All sequence data/EL_IMQ.csv")))
#All3=list(csv.reader(open("/Users/patrick/Google Drive/Research/EarlyLate/2013 All sequence data/EL_Fst_all3.csv")))
#
#BR_only_count=len(BR_only)
#IM_only_count=len(IM_only)
#Q_only_count=len(Q_only)
#IMBR_count=len(IM_BR)
#QBR_count=len(Q_BR)
#IMQ_count=len(IM_Q)
#All3_count=len(All3)

#        if br_count > popT:
#            br_rf=(bre_ref+brl_ref)/(bre_ref+bre_alt+brl_ref+brl_alt)
#            if bre_count> bulkT and brl_count > bulkT:
#                bre_rf=bre_ref/(bre_ref+bre_alt) #Browder Ridge Early_Reference Frequency
#                brl_rf=brl_ref/(brl_ref+brl_alt)
#                br_rfd=(bre_rf-brl_rf)
#                if (bre_rf + brl_rf)==0.0:
#                    BR_fixed_alt+=1
#                elif(bre_rf + brl_rf)==2.0:
#                    BR_fixed_ref+=1
#                else:
#                    BR_poly+=1
#        if im_count > popT:
#            im_rf=(ime_ref+iml_ref)/(ime_ref+ime_alt+iml_ref+iml_alt)
#            if ime_count> bulkT and iml_count > bulkT:
#                ime_rf=ime_ref/(ime_ref+ime_alt) #Browder Ridge Early_Reference Frequency
#                iml_rf=iml_ref/(iml_ref+iml_alt)
#                im_rfd=(ime_rf-iml_rf)
#                if (ime_rf + iml_rf)==0.0:
#                    IM_fixed_alt+=1
#                elif(ime_rf + iml_rf)==2.0:
#                    IM_fixed_ref+=1
#                else:
#                    IM_poly+=1
#        if q_count > popT:
#            q_rf=(qe_ref+ql_ref)/(qe_ref+qe_alt+ql_ref+ql_alt)
#            if qe_count> bulkT and ql_count > bulkT:
#                qe_rf=qe_ref/(qe_ref+qe_alt) #Browder Ridge Early_Reference Frequency
#                ql_rf=ql_ref/(ql_ref+ql_alt)
#                q_rfd=(qe_rf-ql_rf)
#                if (qe_rf + ql_rf)==0.0:
#                    Q_fixed_alt+=1
#                elif(qe_rf + ql_rf)==2.0:
#                    Q_fixed_ref+=1
#                else:
#                    Q_poly+=1