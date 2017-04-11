# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 20:03:26 2014
1/5/15 should be calculating allele frequencies across years correctly.  Reference frequency should be average of frequency in each year, instead of pooled read counts across years divided by total read counts.
        Also, removed BR_14 data from entering calculations of Z_cum.
@author: patrick
"""

from scipy.stats import *
import csv
import math

popT=20.0 #Number of reads required to analye a site for a particular bulk
bulkT=10.0
BS_BR13=0.06
BS_IM13=0.045
BS_Q13=0.045
BS_BR14=0.05
BS_IM14=0.045
BS_Q14=0.045
BS_BR=0.03
BS_IM=0.03
BS_Q=0.03

bulkSize=100

Zs_q13=[]
Zs_br13=[]
Zs_im13=[]

Zs_q14=[]
Zs_br14=[]
Zs_im14=[]

Zs_q=[]
Zs_br=[]
Zs_im=[]

#BR13_only=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_BR13_only.csv","wb")
#IM13_only=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_IM13_only.csv","wb")
#Q13_only=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_Q13_only.csv","wb")
#IM_BR13=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_IMBR13.csv","wb")
#Q_BR13=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_QBR13.csv","wb")
#IM_Q13=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_IMQ13.csv","wb")
#newcsv13=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_Fst_all3_13.csv","wb")
#
#BR14_only=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_BR14_only.csv","wb")
#IM14_only=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_IM14_only.csv","wb")
#Q14_only=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_Q14_only.csv","wb")
#IM_BR14=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_IMBR14.csv","wb")
#Q_BR14=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_QBR14.csv","wb")
#IM_Q14=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_IMQ14.csv","wb")
#newcsv14=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_Fst_all3_14.csv","wb")

BR_only=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_BR_only_new.csv","wb")
IM_only=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_IM_only_new.csv","wb")
Q_only=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_Q_only_new.csv","wb")
IM_BR=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_IMBR_new.csv","wb")
Q_BR=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_QBR_new.csv","wb")
IM_Q=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_IMQ_new.csv","wb")
newcsv=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_Fst_all3_bothyears_new.csv","wb")




#allwriter13=csv.writer(newcsv13,delimiter=",",dialect='excel')
#allwriter14=csv.writer(newcsv14,delimiter=",",dialect='excel')
#BR13writer=csv.writer(BR13_only,delimiter=",",dialect='excel')
#IM13writer=csv.writer(IM13_only,delimiter=",",dialect='excel')
#Q13writer=csv.writer(Q13_only,delimiter=",",dialect='excel')
#IMBR13writer=csv.writer(IM_BR13,delimiter=",",dialect='excel')
#QBR13writer=csv.writer(Q_BR13,delimiter=",",dialect='excel')
#IMQ13writer=csv.writer(IM_Q13,delimiter=",",dialect='excel')
#
#BR14writer=csv.writer(BR14_only,delimiter=",",dialect='excel')
#IM14writer=csv.writer(IM14_only,delimiter=",",dialect='excel')
#Q14writer=csv.writer(Q14_only,delimiter=",",dialect='excel')
#IMBR14writer=csv.writer(IM_BR14,delimiter=",",dialect='excel')
#QBR14writer=csv.writer(Q_BR14,delimiter=",",dialect='excel')
#IMQ14writer=csv.writer(IM_Q14,delimiter=",",dialect='excel')

allwriter=csv.writer(newcsv,delimiter=",",dialect='excel')
BRwriter=csv.writer(BR_only,delimiter=",",dialect='excel')
IMwriter=csv.writer(IM_only,delimiter=",",dialect='excel')
Qwriter=csv.writer(Q_only,delimiter=",",dialect='excel')
IMBRwriter=csv.writer(IM_BR,delimiter=",",dialect='excel')
QBRwriter=csv.writer(Q_BR,delimiter=",",dialect='excel')
IMQwriter=csv.writer(IM_Q,delimiter=",",dialect='excel')

#BR13writer.writerow(["scaff","pos","bre_ref","bre_alt","brl_ref","brl_alt","Z_br","p_br"])
#IM13writer.writerow(["scaff","pos","ime_ref","ime_alt","iml_ref","iml_alt","Z_im","p_im"])
#Q13writer.writerow(["scaff","pos","qe_ref","qe_alt","ql_ref","ql_alt","Z_q","p_q"])
#IMBR13writer.writerow(["scaff","pos","bre_ref","bre_alt","brl_ref","brl_alt","ime_ref","ime_alt","iml_ref","iml_alt","Z_br","p_br","Z_im","p_im","Z_cum","p_cum","Fst"])
#IMQ13writer.writerow(["scaff","pos","ime_ref","ime_alt","iml_ref","iml_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_im","p_im","Z_q","p_q","Z_cum","p_cum","Fst"])
#QBR13writer.writerow(["scaff","pos","bre_ref","bre_alt","brl_ref","brl_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br","p_br","Z_q","p_q","Z_cum","p_cum","Fst"])
#allwriter13.writerow(["scaff","pos","bre13_ref","bre13_alt","brl13_ref","brl13_alt","bre14_ref","bre14_alt","ime13_ref","ime13_alt","iml13_ref","iml13_alt","ime14_ref","ime14_alt","iml14_ref","iml14_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe13_ref","qe13_alt","ql13_ref","ql13_alt","qe14_ref","qe14_alt","ql14_ref","ql14_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br13","p_br13","Z_br14","p_br14","Z_br","p_br","Z_im13","p_im13","Z_im14","p_im14","Z_im","p_im","Z_q13","p_q13","Z_q14","p_q14","Z_q","p_q","pop_num","Z_cum","p_Z_cum","Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","q_rfd_bnYears","qe_rfd_bnYears","ql_rfd_bnYears","im_rfd_bnYears","ime_rfd_bnYears","iml_rfd_bnYears","br_rfd_bnYears","bre_rfd_bnYears"])
#allwriter14.writerow(["scaff","pos","bre13_ref","bre13_alt","brl13_ref","brl13_alt","bre14_ref","bre14_alt","ime13_ref","ime13_alt","iml13_ref","iml13_alt","ime14_ref","ime14_alt","iml14_ref","iml14_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe13_ref","qe13_alt","ql13_ref","ql13_alt","qe14_ref","qe14_alt","ql14_ref","ql14_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br13","p_br13","Z_br14","p_br14","Z_br","p_br","Z_im13","p_im13","Z_im14","p_im14","Z_im","p_im","Z_q13","p_q13","Z_q14","p_q14","Z_q","p_q","pop_num","Z_cum","p_Z_cum","Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","q_rfd_bnYears","qe_rfd_bnYears","ql_rfd_bnYears","im_rfd_bnYears","ime_rfd_bnYears","iml_rfd_bnYears","br_rfd_bnYears","bre_rfd_bnYears"])
#
#BR14writer.writerow(["scaff","pos","bre_ref","bre_alt","brl_ref","brl_alt","Z_br","p_br"])
#IM14writer.writerow(["scaff","pos","ime_ref","ime_alt","iml_ref","iml_alt","Z_im","p_im"])
#Q14writer.writerow(["scaff","pos","qe_ref","qe_alt","ql_ref","ql_alt","Z_q","p_q"])
#IMBR14writer.writerow(["scaff","pos","bre_ref","bre_alt","brl_ref","brl_alt","ime_ref","ime_alt","iml_ref","iml_alt","Z_br","p_br","Z_im","p_im","Z_cum","p_cum","Fst"])
#IMQ14writer.writerow(["scaff","pos","ime_ref","ime_alt","iml_ref","iml_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_im","p_im","Z_q","p_q","Z_cum","p_cum","Fst"])
#QBR14writer.writerow(["scaff","pos","bre_ref","bre_alt","brl_ref","brl_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br","p_br","Z_q","p_q","Z_cum","p_cum","Fst"])

allwriter.writerow(["scaff","pos","bre13_ref","bre13_alt","brl13_ref","brl13_alt","bre14_ref","bre14_alt","ime13_ref","ime13_alt","iml13_ref","iml13_alt","ime14_ref","ime14_alt","iml14_ref","iml14_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe13_ref","qe13_alt","ql13_ref","ql13_alt","qe14_ref","qe14_alt","ql14_ref","ql14_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br13","p_br13","Z_br14","p_br14","Z_br","p_br","Z_im13","p_im13","Z_im14","p_im14","Z_im","p_im","Z_q13","p_q13","Z_q14","p_q14","Z_q","p_q","pop_num","Z_cum","p_Z_cum","Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","q_rfd_bnYears","qe_rfd_bnYears","ql_rfd_bnYears","im_rfd_bnYears","ime_rfd_bnYears","iml_rfd_bnYears","br_rfd_bnYears","bre_rfd_bnYears"])
BRwriter.writerow(["scaff","pos","bre13_ref","bre13_alt","brl13_ref","brl13_alt","bre14_ref","bre14_alt","ime13_ref","ime13_alt","iml13_ref","iml13_alt","ime14_ref","ime14_alt","iml14_ref","iml14_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe13_ref","qe13_alt","ql13_ref","ql13_alt","qe14_ref","qe14_alt","ql14_ref","ql14_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br13","p_br13","Z_br14","p_br14","Z_br","p_br","Z_im13","p_im13","Z_im14","p_im14","Z_im","p_im","Z_q13","p_q13","Z_q14","p_q14","Z_q","p_q","pop_num","Z_cum","p_Z_cum","Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","q_rfd_bnYears","qe_rfd_bnYears","ql_rfd_bnYears","im_rfd_bnYears","ime_rfd_bnYears","iml_rfd_bnYears","br_rfd_bnYears","bre_rfd_bnYears"])
IMwriter.writerow(["scaff","pos","bre13_ref","bre13_alt","brl13_ref","brl13_alt","bre14_ref","bre14_alt","ime13_ref","ime13_alt","iml13_ref","iml13_alt","ime14_ref","ime14_alt","iml14_ref","iml14_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe13_ref","qe13_alt","ql13_ref","ql13_alt","qe14_ref","qe14_alt","ql14_ref","ql14_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br13","p_br13","Z_br14","p_br14","Z_br","p_br","Z_im13","p_im13","Z_im14","p_im14","Z_im","p_im","Z_q13","p_q13","Z_q14","p_q14","Z_q","p_q","pop_num","Z_cum","p_Z_cum","Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","q_rfd_bnYears","qe_rfd_bnYears","ql_rfd_bnYears","im_rfd_bnYears","ime_rfd_bnYears","iml_rfd_bnYears","br_rfd_bnYears","bre_rfd_bnYears"])
Qwriter.writerow(["scaff","pos","bre13_ref","bre13_alt","brl13_ref","brl13_alt","bre14_ref","bre14_alt","ime13_ref","ime13_alt","iml13_ref","iml13_alt","ime14_ref","ime14_alt","iml14_ref","iml14_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe13_ref","qe13_alt","ql13_ref","ql13_alt","qe14_ref","qe14_alt","ql14_ref","ql14_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br13","p_br13","Z_br14","p_br14","Z_br","p_br","Z_im13","p_im13","Z_im14","p_im14","Z_im","p_im","Z_q13","p_q13","Z_q14","p_q14","Z_q","p_q","pop_num","Z_cum","p_Z_cum","Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","q_rfd_bnYears","qe_rfd_bnYears","ql_rfd_bnYears","im_rfd_bnYears","ime_rfd_bnYears","iml_rfd_bnYears","br_rfd_bnYears","bre_rfd_bnYears"])
IMBRwriter.writerow(["scaff","pos","bre13_ref","bre13_alt","brl13_ref","brl13_alt","bre14_ref","bre14_alt","ime13_ref","ime13_alt","iml13_ref","iml13_alt","ime14_ref","ime14_alt","iml14_ref","iml14_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe13_ref","qe13_alt","ql13_ref","ql13_alt","qe14_ref","qe14_alt","ql14_ref","ql14_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br13","p_br13","Z_br14","p_br14","Z_br","p_br","Z_im13","p_im13","Z_im14","p_im14","Z_im","p_im","Z_q13","p_q13","Z_q14","p_q14","Z_q","p_q","pop_num","Z_cum","p_Z_cum","Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","q_rfd_bnYears","qe_rfd_bnYears","ql_rfd_bnYears","im_rfd_bnYears","ime_rfd_bnYears","iml_rfd_bnYears","br_rfd_bnYears","bre_rfd_bnYears"])
IMQwriter.writerow(["scaff","pos","bre13_ref","bre13_alt","brl13_ref","brl13_alt","bre14_ref","bre14_alt","ime13_ref","ime13_alt","iml13_ref","iml13_alt","ime14_ref","ime14_alt","iml14_ref","iml14_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe13_ref","qe13_alt","ql13_ref","ql13_alt","qe14_ref","qe14_alt","ql14_ref","ql14_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br13","p_br13","Z_br14","p_br14","Z_br","p_br","Z_im13","p_im13","Z_im14","p_im14","Z_im","p_im","Z_q13","p_q13","Z_q14","p_q14","Z_q","p_q","pop_num","Z_cum","p_Z_cum","Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","q_rfd_bnYears","qe_rfd_bnYears","ql_rfd_bnYears","im_rfd_bnYears","ime_rfd_bnYears","iml_rfd_bnYears","br_rfd_bnYears","bre_rfd_bnYears"])
QBRwriter.writerow(["scaff","pos","bre13_ref","bre13_alt","brl13_ref","brl13_alt","bre14_ref","bre14_alt","ime13_ref","ime13_alt","iml13_ref","iml13_alt","ime14_ref","ime14_alt","iml14_ref","iml14_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe13_ref","qe13_alt","ql13_ref","ql13_alt","qe14_ref","qe14_alt","ql14_ref","ql14_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br13","p_br13","Z_br14","p_br14","Z_br","p_br","Z_im13","p_im13","Z_im14","p_im14","Z_im","p_im","Z_q13","p_q13","Z_q14","p_q14","Z_q","p_q","pop_num","Z_cum","p_Z_cum","Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","q_rfd_bnYears","qe_rfd_bnYears","ql_rfd_bnYears","im_rfd_bnYears","ime_rfd_bnYears","iml_rfd_bnYears","br_rfd_bnYears","bre_rfd_bnYears"])

with open("/Users/patrick/Google Drive/Research/EarlyLate/counts_test.txt","rb") as sites_file:  
    for i,site in enumerate(sites_file):
        if i==0:
            print site
        if i>0:        
            site=site.strip("\n")
            site=site.split("\t")
            #site=site[0:16]
            scaff,pos,ref,alt,bre14_ref,bre14_alt,bre13_ref,bre13_alt,brl13_ref,brl13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,bre_ref,bre_alt,brl_ref,brl_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe_ref,qe_alt,ql_ref,ql_alt=site  
            
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
            
            bre13_ref=float(bre13_ref)
            bre13_alt=float(bre13_alt)
            brl13_ref=float(brl13_ref)
            brl13_alt=float(brl13_alt)
            ime13_ref=float(ime13_ref)
            ime13_alt=float(ime13_alt)
            iml13_ref=float(iml13_ref)
            iml13_alt=float(iml13_alt)
            qe13_ref=float(qe13_ref)
            qe13_alt=float(qe13_alt)
            ql13_ref=float(ql13_ref)
            ql13_alt=float(ql13_alt)
            
            bre14_ref=float(bre14_ref)
            bre14_alt=float(bre14_alt)
            ime14_ref=float(ime14_ref)
            ime14_alt=float(ime14_alt)
            iml14_ref=float(iml14_ref)
            iml14_alt=float(iml14_alt)
            qe14_ref=float(qe14_ref)
            qe14_alt=float(qe14_alt)
            ql14_ref=float(ql14_ref)
            ql14_alt=float(ql14_alt)
            
            bre13_count=bre13_ref+bre13_alt
            brl13_count=brl13_ref+brl13_alt
            ime13_count=ime13_ref+ime13_alt
            iml13_count=iml13_ref+iml13_alt
            qe13_count=qe13_ref+qe13_alt
            ql13_count=ql13_ref+ql13_alt
            br13_count=bre13_count+brl13_count
            im13_count=ime13_count+iml13_count
            q13_count=qe13_count+ql13_count
            
            bre14_count=bre14_ref+bre14_alt
            ime14_count=ime14_ref+ime14_alt
            iml14_count=iml14_ref+iml14_alt
            qe14_count=qe14_ref+qe14_alt
            ql14_count=ql14_ref+ql14_alt
            br14_count=bre14_count
            im14_count=ime14_count+iml14_count
            q14_count=qe14_count+ql14_count
            
            bre_count=bre_ref+bre_alt
            brl_count=brl_ref+brl_alt
            ime_count=ime_ref+ime_alt
            iml_count=iml_ref+iml_alt
            qe_count=qe_ref+qe_alt
            ql_count=ql_ref+ql_alt
            br_count=bre_count+brl_count
            im_count=ime_count+iml_count
            q_count=qe_count+ql_count
            
            br13_rf=None            
            bre13_rf=None 
            brl13_rf=None
            br13_rfd=None
            im13_rf=None
            ime13_rf=None 
            iml13_rf=None
            im13_rfd=None
            q13_rf=None
            qe13_rf=None 
            ql13_rf=None
            q13_rfd=None
            
            br14_rf=None            
            bre14_rf=None 
            brl14_rf=None
            br14_rfd=None
            im14_rf=None
            ime14_rf=None 
            iml14_rf=None
            im14_rfd=None
            q14_rf=None
            qe14_rf=None 
            ql14_rf=None
            q14_rfd=None
            
            br_rf=None            
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
            ####
            Fst_all_samples=None
            Fst_all=None
            Fst_pooled_pops_within_yr=None
            Fst_pooled_pops_across_yr=None
            
            Fst_IM_BR13=None
            Fst_IM_Q13=None
            Fst_BR_Q13=None
            Fst_IM_Q14=None
            Fst_IM_BR=None
            Fst_IM_Q=None
            Fst_BR_Q=None
            
            ssfd=None
            
            Z_br13=None
            Z_im13=None
            Z_q13=None
            Z_br14=None
            Z_im14=None
            Z_q14=None
            Z_br=None
            Z_im=None
            Z_q=None
            
            at_bre_rf=None
            at_brl_rf=None
            var_at_br_rfd=None
            at_ime_rf=None
            at_iml_rf=None
            var_at_im_rfd=None        
            at_qe_rf=None
            at_ql_rf=None
            var_at_q_rfd=None
            p_br=None
            p_im=None
            p_q=None
            
            at_bre13_rf=None
            at_brl13_rf=None
            var_at_br13_rfd=None
            at_ime13_rf=None
            at_iml13_rf=None
            var_at_im13_rfd=None        
            at_qe13_rf=None
            at_ql13_rf=None
            var_at_q13_rfd=None
            p_br13=None
            p_im13=None
            p_q13=None
            
            at_bre14_rf=None
            at_brl14_rf=None
            var_at_br14_rfd=None
            at_ime14_rf=None
            at_iml14_rf=None
            var_at_im14_rfd=None        
            at_qe14_rf=None
            at_ql14_rf=None
            var_at_q14_rfd=None
            p_br14=None
            p_im14=None
            p_q14=None
            
            pop_num=None
            Z_cum=None
            p_Z_cum=None
            
            q_rfd_bnYears=None
            qe_rfd_bnYears=None
            ql_rfd_bnYears=None
            im_rfd_bnYears=None
            ime_rfd_bnYears=None
            iml_rfd_bnYears=None
            br_rfd_bnYears=None
            bre_rfd_bnYears=None
                    
            
            if bre13_count> bulkT and brl13_count > bulkT:
                br13_rf=(bre13_ref+brl13_ref)/(bre13_ref+bre13_alt+brl13_ref+brl13_alt)                
                bre13_rf=bre13_ref/(bre13_ref+bre13_alt) #Browder Ridge Early_Reference Frequency
                brl13_rf=brl13_ref/(brl13_ref+brl13_alt)
                br13_rfd=(bre13_rf-brl13_rf)
                if (bre13_rf + brl13_rf)==0.0 or (bre13_rf + brl13_rf)==2.0:
                    Z_br13=0.0
                else:
                    at_brl13_rf=math.asin(brl13_rf**0.5)*2
                    at_bre13_rf=math.asin(bre13_rf**0.5)*2 #Angular transformation of bre_rf            
                    var_at_br13_rfd=(1/(bulkSize))+(1/bre13_count)+(1/brl13_count)+BS_BR13
                    Z_br13=(at_bre13_rf - at_brl13_rf)/(var_at_br13_rfd**0.5)
                    Zs_br13.append(Z_br13)
                    if Z_br13 < 0.0:
                        p_br13=norm.cdf(Z_br13)*2
                    else:
                        p_br13=(1.0-norm.cdf(Z_br13))*2
                            
            if bre14_count> bulkT:
                br14_rf=(bre14_ref)/(bre14_ref+bre14_alt)                
                bre14_rf=bre14_ref/(bre14_ref+bre14_alt) #Browder Ridge Early_Reference Frequency
                if (bre14_rf)==0.0 or (bre14_rf)==1.0:
                    Z_br14=0.0
                else:
                    at_bre14_rf=math.asin(bre14_rf**0.5)*2 #Angular transformation of bre_rf            
                    var_at_br14_rfd=(1/(bulkSize))+(1/bre14_count)+BS_BR14
                    Z_br14=(at_bre14_rf)/(var_at_br14_rfd**0.5)
                    Zs_br14.append(Z_br14)
                    if Z_br14 < 0.0:
                        p_br14=norm.cdf(Z_br14)*2
                    else:
                        p_br14=(1.0-norm.cdf(Z_br14))*2
                            
            if bre13_count > bulkT and bre14_count > bulkT and brl13_count > bulkT:
                br_rf=(br14_rf+br13_rf)/2                
                bre_rf=(bre14_rf+bre13_rf)/2
                brl_rf=brl13_rf
                br_rfd=(bre_rf-brl_rf)
                br_rfd_bnYears= br14_rf - br13_rf
                bre_rfd_bnYears= bre14_rf - bre13_rf
                if (bre_rf + brl_rf)==0.0 or (bre_rf + brl_rf)==2.0:
                    Z_br=0.0
                else:
                    at_brl_rf=math.asin(brl_rf**0.5)*2
                    at_bre_rf=math.asin(bre_rf**0.5)*2 #Angular transformation of bre_rf            
                    var_at_br_rfd=(1/(bulkSize))+(1/bre13_count)+(1/brl13_count)+(1/bre14_count)+BS_BR
                    Z_br=(at_bre_rf - at_brl_rf)/(var_at_br_rfd**0.5)
                    Zs_br.append(Z_br)
                    if Z_br < 0.0:
                        p_br=norm.cdf(Z_br)*2
                    else:
                        p_br=(1.0-norm.cdf(Z_br))*2
                
            if ime13_count> bulkT and iml13_count > bulkT:
                im13_rf=((ime13_ref)+(iml13_ref))/((ime13_ref)+(ime13_alt)+(iml13_ref)+(iml13_alt))                
                ime13_rf=(ime13_ref)/((ime13_ref)+(ime13_alt))
                iml13_rf=(iml13_ref)/((iml13_ref)+(iml13_alt))
                im13_rfd=(ime13_rf-iml13_rf)
                if (ime13_rf + iml13_rf)==0.0 or (ime13_rf + iml13_rf)==2.0:
                    Z_im13=0.0
                else:
                    at_iml13_rf=math.asin(iml13_rf**0.5)*2                    
                    at_ime13_rf=math.asin(ime13_rf**0.5)*2
                    var_at_im13_rfd=(1/(bulkSize))+(1/ime13_count)+(1/iml13_count)+BS_IM13
                    Z_im13=(at_ime13_rf - at_iml13_rf)/(var_at_im13_rfd**0.5)
                    Zs_im13.append(Z_im13)
                    if Z_im13 < 0.0:
                        p_im13=norm.cdf(Z_im13)*2
                    else:
                        p_im13=(1.0-norm.cdf(Z_im13))*2   
                            
            if ime14_count> bulkT and iml14_count > bulkT:
                im14_rf=((ime14_ref)+(iml14_ref))/((ime14_ref)+(ime14_alt)+(iml14_ref)+(iml14_alt))                
                ime14_rf=(ime14_ref)/((ime14_ref)+(ime14_alt))
                iml14_rf=(iml14_ref)/((iml14_ref)+(iml14_alt))
                im14_rfd=(ime14_rf-iml14_rf)
                if (ime14_rf + iml14_rf)==0.0 or (ime14_rf + iml14_rf)==2.0:
                    Z_im14=0.0
                else:
                    at_iml14_rf=math.asin(iml14_rf**0.5)*2                    
                    at_ime14_rf=math.asin(ime14_rf**0.5)*2
                    var_at_im14_rfd=(1/(bulkSize))+(1/ime14_count)+(1/iml14_count)+BS_IM14
                    Z_im14=(at_ime14_rf - at_iml14_rf)/(var_at_im14_rfd**0.5)
                    Zs_im14.append(Z_im14)
                    if Z_im14 < 0.0:
                        p_im14=norm.cdf(Z_im14)*2
                    else:
                        p_im14=(1.0-norm.cdf(Z_im14))*2 

            if ime13_count > bulkT and ime14_count > bulkT and iml13_count > bulkT and iml14_count > bulkT:
                im_rf=(im14_rf+im13_rf)/2                
                ime_rf=(ime14_rf+ime13_rf)/2
                iml_rf=(iml14_rf+iml13_rf)/2
                im_rfd=(ime_rf-iml_rf)
                im_rfd_bnYears= im14_rf - im13_rf
                ime_rfd_bnYears= ime14_rf - ime13_rf
                iml_rfd_bnYears= iml14_rf - iml13_rf
                if (ime_rf + iml_rf)==0.0 or (ime_rf + iml_rf)==2.0:
                    Z_im=0.0
                else:
                    at_iml_rf=math.asin(iml_rf**0.5)*2                    
                    at_ime_rf=math.asin(ime_rf**0.5)*2
                    var_at_im_rfd=(1/(bulkSize))+(1/ime13_count)+(1/iml13_count)+(1/ime14_count)+(1/iml14_count)+BS_IM
                    Z_im=(at_ime_rf - at_iml_rf)/(var_at_im_rfd**0.5)
                    Zs_im.append(Z_im)
                    if Z_im < 0.0:
                        p_im=norm.cdf(Z_im)*2
                    else:
                        p_im=(1.0-norm.cdf(Z_im))*2                             
                        
            if qe13_count> bulkT and ql13_count > bulkT:
                q13_rf=((qe13_ref)+(ql13_ref))/((qe13_ref)+(qe13_alt)+(ql13_ref)+(ql13_alt))                
                qe13_rf=(qe13_ref)/((qe13_count))                
                ql13_rf=(ql13_ref)/((ql13_ref)+(ql13_alt))
                q13_rfd=(qe13_rf-ql13_rf)
                if (qe13_rf + ql13_rf)==0.0 or (qe13_rf + ql13_rf)==2.0:
                    Z_q13=0.0
                else:
                    at_qe13_rf=math.asin(qe13_rf**0.5)*2
                    at_ql13_rf=math.asin(ql13_rf**0.5)*2
                    var_at_q13_rfd=(1/(bulkSize))+(1/qe13_count)+(1/ql13_count)+BS_Q13
                    Z_q13=(at_qe13_rf - at_ql13_rf)/(var_at_q13_rfd**0.5)
                    Zs_q13.append(Z_q13)
                    if Z_q13 < 0.0:
                        p_q13=norm.cdf(Z_q13)*2
                    else:
                        p_q13=(1.0-norm.cdf(Z_q13))*2
                            
            if qe14_count> bulkT and ql14_count > bulkT:
                q14_rf=((qe14_ref)+(ql14_ref))/((qe14_ref)+(qe14_alt)+(ql14_ref)+(ql14_alt))                
                qe14_rf=(qe14_ref)/((qe14_count))                
                ql14_rf=(ql14_ref)/((ql14_ref)+(ql14_alt))
                q14_rfd=(qe14_rf-ql14_rf)
                if (qe14_rf + ql14_rf)==0.0 or (qe14_rf + ql14_rf)==2.0:
                    Z_q14=0.0
                else:
                    at_qe14_rf=math.asin(qe14_rf**0.5)*2
                    at_ql14_rf=math.asin(ql14_rf**0.5)*2
                    var_at_q14_rfd=(1/(bulkSize))+(1/qe14_count)+(1/ql14_count)+BS_Q14
                    Z_q14=(at_qe14_rf - at_ql14_rf)/(var_at_q14_rfd**0.5)
                    Zs_q14.append(Z_q14)
                    if Z_q14 < 0.0:
                        p_q14=norm.cdf(Z_q14)*2
                    else:
                        p_q14=(1.0-norm.cdf(Z_q14))*2

            if qe13_count > bulkT and qe14_count > bulkT and ql13_count > bulkT and ql14_count > bulkT:
                q_rf=(q14_rf+q13_rf)/2
                qe_rf=(qe14_rf+qe13_rf)/2
                ql_rf=(ql14_rf+ql13_rf)/2
                q_rfd=(qe_rf-ql_rf)
                q_rfd_bnYears= q14_rf - q13_rf
                qe_rfd_bnYears= qe14_rf - qe13_rf
                ql_rfd_bnYears= ql14_rf - ql13_rf
                if (qe_rf + ql_rf)==0.0 or (qe_rf + ql_rf)==2.0:
                    Z_q=0.0
                else:
                    at_qe_rf=math.asin(qe_rf**0.5)*2
                    at_ql_rf=math.asin(ql_rf**0.5)*2
                    var_at_q_rfd=(1/(bulkSize))+(1/qe13_count)+(1/ql13_count)+(1/qe14_count)+(1/ql14_count)+BS_Q
                    Z_q=(at_qe_rf - at_ql_rf)/(var_at_q_rfd**0.5)
                    Zs_q.append(Z_q)
                    if Z_q < 0.0:
                        p_q=norm.cdf(Z_q)*2
                    else:
                        p_q=(1.0-norm.cdf(Z_q))*2                            
            
            if Z_q13!=None and Z_br13!=None and Z_im13!=None:
                pop_num=3
                Z_cum13=Z_br13+Z_im13+Z_q13
                if Z_cum13 < 0.0:
                    p_Z_cum13=norm.cdf(Z_cum13,scale=(pop_num**0.5))*2
                else:
                    p_Z_cum13=(1.0-norm.cdf(Z_cum13,scale=(pop_num**0.5)))*2
            if Z_q13==None and Z_br13!=None and Z_im13!=None:
                pop_num=2
                Z_cum13=Z_br13+Z_im13
                if Z_cum13 < 0.0:
                    p_Z_cum13=norm.cdf(Z_cum13,scale=(pop_num**0.5))*2
                else:
                    p_Z_cum13=(1.0-norm.cdf(Z_cum13,scale=(pop_num**0.5)))*2
            if Z_q13!=None and Z_br13==None and Z_im13!=None:
                pop_num=2
                Z_cum13=Z_q13+Z_im13
                if Z_cum13 < 0.0:
                    p_Z_cum13=norm.cdf(Z_cum13,scale=(pop_num**0.5))*2
                else:
                    p_Z_cum13=(1.0-norm.cdf(Z_cum13,scale=(pop_num**0.5)))*2            
            if Z_q13!=None and Z_br13!=None and Z_im13==None:
                pop_num=2
                Z_cum13=Z_br13+Z_q13
                if Z_cum13 < 0.0:
                    p_Z_cum13=norm.cdf(Z_cum13,scale=(pop_num**0.5))*2
                else:
                    p_Z_cum13=(1.0-norm.cdf(Z_cum13,scale=(pop_num**0.5)))*2

            if Z_q14!=None and Z_im14!=None:
                pop_num=2
                Z_cum14=Z_q14+Z_im14
                if Z_cum14 < 0.0:
                    p_Z_cum14=norm.cdf(Z_cum14,scale=(pop_num**0.5))*2
                else:
                    p_Z_cum14=(1.0-norm.cdf(Z_cum14,scale=(pop_num**0.5)))*2

            if Z_q!=None and Z_br13!=None and Z_im!=None:
                pop_num=3
                Z_cum=Z_br13+Z_im+Z_q
                if Z_cum < 0.0:
                    p_Z_cum=norm.cdf(Z_cum,scale=(pop_num**0.5))*2
                else:
                    p_Z_cum=(1.0-norm.cdf(Z_cum,scale=(pop_num**0.5)))*2
            if Z_q==None and Z_br!=None and Z_im!=None:
                pop_num=2
                Z_cum=Z_br+Z_im
                if Z_cum < 0.0:
                    p_Z_cum=norm.cdf(Z_cum,scale=(pop_num**0.5))*2
                else:
                    p_Z_cum=(1.0-norm.cdf(Z_cum,scale=(pop_num**0.5)))*2
            if Z_q!=None and Z_br==None and Z_im!=None:
                pop_num=2
                Z_cum=Z_q+Z_im
                if Z_cum < 0.0:
                    p_Z_cum=norm.cdf(Z_cum,scale=(pop_num**0.5))*2
                else:
                    p_Z_cum=(1.0-norm.cdf(Z_cum,scale=(pop_num**0.5)))*2            
            if Z_q!=None and Z_br!=None and Z_im==None:
                pop_num=2
                Z_cum=Z_br+Z_q
                if Z_cum < 0.0:
                    p_Z_cum=norm.cdf(Z_cum,scale=(pop_num**0.5))*2
                else:
                    p_Z_cum=(1.0-norm.cdf(Z_cum,scale=(pop_num**0.5)))*2
                    
            ####
                    """frequency difference between years, frequency difference between bulks across years."""
                    
            if all(i > bulkT for i in [bre13_count,bre14_count,brl13_count,ime13_count,iml13_count,ime14_count,iml14_count,qe13_count,ql13_count,qe14_count,ql14_count]) == True:
                Fst_all=tvar((im_rf,q_rf,br_rf))/((im_rf+q_rf+br_rf)/3)*(1-((im_rf+q_rf+br_rf)/3))
                Fst_IM_BR=tvar((im_rf,br_rf))/((im_rf+br_rf)/2)*(1-((im_rf+br_rf)/2))
                Fst_IM_Q=tvar((im_rf,q_rf))/((im_rf+q_rf)/2)*(1-((im_rf+q_rf)/2))
                Fst_BR_Q=tvar((q_rf,br_rf))/((q_rf+br_rf)/2)*(1-((q_rf+br_rf)/2))
                if bre_count > bulkT and brl_count> bulkT and ime_count> bulkT and iml_count> bulkT and qe_count> bulkT and ql_count > bulkT:
                    ssfd=(br_rfd+im_rfd+q_rfd)/3
            elif all(i > bulkT for i in [bre13_count,bre14_count,brl13_count,ime13_count,iml13_count,ime14_count,iml14_count]) == True:
                Fst_IM_BR=tvar((im_rf,br_rf))/((im_rf+br_rf)/2)*(1-((im_rf+br_rf)/2))
                ssfd=(br_rfd+im_rfd)/2              
            elif all(i > bulkT for i in [bre13_count,bre14_count,brl13_count,qe13_count,ql13_count,qe14_count,ql14_count]) == True:
                Fst_BR_Q=tvar((q_rf,br_rf))/((q_rf+br_rf)/2)*(1-((q_rf+br_rf)/2))
                ssfd=(br_rfd+q_rfd)/2
            elif all(i > bulkT for i in [ime13_count,iml13_count,ime14_count,iml14_count,qe13_count,ql13_count,qe14_count,ql14_count]) == True:
                Fst_IM_Q=tvar((im_rf,q_rf))/((im_rf+q_rf)/2)*(1-((im_rf+q_rf)/2))  
                ssfd=(im_rfd+q_rfd)/2
            ###          
            if Z_q!=None and Z_br!=None and Z_im!=None:
                allwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_BR,Fst_IM_Q,Fst_BR_Q,q_rfd_bnYears,qe_rfd_bnYears,ql_rfd_bnYears,im_rfd_bnYears,ime_rfd_bnYears,iml_rfd_bnYears,br_rfd_bnYears,bre_rfd_bnYears])
            if Z_q==None and Z_br!=None and Z_im!=None:
                IMBRwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_BR,Fst_IM_Q,Fst_BR_Q,q_rfd_bnYears,qe_rfd_bnYears,ql_rfd_bnYears,im_rfd_bnYears,ime_rfd_bnYears,iml_rfd_bnYears,br_rfd_bnYears,bre_rfd_bnYears])
            if Z_q!=None and Z_br==None and Z_im!=None:
                IMQwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_BR,Fst_IM_Q,Fst_BR_Q,q_rfd_bnYears,qe_rfd_bnYears,ql_rfd_bnYears,im_rfd_bnYears,ime_rfd_bnYears,iml_rfd_bnYears,br_rfd_bnYears,bre_rfd_bnYears])
            if Z_q!=None and Z_br!=None and Z_im==None:
                QBRwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_BR,Fst_IM_Q,Fst_BR_Q,q_rfd_bnYears,qe_rfd_bnYears,ql_rfd_bnYears,im_rfd_bnYears,ime_rfd_bnYears,iml_rfd_bnYears,br_rfd_bnYears,bre_rfd_bnYears])
            if Z_q!=None and Z_br==None and Z_im==None:
                Qwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_BR,Fst_IM_Q,Fst_BR_Q,q_rfd_bnYears,qe_rfd_bnYears,ql_rfd_bnYears,im_rfd_bnYears,ime_rfd_bnYears,iml_rfd_bnYears,br_rfd_bnYears,bre_rfd_bnYears])
            if Z_q==None and Z_br!=None and Z_im==None:
                BRwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_BR,Fst_IM_Q,Fst_BR_Q,q_rfd_bnYears,qe_rfd_bnYears,ql_rfd_bnYears,im_rfd_bnYears,ime_rfd_bnYears,iml_rfd_bnYears,br_rfd_bnYears,bre_rfd_bnYears])
            if Z_q==None and Z_br==None and Z_im!=None:            
                IMwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_BR,Fst_IM_Q,Fst_BR_Q,q_rfd_bnYears,qe_rfd_bnYears,ql_rfd_bnYears,im_rfd_bnYears,ime_rfd_bnYears,iml_rfd_bnYears,br_rfd_bnYears,bre_rfd_bnYears])

#            if Z_q13!=None and Z_br13!=None and Z_im13!=None and (Z_q==None or Z_im==None or Z_br==None):
#                allwriter13.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_BR,Fst_IM_Q,Fst_BR_Q])
#            if Z_q==None and Z_br!=None and Z_im!=None:
#                IMBRwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,pop_num,Z_cum,p_Z_cum,Fst_IM_BR])
#            if Z_q!=None and Z_br==None and Z_im!=None:
#                IMQwriter.writerow([scaff,pos,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_Q])
#            if Z_q!=None and Z_br!=None and Z_im==None:
#                QBRwriter.writerow([scaff,pos,,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_BR_Q])
#            if Z_q!=None and Z_br==None and Z_im==None:
#                Qwriter.writerow([scaff,pos,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_q,p_q])
#            if Z_q==None and Z_br!=None and Z_im==None:
#                BRwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,bre_ref,bre_alt,brl_ref,brl_alt,Z_br13,p_br13,Z_br,p_br])
#            if Z_q==None and Z_br==None and Z_im!=None:            
#                IMwriter.writerow([scaff,pos,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im])
#
#            if Z_q!=None and Z_br!=None and Z_im!=None:
#                allwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_BR,Fst_IM_Q,Fst_BR_Q])
#            if Z_q==None and Z_br!=None and Z_im!=None:
#                IMBRwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,pop_num,Z_cum,p_Z_cum,Fst_IM_BR])
#            if Z_q!=None and Z_br==None and Z_im!=None:
#                IMQwriter.writerow([scaff,pos,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_IM_Q])
#            if Z_q!=None and Z_br!=None and Z_im==None:
#                QBRwriter.writerow([scaff,pos,,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_br13,p_br13,Z_br14,p_br14,Z_br,p_br,Z_q13,p_q13,Z_q14,p_q14,Z_q,p_q,pop_num,Z_cum,p_Z_cum,Fst_BR_Q])
#            if Z_q!=None and Z_br==None and Z_im==None:
#                Qwriter.writerow([scaff,pos,qe13_ref,qe13_alt,ql13_ref,ql13_alt,qe14_ref,qe14_alt,ql14_ref,ql14_alt,qe_ref,qe_alt,ql_ref,ql_alt,Z_q,p_q])
#            if Z_q==None and Z_br!=None and Z_im==None:
#                BRwriter.writerow([scaff,pos,bre13_ref,bre13_alt,brl13_ref,brl13_alt,bre14_ref,bre14_alt,bre_ref,bre_alt,brl_ref,brl_alt,Z_br13,p_br13,Z_br,p_br])
#            if Z_q==None and Z_br==None and Z_im!=None:            
#                IMwriter.writerow([scaff,pos,ime13_ref,ime13_alt,iml13_ref,iml13_alt,ime14_ref,ime14_alt,iml14_ref,iml14_alt,ime_ref,ime_alt,iml_ref,iml_alt,Z_im13,p_im13,Z_im14,p_im14,Z_im,p_im])
            
    print "Zs_Quarry IQR = ", mstats.mquantiles(Zs_q,prob=[0.25,0.5,0.75])
    print "Zs_BR IQR = ", mstats.mquantiles(Zs_br,prob=[0.25,0.5,0.75])
    print "Zs_IM IQR = ", mstats.mquantiles(Zs_im,prob=[0.25,0.5,0.75])
    print "Zs_Quarry13 IQR = ", mstats.mquantiles(Zs_q13,prob=[0.25,0.5,0.75])
    print "Zs_BR13 IQR = ", mstats.mquantiles(Zs_br13,prob=[0.25,0.5,0.75])
    print "Zs_IM13 IQR = ", mstats.mquantiles(Zs_im13,prob=[0.25,0.5,0.75]) 
    print "Zs_Quarry14 IQR = ", mstats.mquantiles(Zs_q14,prob=[0.25,0.5,0.75])
    print "Zs_BR14 IQR = ", mstats.mquantiles(Zs_br14,prob=[0.25,0.5,0.75])
    print "Zs_IM14 IQR = ", mstats.mquantiles(Zs_im14,prob=[0.25,0.5,0.75])
     
          #Fst_IM_BR,Fst_IM_Q,Fst_BR_Q,Fst_all,ssfd     
          #"Fst_IM_BR","Fst_IM_Q","Fst_BR_Q","Fst_all","ssfd"
              
                        
                
                