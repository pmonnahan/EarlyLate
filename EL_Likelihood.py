# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:37:06 2015

Important Note:  I have stupid, inconsistent variable names.  For the likelihood, ll_Y_im is the likelihood in which bulk is constrained (i.e. year is allowed to vary, which will be compared to a model in which both yeah AND bulk car vary...hence the Y in lrt_Y)
                I changed the letters between the likelihoods and the likelihood ratios because if ll_Y_br (year allowed to vary) is compared to a model in which year and bulk are unconstrained (ll_C = 0.0), this provides a test of the effect of the bulk. 

@author: patrick
  """
from scipy.stats import *
import csv
import math
import time
import timeit

start = timeit.default_timer()

INPUT_FILE="/Users/patrick/Documents/Research/EarlyLate/counts_massiveJan16_AF250.txt"
OUTDIR="/Volumes/TOSHIBA EXT/EarlyLate/"
#OUTDIR="/Users/patrick/Documents/Research/EarlyLate/"

timestr = time.strftime("%Y%m%d-%H%M")
file1=0
file2=0
file3=0
file4=0
file5=0
file6=1

#filters
p_min1=0.05
p_max1=0.95
p_min=2*math.asin(p_min1**0.5)
p_max=2*math.asin(p_max1**0.5)
min_cov=10
max_cov=250
years="2013"
window_size=3000000

   
FDR=0.1



jjj=0
if file1 !=0:
    EL_Likelihoods=open(OUTDIR+"EL_Likelihoods_" + years + timestr + ".csv","wb")
    like=csv.writer(EL_Likelihoods,delimiter=",",dialect='excel')
    like.writerow(["scaff","pos","BRE13_R","BRE13_A","BRL13_R","BRL13_A","BRE14_R","BRE14_A","IME13_R","IME13_A","IML13_R","IML13_A","IME14_R","IME14_A","IML14_R","IML14_A","QE13_R","QE13_A","QL13_R","QL13_A","QE14_R","QE14_A","QL14_R","QL14_A","ll_S","ll_Y","ll_P","ll_B","ll_YP","ll_YB","ll_PB","ll_C","lrt_Y","lrt_P","lrt_B","p_Y","p_P","p_B"])
    #like.writerow(["scaff","pos","BRE13_R","BRE13_A","BRL13_R","BRL13_A","BRE14_R","BRE14_A","IME13_R","IME13_A","IML13_R","IML13_A","IME14_R","IME14_A","IML14_R","IML14_A","QE13_R","QE13_A","QL13_R","QL13_A","QE14_R","QE14_A","QL14_R","QL14_A","ll_S","ll_Y","ll_P","ll_B","ll_YP","ll_YB","ll_PB","lrt_Y","p_Y","lrt_y","p_y","lrt_P","p_P","lrt_p","p_p","lrt_B","p_B","lrt_b","p_b"])
if file2 != 0:
    out2=open(OUTDIR+"EL_Likelihoods_SimplePops_test" + timestr + ".csv","wb")
    out2x=csv.writer(out2,delimiter=",",dialect='excel')
    out2x.writerow(["scaff","pos","BRE13_R","BRE13_A","BRL13_R","BRL13_A","BRE14_R","BRE14_A","IME13_R","IME13_A","IML13_R","IML13_A","IME14_R","IME14_A","IML14_R","IML14_A","QE13_R","QE13_A","QL13_R","QL13_A","QE14_R","QE14_A","QL14_R","QL14_A","ll_S","ll_P","X2","p_p"])
if file3 != 0:
    out3=open(OUTDIR+"EL_ZandVs_" + timestr + ".csv","wb")
    out3x=csv.writer(out3,delimiter=",",dialect='excel')
    out3x.writerow(["scaff","pos","BRE13_R","BRE13_A","BRL13_R","BRL13_A","BRE14_R","BRE14_A","IME13_R","IME13_A","IML13_R","IML13_A","IME14_R","IME14_A","IML14_R","IML14_A","QE13_R","QE13_A","QL13_R","QL13_A","QE14_R","QE14_A","QL14_R","QL14_A","zbre13","zbrl13","zbre14","zime13","ziml13","zime14","ziml14","zqe13","zql13","zqe14","zql14","vbre13","vbrl13","vbre14","vime13","viml13","vime14","viml14","vqe13","vql13","vqe14","vql14"])
               
if file4 != 0:
    out4=open(OUTDIR+"EL_Likelihoods_PerPop_" + years + timestr + ".csv","wb")
    out4x=csv.writer(out4,delimiter=",",dialect='excel')    
    out4x.writerow(["scaff","pos","BRE13_R","BRE13_A","BRL13_R","BRL13_A","BRE14_R","BRE14_A","IME13_R","IME13_A","IML13_R","IML13_A","IME14_R","IME14_A","IML14_R","IML14_A","QE13_R","QE13_A","QL13_R","QL13_A","QE14_R","QE14_A","QL14_R","QL14_A","lrt_Y_br","df_Y_br","p_Y_br","lrt_B_br","df_B_br","p_B_br","lrt_Y_im","df_Y_im","p_Y_im","lrt_B_im","df_B_im","p_B_im","lrt_Y_q","df_Y_q","p_Y_q","lrt_B_q","df_B_q","p_B_q","lrt_YCB_im","df_YCB_im","p_YCB_im","lrt_YCB_q","df_YCB_q","p_YCB_q"])
if file5 != 0:
    out5BR=open(OUTDIR+"EL_Likelihoods_Windows_BR" + timestr + ".csv","wb")
    out5xBR=csv.writer(out5BR,delimiter=",",dialect='excel')    
    out5xBR.writerow(["scaff","window","start_bp","wind_sites","dist","numYtests","numBtests","sigY","sigB","BRE13_R","BRE13_A","BRL13_R","BRL13_A","BRE14_R","BRE14_A","lrt_Y","lrt_B","p_Y","p_B"])
    out5IM=open("/Users/patrick/Documents/Research/EarlyLate/EL_Likelihoods_Windows_IM" + timestr + ".csv","wb")
    out5xIM=csv.writer(out5IM,delimiter=",",dialect='excel')    
    out5xIM.writerow(["scaff","window","start_bp","wind_sites","dist","numYtests","numBtests","numYCBtests","sigY","sigB","sigYCB","IME13_R","IME13_A","IML13_R","IML13_A","IME14_R","IME14_A","IML14_R","IML14_A","lrt_Y","lrt_B","lrt_YCB","p_Y","p_B","p_YCB"])
    out5Q=open("/Users/patrick/Documents/Research/EarlyLate/EL_Likelihoods_Windows_Q" + timestr + ".csv","wb")
    out5xQ=csv.writer(out5Q,delimiter=",",dialect='excel')    
    out5xQ.writerow(["scaff","window","start_bp","wind_sites","dist","numYtests","numBtests","numYCBtests","sigY","sigB","sigYCB","QE13_R","QE13_A","QL13_R","QL13_A","QE14_R","QE14_A","QL14_R","QL14_A","lrt_Y","lrt_B","lrt_YCB","p_Y","p_B","p_YCB"])
if file6!=0:
    out6=open(OUTDIR+"EL_Likelihoods_PerPopSigCounts_3MbWindows" + timestr + ".csv","wb")
    out6x=csv.writer(out6,delimiter=",",dialect='excel')  
    out6x.writerow(["scaff","window","start_bp","wind_sites","dist","sigBim","sigBq","sigBbr"])


#USED ONLY IN WINDOW CALCULATIONS FOR NUM SIG SNPS IN A WINDOW...BASED ON PREVIOUS RUN AND FDR CALC
p_Y_im_cutoff = 0.000037
p_YCB_im_cutoff = 0.000007
p_B_im_cutoff = 4.62448137005e-06
p_Y_br_cutoff = 0.000010
p_B_br_cutoff = 1.60241798643e-06
p_Y_q_cutoff = 0.000106
p_B_q_cutoff = 0.000138586071041
p_YCB_q_cutoff = 0.000017


#BS variance factor : Paired + Single
vb1=0.034805286266108827198
vb2=0.036184480292554409286
vb3=0.022731757864254985291

vi1=0.0137298312949835840668 #With new data the variance calculations were similar to old data except for 2014...Early var was similar to late var and vice versa.
vi2=0.0184233013099333840790 #Should run through calc with old data to confirm that there wasn't an error
vi3=0.0066893073640574111927
vi4=0.0107630749671819188340

vq1=0.0152782122396757438776
vq2=0.0212079200267419748505
vq3=0.0082510518410872438211
vq4=0.0130222983597531889038

with open(INPUT_FILE,"rb") as sites_file:   
    sites=0
    fixed=0
    filtered=0
    dist=0
    pos=0
    wind_sites=0
    wind_df=0
    pvals_Y_br=[]
    pvals_Y_im=[]
    pvals_Y_q=[]
    pvals_B_br=[]
    pvals_B_im=[]
    pvals_B_q=[]
    pvals_YCB_br=[]
    pvals_YCB_im=[]
    pvals_YCB_q=[]
    pvals_Y=[]
    pvals_P=[]
    pvals_B=[]
    pvals_Y_windows=[]
    pvals_P_windows=[]
    pvals_B_windows=[]
    exceptions=0
    for i,site in enumerate(sites_file):
#        if i == 0:
#            print site
        if i>0:
            site=site.strip("\n")
            site=site.split("\t")
            scaff=site[0]
            scaff=scaff
            pos=float(site[1])
            site=[float(a) for a in site[4:]]
            bre13r,bre13a,brl13r,brl13a,ime14r,ime14a,iml14r,iml14a,ime13r,ime13a,iml13r,iml13a,qe14r,qe14a,ql14r,ql14a,qe13r,qe13a,ql13r,ql13a=site
            
            bre14r=0.0
            bre14a=0.0
            xbre14=-99.0
            vbre14=-99.0
            
            if bre13r+bre13a <= min_cov or bre13r+bre13a >= max_cov or years == "2014":
                xbre13=-99.0
                vbre13=-99.0
            else:
                xbre13=2*math.asin((bre13r/(bre13r+bre13a))**0.5)
                vbre13=(vb1+(1/(bre13r+bre13a)))
                
            if brl13r+brl13a <= min_cov or brl13r+brl13a >= max_cov or years == "2014":
                xbrl13=-99.0
                vbrl13=-99.0
            else:
                xbrl13=2*math.asin((brl13r/(brl13r+brl13a))**0.5)
                vbrl13=(vb2+(1/(brl13r+brl13a)))
                
            if ime13r+ime13a <= min_cov or ime13r+ime13a >= max_cov or years == "2014":
                xime13=-99.0
                vime13=-99.0
            else:
                xime13=2*math.asin((ime13r/(ime13r+ime13a))**0.5)
                vime13=(vi1+(1/(ime13r+ime13a)))
                
            if iml13r+iml13a <= min_cov or iml13r+iml13a >= max_cov or years == "2014":
                ximl13=-99.0
                viml=-99.0
            else:
                ximl13=2*math.asin((iml13r/(iml13r+iml13a))**0.5)
                viml13=(vi2+(1/(iml13r+iml13a)))
                
            if ime14r+ime14a <= min_cov or ime14r+ime14a >= max_cov or years == "2013":
                xime14=-99.0
                vime14=-99.0
            else:
                xime14=2*math.asin((ime14r/(ime14r+ime14a))**0.5)
                vime14=(vi3+(1/(ime14r+ime14a)))
                
            if iml14r+iml14a <= min_cov or iml14r+iml14a >= max_cov or years == "2013":
                ximl14=-99.0
                viml14=-99.0
            else:
                ximl14=2*math.asin((iml14r/(iml14r+iml14a))**0.5)
                viml14=(vi4+(1/(iml14r+iml14a)))
                
            if qe13r+qe13a <= min_cov or qe13r+qe13a >= max_cov or years == "2014":
                xqe13=-99.0
                vqe13=-99.0
            else:
                xqe13=2*math.asin((qe13r/(qe13r+qe13a))**0.5)
                vqe13=(vq1+(1/(qe13r+qe13a)))
                
                
            if ql13r+ql13a <= min_cov or ql13r+ql13a >= max_cov or years == "2014":
                xql13=-99.0
                vql13=-99.0
            else:
                xql13=2*math.asin((ql13r/(ql13r+ql13a))**0.5)
                vql13=(vq2+(1/(ql13r+ql13a)))
                
            if qe14r+qe14a <= min_cov or qe14r+qe14a >= max_cov or years == "2013":
                xqe14=-99.0
                vqe14=-99.0
            else:
                xqe14=2*math.asin((qe14r/(qe14r+qe14a))**0.5)
                vqe14=(vq3+(1/(qe14r+qe14a)))
            
            if ql14r+ql14a <= min_cov or ql14r+ql14a >= max_cov or years == "2013":
                xql14=-99.0
                vql14=-99.0
            else:
                xql14=2*math.asin((ql14r/(ql14r+ql14a))**0.5)
                vql14=(vq4 +(1/(ql14r+ql14a)))                
                

            if (all(k==0.0 for k in [xime13,ximl13,xime14,ximl14,xbre13,xbrl13,xbre14,xqe13,xql13,xqe14,xql14]) or all(k==math.pi for k in [xime13,ximl13,xime14,ximl14,xbre13,xbrl13,xbre14,xqe13,xql13,xqe14,xql14])):
                fixed+=1
                filtered+=1
                
            else:
                sites+=1
                if sites%10000==0:
                    print sites
                    
                    
                """Model 1 = Simplest Model.  All bulks share one mean."""                        
                x=[xbre13,xbrl13,xbre14,xime13,ximl13,xime14,ximl14,xqe13,xql13,xqe14,xql14]
                v=[vbre13,vbrl13,vbre14,vime13,viml13,vime14,viml14,vqe13,vql13,vqe14,vql14]
                w=[1/vbre13,1/vbrl13,1/vbre14,1/vime13,1/viml13,1/vime14,1/viml14,1/vqe13,1/vql13,1/vqe14,1/vql14]            
                W=sum(w)            
                uS=sum([xx*(ww/W) for xx,ww in zip(x,w)])
                ll_S=0.0
                for j,xx in enumerate(x):
                    ll_S+=(-(xx-uS)**2)/(2*v[j])
                
                """Model 2 = Year effect.  No differences among populations or bulk"""
                x13=[xbre13,xbrl13,xime13,ximl13,xqe13,xql13]
                x14=[xbre14,xime14,ximl14,xqe14,xql14]
                v13=[vbre13,vbrl13,vime13,viml13,vqe13,vql13]
                v14=[vbre14,vime14,viml14,vqe14,vql14]  
                w13=[1/vbre13,1/vbrl13,1/vime13,1/viml13,1/vqe13,1/vql13]
                w14=[1/vbre14,1/vime14,1/viml14,1/vqe14,1/vql14] 
                u13=sum([xx*(ww/sum(w13)) for xx,ww in zip(x13,w13)])
                u14=sum([xx*(ww/sum(w14)) for xx,ww in zip(x14,w14)])
                
                ll_Y=0.0
                for j,xx in enumerate(x13):
                    ll_Y+=(-(xx-u13)**2)/(2*v13[j])
                for j,xx in enumerate(x14):
                    ll_Y+=(-(xx-u14)**2)/(2*v14[j])
                
                """Model 3 = Population effect.  No differences across years or bulk"""
                ll_P=0.0
                xbr=[xbre13,xbrl13,xbre14]
                vbr=[vbre13,vbrl13,vbre14]
                wbr=[1/vbre13,1/vbrl13,1/vbre14]
                xim=[xime13,ximl13,xime14,ximl14]
                vim=[vime13,viml13,vime14,viml14]
                wim=[1/vime13,1/viml13,1/vime14,1/viml14]
                xq=[xqe13,xql13,xqe14,xql14]
                vq=[vqe13,vql13,vqe14,vql14]
                wq=[1/vqe13,1/vql13,1/vqe14,1/vql14]
                ubr=sum([xx*(ww/sum(wbr)) for xx,ww in zip(xbr,wbr)])
                uim=sum([xx*(ww/sum(wim)) for xx,ww in zip(xim,wim)])
                uq=sum([xx*(ww/sum(wq)) for xx,ww in zip(xq,wq)])
                
                for j,xx in enumerate(xbr):
                    ll_P+=(-(xx-ubr)**2)/(2*vbr[j])
                for j,xx in enumerate(xim):
                    ll_P+=(-(xx-uim)**2)/(2*vim[j])
                for j,xx in enumerate(xq):
                    ll_P+=(-(xx-uq)**2)/(2*vq[j])
                
                """Model 4 = Bulk effect.  No differences across year or pops."""
                ll_B=0.0
                xe=[xbre13,xbre14,xime13,xime14,xqe13,xqe14]
                xl=[xbrl13,ximl13,ximl14,xql13,xql14]
                ve=[vbre13,vbre14,vime13,vime14,vqe13,vqe14]
                vl=[vbrl13,viml13,viml14,vql13,vql14]
                we=[1/vbre13,1/vbre14,1/vime13,1/vime14,1/vqe13,1/vqe14]
                wl=[1/vbrl13,1/viml13,1/viml14,1/vql13,1/vql14]
                ue=sum([xx*(ww/sum(we)) for xx,ww in zip(xe,we)])
                ul=sum([xx*(ww/sum(wl)) for xx,ww in zip(xl,wl)])
                
                for j,xx in enumerate(xe):
                    ll_B+=(-(xx-ue)**2)/(2*ve[j])
                for j,xx in enumerate(xl):
                    ll_B+=(-(xx-ul)**2)/(2*vl[j])
                
                """Model 5 = Year and Pop effect.  No differences among bulks.""" 
                    
                ll_YP=0.0
                ll_Y_br=0.0
                ll_Y_im=0.0
                ll_Y_q=0.0
                xbr13=[xbre13,xbrl13]
                xbr14=[xbre14]
                vbr13=[vbre13,vbrl13]
                vbr14=[vbre14]
                wbr13=[1/vbre13,1/vbrl13]
                wbr14=[1/vbre14]
                xim13=[xime13,ximl13]
                xim14=[xime14,ximl14]
                vim13=[vime13,viml13]
                vim14=[vime14,viml14]
                wim13=[1/vime13,1/viml13]
                wim14=[1/vime14,1/viml14]
                xq13=[xqe13,xql13]
                xq14=[xqe14,xql14]
                vq13=[vqe13,vql13]
                vq14=[vqe14,vql14]
                wq13=[1/vqe13,1/vql13]
                wq14=[1/vqe14,1/vql14]
                ubr13=sum([xx*(ww/sum(wbr13)) for xx,ww in zip(xbr13,wbr13)])
                ubr14=sum([xx*(ww/sum(wbr14)) for xx,ww in zip(xbr14,wbr14)])
                uim13=sum([xx*(ww/sum(wim13)) for xx,ww in zip(xim13,wim13)])
                uim14=sum([xx*(ww/sum(wim14)) for xx,ww in zip(xim14,wim14)])
                uq13=sum([xx*(ww/sum(wq13)) for xx,ww in zip(xq13,wq13)])            
                uq14=sum([xx*(ww/sum(wq14)) for xx,ww in zip(xq14,wq14)]) 
                
                df_B_br=1
                df_B_im=2
                df_B_q=2                    
                df_B=5

                if ((xbre13 < p_min and xbrl13 < p_min) or (xbre13 > p_max and xbrl13 > p_max) or bre13r+bre13a <= min_cov or brl13r+brl13a <= min_cov or bre13r+bre13a >= max_cov or brl13r+brl13a >= max_cov):
                    df_B_br-=1
                    df_B+=(-1)
                else:
                    for j,xx in enumerate(xbr13):                   
                        ll_Y_br+=(-(xx-ubr13)**2)/(2*vbr13[j])
                        ll_YP+=(-(xx-ubr13)**2)/(2*vbr13[j])
                        
                if ((xime13 < p_min and ximl13 < p_min) or (xime13 > p_max and ximl13 > p_max) or ime13r+ime13a <= min_cov or iml13r+iml13a <= min_cov or ime13r+ime13a >= max_cov or iml13r+iml13a >= max_cov):
                    df_B_im-=1
                    df_B+=(-1)
                else:
                    for j,xx in enumerate(xim13):
                        if xx==-99.0:
                            print "fuck"
                        ll_Y_im+=(-(xx-uim13)**2)/(2*vim13[j])
                        ll_YP+=(-(xx-uim13)**2)/(2*vim13[j])
                        
                if ((xqe13 < p_min and xql13 < p_min) or (xqe13 > p_max and xql13 > p_max) or qe13r+qe13a <= min_cov or ql13r+ql13a <= min_cov or qe13r+qe13a >= max_cov or ql13r+ql13a >= max_cov):
                    df_B_q-=1 
                    df_B+=(-1)
                else:
                    for j,xx in enumerate(xq13):
                        ll_Y_q+=(-(xx-uq13)**2)/(2*vq13[j])
                        ll_YP+=(-(xx-uq13)**2)/(2*vq13[j])
                    
                if ((xime14 < p_min and ximl14 < p_min) or (xime14 > p_max and ximl14 > p_max) or ime14r+ime14a <= min_cov or iml14r+iml14a <= min_cov or ime14r+ime14a >= max_cov or iml14r+iml14a >= max_cov):
                    df_B_im-=1
                    df_B+=(-1)
                else:
                    for j,xx in enumerate(xim14):
                        ll_Y_im+=(-(xx-uim14)**2)/(2*vim14[j])
                        ll_YP+=(-(xx-uim14)**2)/(2*vim14[j])
                        
                if ((xqe14 < p_min and xql14 < p_min) or (xqe14 > p_max and xql14 > p_max) or qe14r+qe14a <= min_cov or ql14r+ql14a <= min_cov or qe14r+qe14a >= max_cov or ql14r+ql14a >= max_cov):
                    df_B_q-=1
                    df_B+=(-1)
                else:
                    for j,xx in enumerate(xq14):                    
                        ll_Y_q+=(-(xx-uq14)**2)/(2*vq14[j])
                        ll_YP+=(-(xx-uq14)**2)/(2*vq14[j])
                                        
                
                """Model 6 = Year and bulk effect.  No differences across pops."""
                ll_YB=0.0
                xe13=[xbre13,xime13,xqe13]
                xe14=[xbre14,xime14,xqe14]
                xl13=[xbrl13,ximl13,xql13]
                xl14=[ximl14,xql14]
                ve13=[vbre13,vime13,vqe13]
                ve14=[vbre14,vime14,vqe14]
                vl13=[vbrl13,viml13,vql13]
                vl14=[viml14,vql14]
                we13=[1/vbre13,1/vime13,1/vqe13]
                we14=[1/vbre14,1/vime14,1/vqe14]
                wl13=[1/vbrl13,1/viml13,1/vql13]
                wl14=[1/viml14,1/vql14]
                ue13=sum([xx*(ww/sum(we13)) for xx,ww in zip(xe13,we13)])
                ue14=sum([xx*(ww/sum(we14)) for xx,ww in zip(xe14,we14)])
                ul13=sum([xx*(ww/sum(wl13)) for xx,ww in zip(xl13,wl13)])
                ul14=sum([xx*(ww/sum(wl14)) for xx,ww in zip(xl14,wl14)])
                
                df_P = 7
                
                if (any(k <= min_cov for k in [ime14r+ime14a,iml14r+iml14a,qe14r+qe14a,ql14r+ql14a])) or (any(k >= max_cov for k in [ime14r+ime14a,iml14r+iml14a,qe14r+qe14a,ql14r+ql14a])):
                    df_P -= 3
                else:
                    for j,xx in enumerate(xe14):
                        ll_YB+=(-(xx-ue14)**2)/(2*ve14[j])
                    for j,xx in enumerate(xl14):
                        ll_YB+=(-(xx-ul14)**2)/(2*vl14[j])
                if (any(k <= min_cov for k in [ime13r+ime13a,iml13r+iml13a,qe13r+qe13a,ql13r+ql13a,bre13r+bre13a,brl13r+brl13a])) or (any(k < max_cov for k in [ime13r+ime13a,iml13r+iml13a,qe13r+qe13a,ql13r+ql13a,bre13r+bre13a,brl13r+brl13a])):
                    df_P -= 4
                else:                                    
                    for j,xx in enumerate(xe13):
                        ll_YB+=(-(xx-ue13)**2)/(2*ve13[j])
                    for j,xx in enumerate(xl13):
                        ll_YB+=(-(xx-ul13)**2)/(2*vl13[j])
                
                """Model 6b = Year effect and consistent effect of bulk"""
                ll_YCB_im=0.0
                ll_YCB_q=0.0
                uim13 = -((-viml13*xime13 - vime13*ximl13 - vime14*ximl13 - viml14*ximl13 + viml13*xime14 - viml13*ximl14)/(vime13 + viml13 + vime14 + viml14))
                uim14 = -((viml14*xime13 - viml14*ximl13 - viml14*xime14 - vime13*ximl14 - viml13*ximl14 - vime14*ximl14)/(vime13 + viml13 + vime14 + viml14))                
                uq13 = -((-vql13*xqe13 - vqe13*xql13 - vqe14*xql13 - vql14*xql13 + vql13*xqe14 - vql13*xql14)/(vqe13 + vql13 + vqe14 + vql14))
                uq14 = -((vql14*xqe13 - vql14*xql13 - vql14*xqe14 - vqe13*xql14 - vql13*xql14 - vqe14*xql14)/(vqe13 + vql13 + vqe14 + vql14))                                
                aim = -((-vime14*xime13 - viml14*xime13 + vime14*ximl13 + viml14*ximl13 - vime13*xime14 - viml13*xime14 + 
                vime13*ximl14 + viml13*ximl14)/(vime13 + viml13 + vime14 + viml14))
                aq = -((-vqe14*xqe13 - vql14*xqe13 + vqe14*xql13 + vql14*xql13 - vqe13*xqe14 - vql13*xqe14 + 
                vqe13*xql14 + vql13*xql14)/(vqe13 + vql13 + vqe14 + vql14))
                
                df_YCB_br=1
                df_YCB_im=1
                df_YCB_q=1                    
                df_YCB=5
                        
                if (xime13 < p_min and ximl13 < p_min) or (xime13 > p_max and ximl13 > p_max) or (xime14 < p_min and ximl14 < p_min) or (xime14 > p_max and ximl14 > p_max) or (any(k <= min_cov for k in [ime13r+ime13a,iml13r+iml13a,ime14r+ime14a,iml14r+iml14a])) or (any(k >= max_cov for k in [ime13r+ime13a,iml13r+iml13a,ime14r+ime14a,iml14r+iml14a])):
                    df_YCB_im-=1
                    df_YCB+=(-1)
                else:
                    ll_YCB_im+=(-(xime13-(uim13+aim))**2)/(2*vime13)
                    ll_YCB_im+=(-(xime14-(uim14+aim))**2)/(2*vime14)
                    ll_YCB_im+=(-(ximl13-(uim13))**2)/(2*viml13)
                    ll_YCB_im+=(-(ximl14-(uim14))**2)/(2*viml14)
                        
                if (xqe13 < p_min and xql13 < p_min) or (xqe13 > p_max and xql13 > p_max) or (xqe14 < p_min and xql14 < p_min) or (xqe14 > p_max and xql14 > p_max) or (any(k <= min_cov for k in [qe13r+qe13a,ql13r+ql13a,qe14r+qe14a,ql14r+ql14a])) or (any(k >= max_cov for k in [qe13r+qe13a,ql13r+ql13a,qe14r+qe14a,ql14r+ql14a])):
                    df_YCB_q+=(-1) 
                    df_YCB+=(-1)
                else:
                    ll_YCB_q+=(-(xqe13-(uq13+aq))**2)/(2*vqe13)
                    ll_YCB_q+=(-(xqe14-(uq14+aq))**2)/(2*vqe14)
                    ll_YCB_q+=(-(xql13-(uq13))**2)/(2*vql13)
                    ll_YCB_q+=(-(xql14-(uq14))**2)/(2*vql14)

    
                """Model 7 = Pop and bulk effect.  No effect of year."""
                ll_PB=0.0
                ll_B_br=0.0
                ll_B_im=0.0
                ll_B_q=0.0
                xbre=[xbre13,xbre14]
                xbrl=[xbrl13]
                vbre=[vbre13,vbre14]
                vbrl=[vbrl13]
                wbre=[1/vbre13,1/vbre14]
                wbrl=[1/vbrl13]
                ubre=sum([xx*(ww/sum(wbre)) for xx,ww in zip(xbre,wbre)])
                ubrl=sum([xx*(ww/sum(wbrl)) for xx,ww in zip(xbrl,wbrl)])
                xime=[xime13,xime14]
                ximl=[ximl13,ximl14]
                vime=[vime13,vime14]
                viml=[viml13,viml14]
                wime=[1/vime13,1/vime14]
                wiml=[1/viml13,1/viml14]
                uime=sum([xx*(ww/sum(wime)) for xx,ww in zip(xime,wime)])
                uiml=sum([xx*(ww/sum(wiml)) for xx,ww in zip(ximl,wiml)])
                xqe=[xqe13,xqe14]
                xql=[xql13,xql14]
                vqe=[vqe13,vqe14]
                vql=[vql13,vql14]
                wqe=[1/vqe13,1/vqe14]
                wql=[1/vql13,1/vql14]
                uqe=sum([xx*(ww/sum(wqe)) for xx,ww in zip(xqe,wqe)])
                uql=sum([xx*(ww/sum(wql)) for xx,ww in zip(xql,wql)])
                
                df_Y=5                                
                df_Y_br=1
                df_Y_im=2
                df_Y_q=2
                    
                if ((xbre13 < p_min and xbrl13 < p_min and xbre14 < p_min) or (xbre13 > p_max and xbrl13 > p_max and xbre14 > p_max) or (any(k <= min_cov for k in [bre13r+bre13a,brl13r+brl13a,bre14r+bre14a])) or (any(k >= max_cov for k in [bre13r+bre13a,brl13r+brl13a,bre14r+bre14a]))):
                    df_Y_br-=1
                    df_Y-=1
                else:
                    for j,xx in enumerate(xbre):
                        ll_B_br+=(-(xx-ubre)**2)/(2*vbre[j])
                        ll_PB+=(-(xx-ubre)**2)/(2*vbre[j])
                    for j,xx in enumerate(xbrl):
                        ll_B_br+=(-(xx-ubrl)**2)/(2*vbrl[j])
                        ll_PB+=(-(xx-ubrl)**2)/(2*vbrl[j])
                        
                if ((xime13 < p_min and ximl13 < p_min and xime14 < p_min and ximl14 < p_min) or (xime13 > p_max and ximl13 > p_max and xime14 > p_max and ximl14 > p_max) or (any(k <= min_cov for k in [ime13r+ime13a,iml13r+iml13a,ime14r+ime14a,iml14r+iml14a])) or (any(k >= max_cov for k in [ime13r+ime13a,iml13r+iml13a,ime14r+ime14a,iml14r+iml14a]))):
                    df_Y_im-=2
                    df_Y-=2
                else:
                    for j,xx in enumerate(xime):
                        ll_B_im+=(-(xx-uime)**2)/(2*vime[j])
                        ll_PB+=(-(xx-uime)**2)/(2*vime[j])
                    for j,xx in enumerate(ximl):
                        ll_B_im+=(-(xx-uiml)**2)/(2*viml[j])
                        ll_PB+=(-(xx-uiml)**2)/(2*viml[j])
                
                if ((xqe13 < p_min and xql13 < p_min and xqe14 < p_min and xql14 < p_min) or (xqe13 > p_max and xql13 > p_max and xqe14 > p_max and xql14 > p_max) or (any(k <= min_cov for k in [qe13r+qe13a,ql13r+ql13a,qe14r+qe14a,ql14r+ql14a])) or (any(k >= max_cov for k in [qe13r+qe13a,ql13r+ql13a,qe14r+qe14a,ql14r+ql14a]))):
                    df_Y_q-=2
                    df_Y-=2
                else:
                    for j,xx in enumerate(xqe):
                        ll_B_q+=(-(xx-uqe)**2)/(2*vqe[j])
                        ll_PB+=(-(xx-uqe)**2)/(2*vqe[j])
                    for j,xx in enumerate(xql):
                        ll_B_q+=(-(xx-uql)**2)/(2*vql[j])
                        ll_PB+=(-(xx-uql)**2)/(2*vql[j]) 
                
                """Model 8 = Most complex model.  All bulks have distinct mean."""
                ll_C=0.0
                
                """tests for each population individually"""                
                if df_Y_br > 0:
                    lrt_Y_br=-2*ll_B_br
                    p_Y_br=chisqprob(lrt_Y_br,df_Y_br)
                    pvals_Y_br.append(p_Y_br)
                else:
                    p_Y_br="-"
                    lrt_Y_br="-"
                if df_B_br > 0:
                    lrt_B_br=-2*ll_Y_br
                    p_B_br=chisqprob(lrt_B_br,df_B_br)
                    pvals_B_br.append(p_B_br)
                else:
                    p_B_br="-"
                    lrt_B_br="-"
                if df_Y_im > 0:
                    lrt_Y_im=-2*ll_B_im
                    p_Y_im=chisqprob(lrt_Y_im,df_Y_im)
                    pvals_Y_im.append(p_Y_im)
                else:
                    p_Y_im="-"
                    lrt_Y_im="-"
                if df_B_im > 0:
                    lrt_B_im=-2*ll_Y_im
                    p_B_im=chisqprob(lrt_B_im,df_B_im)
                    pvals_B_im.append(p_B_im)
                else:
                    p_B_im="-"
                    lrt_B_im="-"
                if df_Y_q > 0:
                    lrt_Y_q=-2*ll_B_q
                    p_Y_q=chisqprob(lrt_Y_q,df_Y_q)
                    pvals_Y_q.append(p_Y_q)
                else:
                    p_Y_q="-"
                    lrt_Y_q="-"
                if df_B_q > 0:
                    lrt_B_q=-2*ll_Y_q
                    p_B_q=chisqprob(lrt_B_q,df_B_q)
                    pvals_B_q.append(p_B_q)
                else:
                    p_B_q="-"
                    lrt_B_q="-"
                if (df_YCB_q > 0 and df_B_q==2):
                    lrt_YCB_q=-2*(ll_YCB_q)
                    p_YCB_q=chisqprob(lrt_YCB_q,df_YCB_q)
                    pvals_YCB_q.append(p_YCB_q)
                else:
                    p_YCB_q="-"
                    lrt_YCB_q="-"
                if (df_YCB_im > 0 and df_B_im==2):
                    lrt_YCB_im=-2*(ll_YCB_im)
                    p_YCB_im=chisqprob(lrt_YCB_im,df_YCB_im)
                    pvals_YCB_im.append(p_YCB_im)
                else:
                    p_YCB_im="-"
                    lrt_YCB_im="-"
                
                
                """test for all pops together"""                
                if df_Y > 0:
                    lrt_Y=-2*(ll_PB-ll_C)
                    p_Y=chisqprob(lrt_Y,df_Y)
                    pvals_Y.append(p_Y)
                else:
                    lrt_Y="-"
                    p_Y="-"
                if df_B > 0:
                    lrt_B=-2*(ll_YP-ll_C)
                    p_B=chisqprob(lrt_B,df_B)
                    pvals_B.append(p_B)
                else:
                    lrt_B="-"
                    p_B="-"
                if df_P > 0:    
                    lrt_P=-2*(ll_YB-ll_C)
                    p_P=chisqprob(lrt_P,df_P)
                    pvals_P.append(p_P)
                else:
                    lrt_P="-"
                    p_P="-"
                
                """Print to files 1 through 4"""
                if file1!=0:
                    like.writerow([scaff,pos,bre13r,bre13a,brl13r,brl13a,bre14r,bre14a,ime13r,ime13a,iml13r,iml13a,ime14r,ime14a,iml14r,iml14a,qe13r,qe13a,ql13r,ql13a,qe14r,qe14a,ql14r,ql14a,ll_S,ll_Y,ll_P,ll_B,ll_YP,ll_YB,ll_PB,ll_C,lrt_Y,lrt_P,lrt_B,p_Y,p_P,p_B])   
                    #like.writerow([scaff,pos,bre13r,bre13a,brl13r,brl13a,bre14r,bre14a,ime13r,ime13a,iml13r,iml13a,ime14r,ime14a,iml14r,iml14a,qe13r,qe13a,ql13r,ql13a,qe14r,qe14a,ql14r,ql14a,ll_S,ll_Y,ll_P,ll_B,ll_YP,ll_YB,ll_PB,lrt_Y,p_Y,lrt_y,p_y,lrt_P,p_P,lrt_p,p_p,lrt_B,p_B,lrt_b,p_b])   
                if file2 != 0:
                    out2x.writerow([scaff,pos,bre13r,bre13a,brl13r,brl13a,bre14r,bre14a,ime13r,ime13a,iml13r,iml13a,ime14r,ime14a,iml14r,iml14a,qe13r,qe13a,ql13r,ql13a,qe14r,qe14a,ql14r,ql14a,ll_S,ll_P,lrt_p,p_p])
                if file3 != 0:
                    out3x.writerow([scaff,pos,bre13r,bre13a,brl13r,brl13a,bre14r,bre14a,ime13r,ime13a,iml13r,iml13a,ime14r,ime14a,iml14r,iml14a,qe13r,qe13a,ql13r,ql13a,qe14r,qe14a,ql14r,ql14a,xbre13,xbrl13,xbre14,xime13,ximl13,xime14,ximl14,xqe13,xql13,xqe14,xql14,vbre13,vbrl13,vbre14,vime13,viml13,vime14,viml14,vqe13,vql13,vqe14,vql14])
                if file4 != 0:
                    out4x.writerow([scaff,pos,bre13r,bre13a,brl13r,brl13a,bre14r,bre14a,ime13r,ime13a,iml13r,iml13a,ime14r,ime14a,iml14r,iml14a,qe13r,qe13a,ql13r,ql13a,qe14r,qe14a,ql14r,ql14a,lrt_Y_br,df_Y_br,p_Y_br,lrt_B_br,df_B_br,p_B_br,lrt_Y_im,df_Y_im,p_Y_im,lrt_B_im,df_B_im,p_B_im,lrt_Y_q,df_Y_q,p_Y_q,lrt_B_q,df_B_q,p_B_q,lrt_YCB_im,df_YCB_im,p_YCB_im,lrt_YCB_q,df_YCB_q,p_YCB_q])
                
                
                """determine windows and print to file"""
                if i==1:
                    start_bp=pos
                    wind_sites=1
                    sigBim=0
                    sigBq=0
                    sigBbr=0
                    sigYim=0
                    sigYq=0
                    sigYbr=0
                    sigYCBim=0
                    sigYCBq=0
                    window=1
                    IME13R=ime13r
                    IME13A=ime13a
                    IML13R=iml13r
                    IML13A=iml13a
                    IME14R=ime14r
                    IME14A=ime14a
                    IML14R=iml14r
                    IML14A=iml14a
                    BRE13R=bre13r
                    BRE13A=bre13a
                    BRL13R=brl13r
                    BRL13A=brl13a
                    BRE14R=bre14r
                    BRE14A=bre14a
                    QE13R=qe13r
                    QE13A=qe13a
                    QL13R=ql13r
                    QL13A=ql13a
                    QE14R=qe14r
                    QE14A=qe14a
                    QL14R=ql14r
                    QL14A=ql14a                   
                                      
                    if lrt_Y_br == "-":
                        LRT_Y_br=0.0
                        numYtests_br=0
                    else:
                        LRT_Y_br=lrt_Y_br
                        wind_sites=1
                        numYtests_br=1
                    if lrt_B_br == "-":
                        LRT_B_br=0.0
                        numBtests_br=0
                    else:
                        LRT_B_br=lrt_B_br
                        wind_sites=1
                        numBtests_br=1
                    if lrt_Y_im == "-" or df_Y_im < 2:
                        LRT_Y_im=0.0
                        numYtests_im=0
                    else:
                        LRT_Y_im=lrt_Y_im
                        wind_sites=1
                        numYtests_im=1
                    if lrt_B_im == "-" or df_B_im < 2:
                        LRT_B_im=0.0
                        numBtests_im=0
                    else:
                        LRT_B_im=lrt_B_im
                        wind_sites=1
                        numBtests_im=1
                    if lrt_Y_q == "-" or df_Y_q < 2:
                        LRT_Y_q=0.0
                        numYtests_q=0
                    else:
                        LRT_Y_q=lrt_Y_q
                        wind_sites=1
                        numYtests_q=1
                    if lrt_B_q == "-" or df_B_q < 2:
                        LRT_B_q=0.0
                        numBtests_q=0
                    else:
                        LRT_B_q=lrt_B_q
                        wind_sites=1
                        numBtests_q=1
                    if lrt_YCB_q == "-":
                        LRT_YCB_q=0.0
                        numYCBtests_q=0
                    else:
                        LRT_YCB_q=lrt_YCB_q
                        numYCBtests_q=1
                    if lrt_YCB_im == "-":
                        LRT_YCB_im=0.0
                        numYCBtests_im=0
                    else:
                        LRT_YCB_im=lrt_YCB_im
                        numYCBtests_im=1
                        
                elif i>1:                      
                    if pos-start_bp>window_size or scaff!=prev_scaff: # or wind_sites == 5
                        if file5 != 0:
                            out5xQ.writerow([scaff,window,start_bp,wind_sites,dist,numYtests_q,numBtests_q,numYCBtests_q,sigYq,sigBq,sigYCBq,QE13R,QE13A,QL13R,QL13A,QE14R,QE14A,QL14R,QL14A,LRT_Y_q,LRT_B_q,LRT_YCB_q])
                            out5xBR.writerow([scaff,window,start_bp,wind_sites,dist,numYtests_br,numBtests_br,sigYbr,sigBbr,BRE13R,BRE13A,BRL13R,BRL13A,BRE14R,BRE14A,LRT_Y_br,LRT_B_br])
                            out5xIM.writerow([scaff,window,start_bp,wind_sites,dist,numYtests_im,numBtests_im,numYCBtests_im,sigYim,sigBim,sigYCBim,IME13R,IME13A,IML13R,IML13A,IME14R,IME14A,IML14R,IML14A,LRT_Y_im,LRT_B_im,LRT_YCB_im])
                        if file6 != 0:
                            if sigBim > 0:
                                sigBim=1
                            if sigBq > 0:
                                sigBq=1
                            if sigBbr > 0:
                                sigBbr=1
                            out6x.writerow([scaff,window,start_bp,wind_sites,dist,sigBim,sigBq,sigBbr])
                        if scaff!=prev_scaff:
                            window=0
                        else:    
                            window+=1
                        start_bp=pos                    
                        wind_sites=1
                        sigBim=0
                        sigBq=0
                        sigBbr=0
                        sigYim=0
                        sigYq=0
                        sigYbr=0
                        sigYCBim=0
                        sigYCBq=0
                        IME13R=ime13r
                        IME13A=ime13a
                        IML13R=iml13r
                        IML13A=iml13a
                        IME14R=ime14r
                        IME14A=ime14a
                        IML14R=iml14r
                        IML14A=iml14a
                        BRE13R=bre13r
                        BRE13A=bre13a
                        BRL13R=brl13r
                        BRL13A=brl13a
                        BRE14R=bre14r
                        BRE14A=bre14a
                        QE13R=qe13r
                        QE13A=qe13a
                        QL13R=ql13r
                        QL13A=ql13a
                        QE14R=qe14r
                        QE14A=qe14a
                        QL14R=ql14r
                        QL14A=ql14a
                            
                        if lrt_Y_im == "-" or df_Y_im < 2:
                            LRT_Y_im=0.0
                            numYtests_im=0
                        else:
                            LRT_Y_im=lrt_Y_im
                            numYtests_im=1
                            if p_Y_im < p_Y_im_cutoff:
                                sigYim+=1
                        if lrt_B_im == "-": #or df_B_im < 2
                            LRT_B_im=0.0
                            numBtests_im=0
                        else:
                            LRT_B_im=lrt_B_im
                            numBtests_im=1
                            if p_B_im < p_B_im_cutoff:
                                sigBim+=1
                        if lrt_YCB_im == "-":
                            LRT_YCB_im=0.0
                            numYCBtests_im=0
                        else:
                            LRT_YCB_im=lrt_YCB_im 
                            numYCBtests_im=1
                            if p_YCB_im < p_YCB_im_cutoff:
                                sigYCBim+=1
                        if lrt_Y_br == "-":
                            LRT_Y_br=0.0
                            numYtests_br=0
                        else:
                            LRT_Y_br=lrt_Y_br
                            numYtests_br=1
                            if p_Y_br < p_Y_br_cutoff:
                                sigYbr+=1
                        if lrt_B_br == "-":
                            LRT_B_br=0.0
                            numBtests_br=0
                        else:
                            LRT_B_br=lrt_B_br                
                            numBtests_br=1
                            if p_B_br < p_B_br_cutoff:
                                sigBbr+=1
                        if lrt_Y_q == "-" or df_Y_q < 2:
                            LRT_Y_q=0.0
                            numYtests_q=0
                        else:
                            LRT_Y_q=lrt_Y_q
                            numYtests_q=1
                            if p_Y_q < p_Y_q_cutoff:
                                sigYq+=1
                        if lrt_B_q == "-": #or df_B_q < 2
                            LRT_B_q=0.0
                            numBtests_q=0
                        else:
                            LRT_B_q=lrt_B_q
                            numBtests_q=1
                            if p_B_q < p_B_q_cutoff:
                                sigBq+=1
                        if lrt_YCB_q == "-":
                            LRT_YCB_q=0.0
                            numYCBtests_q=0
                        else:
                            LRT_YCB_q=lrt_YCB_q
                            numYCBtests_q=1
                            if p_YCB_q < p_YCB_q_cutoff:
                                sigYCBq+=1
                            
                    else:               
                        wind_sites+=1
                        dist=pos-start_bp
                        if p_B_im < prev_p_B_IM:
                            IME13R=ime13r #Only record read info for most significant site.  Best method to represent important allele frequency in a window?
                            IME13A=ime13a
                            IML13R=iml13r
                            IML13A=iml13a
                            IME14R=ime14r
                            IME14A=ime14a
                            IML14R=iml14r
                            IML14A=iml14a
                        if p_B_br < prev_p_B_BR:
                            BRE13R=bre13r
                            BRE13A=bre13a
                            BRL13R=brl13r
                            BRL13A=brl13a
                            BRE14R=bre14r
                            BRE14A=bre14a
                        if p_B_q < prev_p_B_Q:
                            QE13R=qe13r
                            QE13A=qe13a
                            QL13R=ql13r
                            QL13A=ql13a
                            QE14R=qe14r
                            QE14A=qe14a
                            QL14R=ql14r
                            QL14A=ql14a
                                                    
                        if lrt_Y_im == "-": #or df_Y_im < 2
                            pass
                        else:
                            LRT_Y_im+=lrt_Y_im                          
                            numYtests_im+=1
                            if p_Y_im < p_Y_im_cutoff:
                                sigYim+=1
                        if lrt_B_im == "-":
                            pass
                        else:
                            LRT_B_im+=lrt_B_im
                            numBtests_im+=1
                            if p_B_im < p_B_im_cutoff: #or df_B_im < 2
                                sigBim+=1
                        if lrt_YCB_im == "-":
                            pass
                        else:
                            LRT_YCB_im+=lrt_YCB_im
                            numYCBtests_im+=1
                            if p_YCB_im < p_YCB_im_cutoff:
                                sigYCBim+=1
                        if lrt_Y_br == "-":
                            pass
                        else:
                            LRT_Y_br+=lrt_Y_br
                            numYtests_br+=1
                            if p_Y_br < p_Y_br_cutoff:
                                sigYbr+=1
                        if lrt_B_br == "-":
                            pass
                        else:
                            LRT_B_br+=lrt_B_br
                            numBtests_br+=1
                            if p_B_br < p_B_br_cutoff:
                                sigBbr+=1
                        if lrt_Y_q == "-" or df_Y_q < 2:
                            pass
                        else:
                            LRT_Y_q+=lrt_Y_q
                            numYtests_q+=1
                            if p_Y_q < p_Y_q_cutoff:
                                sigYq+=1
                        if lrt_B_q == "-": #or df_B_q < 2
                            pass
                        else:
                            LRT_B_q+=lrt_B_q
                            numBtests_q+=1
                            if p_B_q < p_B_q_cutoff:
                                sigBq+=1
                        if lrt_YCB_q == "-":
                            pass
                        else:
                            LRT_YCB_q+=lrt_YCB_q
                            numYCBtests_q+=1
                            if p_YCB_q < p_YCB_q_cutoff:
                                sigYCBq+=1
                        
            prev_scaff=scaff
            prev_p_B_Q=p_B_q
            prev_p_B_IM=p_B_im
            prev_p_B_BR=p_B_br               
                

                   
    """Determine p-value tresholds given a specified FDR"""
    pvals_Y_br.sort()
    pvals_Y_im.sort()
    pvals_Y_q.sort()
    pvals_B_br.sort()
    pvals_B_im.sort()
    pvals_B_q.sort()
    pvals_YCB_im.sort()
    pvals_YCB_q.sort()
    pvals_Y.sort()
    pvals_P.sort()
    pvals_B.sort()
    pvals_Y_windows.sort()
    pvals_P_windows.sort()
    pvals_B_windows.sort()
    
    p_Y_cutoff=0.0
    p_B_cutoff=0.0
    p_P_cutoff=0.0
    p_Y_br_cutoff=0.0
    p_Y_im_cutoff=0.0
    p_Y_q_cutoff=0.0
    p_B_br_cutoff=0.0
    p_B_im_cutoff=0.0
    p_B_q_cutoff=0.0
    p_YCB_im_cutoff=0.0
    p_YCB_q_cutoff=0.0
    
    
    n=0
    while pvals_Y[n] < ((float(n)+1.0)/float(len(pvals_Y)))*FDR and n < len(pvals_Y):
        p_Y_cutoff=pvals_Y[n]
        n+=1
#    print n
#    print ((float(n)+1)/float(len(pvals_Y)))*FDR
#    print pvals_Y[0:100]
#    print pvals_Y[n]
#    print len(pvals_Y) 
        
    m=0
    while pvals_B[m] < ((float(m)+1.0)/float(len(pvals_B)))*FDR and m < len(pvals_B):
        p_B_cutoff=pvals_B[m]
        m+=1
    
    o=0
    while pvals_P[o] < ((float(o)+1.0)/float(len(pvals_P)))*FDR and o < len(pvals_P):
        p_P_cutoff=pvals_P[o]
        o+=1

    w=0
#    while pvals_Y_br[w] < ((float(w)+1.0)/float(len(pvals_Y_br)))*FDR and w < len(pvals_Y_br):
#        p_Y_br_cutoff=pvals_Y_br[w]
#        w+=1
    x=0
#    while pvals_Y_im[x] < ((float(x)+1.0)/float(len(pvals_Y_im)))*FDR and x < len(pvals_Y_im):
#        p_Y_im_cutoff=pvals_Y_im[x]
#        x+=1
    y=0
#    while y < len(pvals_Y_q) and pvals_Y_q[y] < ((float(y)+1.0)/float(len(pvals_Y_q)))*FDR:
#        p_Y_q_cutoff=pvals_Y_q[y]
#        y+=1
    z=0    
    if years != "2014":

        while pvals_B_br[z] < ((float(z)+1.0)/float(len(pvals_B_br)))*FDR and z < len(pvals_B_br):
            p_B_br_cutoff=pvals_B_br[z]
            z+=1    
    t=0
    while pvals_B_im[t] < ((float(t)+1.0)/float(len(pvals_B_im)))*FDR and t < len(pvals_B_im):
        p_B_im_cutoff=pvals_B_im[t]
        t+=1    
    k=0
    while pvals_B_q[k] < ((float(k)+1.0)/float(len(pvals_B_q)))*FDR and k < len(pvals_B_q):
        p_B_q_cutoff=pvals_B_q[k]
        k+=1
    h=0
    j=0
    if years == "both":
        while pvals_YCB_im[h] < ((float(h)+1.0)/float(len(pvals_YCB_im)))*FDR and h < len(pvals_YCB_im):
            p_YCB_im_cutoff=pvals_YCB_im[h]
            h+=1        
 
        while pvals_YCB_q[j] < ((float(j)+1.0)/float(len(pvals_YCB_q)))*FDR and j < len(pvals_YCB_q):
            p_YCB_q_cutoff=pvals_YCB_q[j]
            j+=1        
        
        
    stop = timeit.default_timer()
    runtime=stop-start
    
    paramfile=open(OUTDIR + "EL_Likelihoods_" + timestr + ".params.txt","w")
    paramfile.write("File ID: " + timestr + "\n")
    paramfile.write("Input: " + INPUT_FILE + "\n")
    paramfile.write("Years Included: " + years + "\n")
    paramfile.write("p_min: " + str(p_min1) + "\n")
    paramfile.write("p_max: " + str(p_max1) + "\n")
    paramfile.write("min_cov: " + str(min_cov) + "\n")
    paramfile.write("max_cov: " + str(max_cov) + "\n")
    paramfile.write("FDR: " + str(FDR) + "\n")
    paramfile.write("Number sites filtered: " + str(filtered) + "\n")
    paramfile.write("Total number of sites: " + str(sites) + "\n")
    paramfile.write("Program run time: " + str(runtime/60) + " minutes \n")
    paramfile.write("p_Y_cutoff = " + str(p_Y_cutoff) + " ; num_sig_tests = " + str(n) + " ; num_tests = " +str(len(pvals_Y)) + "\n")
    paramfile.write("p_P_cutoff = " + str(p_P_cutoff) + " ; num_sig_tests = " + str(o) + " ; num_tests = " +str(len(pvals_P)) + "\n")
    paramfile.write("p_B_cutoff = " + str(p_B_cutoff) + " ; num_sig_tests = " + str(m) + " ; num_tests = " +str(len(pvals_B)) + "\n")
    paramfile.write("p_Y_im_cutoff = " + str(p_Y_im_cutoff) + " ; num_sig_tests = " + str(x) + " ; num_tests = " +str(len(pvals_Y_im)) + "\n")
    paramfile.write("p_B_im_cutoff = " + str(p_B_im_cutoff) + " ; num_sig_tests = " + str(t) + " ; num_tests = " +str(len(pvals_B_im)) + "\n")
    paramfile.write("p_YCB_im_cutoff = " + str(p_YCB_im_cutoff) + " ; num_sig_tests = " + str(h) + " ; num_tests = " +str(len(pvals_YCB_im)) + "\n")
    paramfile.write("p_Y_q_cutoff = " + str(p_Y_q_cutoff) + " ; num_sig_tests = " + str(y) + " ; num_tests = " +str(len(pvals_Y_q)) + "\n")
    paramfile.write("p_B_q_cutoff = " + str(p_B_q_cutoff) + " ; num_sig_tests = " + str(k) + " ; num_tests = " +str(len(pvals_B_q)) + "\n")
    paramfile.write("p_YCB_q_cutoff = " + str(p_YCB_q_cutoff) + " ; num_sig_tests = " + str(j) + " ; num_tests = " +str(len(pvals_YCB_q)) + "\n")
    paramfile.write("p_Y_br_cutoff = " + str(p_Y_br_cutoff) + " ; num_sig_tests = " + str(w) + " ; num_tests = " +str(len(pvals_Y_br)) + "\n")
    paramfile.write("p_B_br_cutoff = " + str(p_B_br_cutoff) + " ; num_sig_tests = " + str(z) + " ; num_tests = " +str(len(pvals_B_br)) + "\n")
    paramfile.close()
    
    print "FDR = ",FDR
    print "p_Y_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_Y_cutoff, n, len(pvals_Y))         
    print "p_B_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_cutoff, m, len(pvals_B))
    print "p_P_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_P_cutoff, o, len(pvals_P))
    print "p_Y_im_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_Y_im_cutoff, x, len(pvals_Y_im))     
    print "p_YCB_im_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_YCB_im_cutoff, h, len(pvals_YCB_im))     
    print "p_B_im_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_im_cutoff, t, len(pvals_B_im))
    print "p_Y_br_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_Y_br_cutoff, w, len(pvals_Y_br) )        
    print "p_B_br_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_br_cutoff, z, len(pvals_B_br))
    print "p_Y_q_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_Y_q_cutoff, y, len(pvals_Y_q)    )     
    print "p_B_q_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_q_cutoff, k, len(pvals_B_q))
    print "p_YCB_q_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_YCB_q_cutoff, j, len(pvals_YCB_q)) 
    print "number of sites fixed for alt = ",fixed
    print "filtered sites = ", filtered
    print "number of sites =", sites
    print "program run time = ", runtime/60
    print "num_lines = ", i
    
    
    
#                        IME13maf=max(1-(ime13r/(ime13r+ime13a)),(ime13r/(ime13r+ime13a)))
#                        IML13maf=max((iml13r/(iml13r+iml13a)),1-(iml13r/(iml13r+iml13a)))
#                        IME14maf=max((ime14r/(ime14r+ime14a)),1-(ime14r/(ime14r+ime14a)))
#                        IML14maf=max((iml14r/(iml14r+iml14a)),1-(iml14r/(iml14r+iml14a)))
#
#
#                    if (ime13r/(ime13r+BRe13a))>0.5:
#                        IME13maf=(ime13r/(ime13r+ime13a))
#                        IML13maf=(iml13r/(iml13r+iml13a))
#                        IME14maf=(ime14r/(ime14r+ime14a))
#                        IML14maf=(iml14r/(iml14r+iml14a))
#                    else:
#                        IME13maf=1-(ime13r/(ime13r+ime13a))
#                        IML13maf=1-(iml13r/(iml13r+iml13a))
#                        IME14maf=1-(ime14r/(ime14r+ime14a))
#                        IML14maf=1-(iml14r/(iml14r+iml14a)) 

        