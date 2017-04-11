# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 21:44:31 2016

@author: patrick
"""

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

INPUT_FILE="/Users/patrick/Documents/Research/EarlyLate/Current_Files/counts_massiveJan16_AF250.txt"
OUTDIR="/Users/patrick/Documents/Research/EarlyLate/"
#OUTDIR="/Volumes/avery/Research/EarlyLate/"

timestr = time.strftime("%Y%m%d-%H%M")

#filters
p_min1=0.05
p_max1=0.95
p_min=2*math.asin(p_min1**0.5)
p_max=2*math.asin(p_max1**0.5)
min_cov=25
max_cov=100
window_size=1000000

#Run Description
CustomMessage=""
   
FDR=0.1




out3=open(OUTDIR+"Fst_13" + timestr + ".csv","wb")
out3x=csv.writer(out3,delimiter=",",dialect='excel')
out3x.writerow(["pime13","piml13","pqe13","pql13","fstIM","fstQ","fstE","fstL"])
out4=open(OUTDIR+"Fst_14" + timestr + ".csv","wb")
out4x=csv.writer(out4,delimiter=",",dialect='excel')
out4x.writerow(["pime14","piml14","pqe14","pql14","fstIM","fstQ","fstE","fstL"])


#BS variance factor : Paired + Single
vb1=0.031943454869850194944 
vb2=0.032900188828741155911
vb3=0.027101898310498844652 

vi1=0.0141406715106648095404 #With new data the variance calculations were similar to old data except for 2014...Early var was similar to late var and vice versa.
vi2=0.0200071498211230672237 #Should run through calc with old data to confirm that there wasn't an error
vi3=0.0066349383482011761726 #Calculated with minimum depth of 25
vi4=0.0110139598993814532418

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
            
            if bre13r+bre13a <= min_cov or bre13r+bre13a >= max_cov:
                xbre13=-99.0
                vbre13=-99.0
            else:
                xbre13=2*math.asin((bre13r/(bre13r+bre13a))**0.5)
                vbre13=(vb1+(1/(bre13r+bre13a)))
                
            if brl13r+brl13a <= min_cov or brl13r+brl13a >= max_cov:
                xbrl13=-99.0
                vbrl13=-99.0
            else:
                xbrl13=2*math.asin((brl13r/(brl13r+brl13a))**0.5)
                vbrl13=(vb2+(1/(brl13r+brl13a)))
                
            if ime13r+ime13a <= min_cov or ime13r+ime13a >= max_cov:
                xime13=-99.0
                vime13=-99.0
            else:
                xime13=2*math.asin((ime13r/(ime13r+ime13a))**0.5)
                vime13=(vi1+(1/(ime13r+ime13a)))
                
            if iml13r+iml13a <= min_cov or iml13r+iml13a >= max_cov:
                ximl13=-99.0
                viml13=-99.0
            else:
                ximl13=2*math.asin((iml13r/(iml13r+iml13a))**0.5)
                viml13=(vi2+(1/(iml13r+iml13a)))
                
            if ime14r+ime14a <= min_cov or ime14r+ime14a >= max_cov:
                xime14=-99.0
                vime14=-99.0
            else:
                xime14=2*math.asin((ime14r/(ime14r+ime14a))**0.5)
                vime14=(vi3+(1/(ime14r+ime14a)))
                
            if iml14r+iml14a <= min_cov or iml14r+iml14a >= max_cov:
                ximl14=-99.0
                viml14=-99.0
            else:
                ximl14=2*math.asin((iml14r/(iml14r+iml14a))**0.5)
                viml14=(vi4+(1/(iml14r+iml14a)))
                
            if qe13r+qe13a <= min_cov or qe13r+qe13a >= max_cov:
                xqe13=-99.0
                vqe13=-99.0
            else:
                xqe13=2*math.asin((qe13r/(qe13r+qe13a))**0.5)
                vqe13=(vq1+(1/(qe13r+qe13a)))
                
                
            if ql13r+ql13a <= min_cov or ql13r+ql13a >= max_cov:
                xql13=-99.0
                vql13=-99.0
            else:
                xql13=2*math.asin((ql13r/(ql13r+ql13a))**0.5)
                vql13=(vq2+(1/(ql13r+ql13a)))
                
            if qe14r+qe14a <= min_cov or qe14r+qe14a >= max_cov:
                xqe14=-99.0
                vqe14=-99.0
            else:
                xqe14=2*math.asin((qe14r/(qe14r+qe14a))**0.5)
                vqe14=(vq3+(1/(qe14r+qe14a)))
            
            if ql14r+ql14a <= min_cov or ql14r+ql14a >= max_cov:
                xql14=-99.0
                vql14=-99.0
            else:
                xql14=2*math.asin((ql14r/(ql14r+ql14a))**0.5)
                vql14=(vq4 +(1/(ql14r+ql14a)))                
                

            if (all(k==0.0 for k in [xime13,ximl13,xime14,ximl14,xbre13,xbrl13,xqe13,xql13,xqe14,xql14]) or all(k==math.pi for k in [xime13,ximl13,xime14,ximl14,xbre13,xbrl13,xqe13,xql13,xqe14,xql14])):
                fixed+=1
                filtered+=1
            elif (all(k==-99.0 for k in [xime13,ximl13,xime14,ximl14,xbre13,xbrl13,xqe13,xql13,xqe14,xql14])):                
                filtered+=1
            
            else:
                sites+=1
                if sites%10000==0:
                    print sites
                    
                    
                """Model 5 = Year and Pop effect.  No differences among bulks.""" 
                    
                ll_Y_br13=0.0
                ll_Y_im13=0.0
                ll_Y_q13=0.0
                ll_Y_im14=0.0
                ll_Y_q14=0.0
                ll_Y_im=0.0
                ll_Y_q=0.0
                ll_Y_13=0.0
                ll_Y_14=0.0
                xbr13=[xbre13,xbrl13]
                vbr13=[vbre13,vbrl13]
                wbr13=[1/vbre13,1/vbrl13]
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
                uim13=sum([xx*(ww/sum(wim13)) for xx,ww in zip(xim13,wim13)])
                uim14=sum([xx*(ww/sum(wim14)) for xx,ww in zip(xim14,wim14)])
                uq13=sum([xx*(ww/sum(wq13)) for xx,ww in zip(xq13,wq13)])            
                uq14=sum([xx*(ww/sum(wq14)) for xx,ww in zip(xq14,wq14)]) 
                
                df_B_br13=1
                df_B_im13=1
                df_B_q13=1
                df_B_im14=1
                df_B_q14=1  
                df_B_im=2
                df_B_q=2 
                df_B_13=2
                df_B_14=2                 

                if ((xbre13 < p_min and xbrl13 < p_min) or (xbre13 > p_max and xbrl13 > p_max) or bre13r+bre13a <= min_cov or brl13r+brl13a <= min_cov or bre13r+bre13a >= max_cov or brl13r+brl13a >= max_cov):
                    df_B_br13-=1
                else:
                    for j,xx in enumerate(xbr13):                   
                        ll_Y_br13+=(-(xx-ubr13)**2)/(2*vbr13[j])
                        
                if ((xime13 < p_min and ximl13 < p_min) or (xime13 > p_max and ximl13 > p_max) or ime13r+ime13a <= min_cov or iml13r+iml13a <= min_cov or ime13r+ime13a >= max_cov or iml13r+iml13a >= max_cov):
                    df_B_im13-=1
                    df_B_im-=1
                    df_B_13-=1
                else:
                    for j,xx in enumerate(xim13):
                        ll_Y_im13+=(-(xx-uim13)**2)/(2*vim13[j])
                    ll_Y_im+=ll_Y_im13
                    ll_Y_13+=ll_Y_im13
                if ((xqe13 < p_min and xql13 < p_min) or (xqe13 > p_max and xql13 > p_max) or qe13r+qe13a <= min_cov or ql13r+ql13a <= min_cov or qe13r+qe13a >= max_cov or ql13r+ql13a >= max_cov):
                    df_B_q13-=1 
                    df_B_q-=1
                    df_B_13-=1
                else:
                    for j,xx in enumerate(xq13):
                        ll_Y_q13+=(-(xx-uq13)**2)/(2*vq13[j])
                    ll_Y_q+=ll_Y_q13
                    ll_Y_13+=ll_Y_q13
                    
                if ((xime14 < p_min and ximl14 < p_min) or (xime14 > p_max and ximl14 > p_max) or ime14r+ime14a <= min_cov or iml14r+iml14a <= min_cov or ime14r+ime14a >= max_cov or iml14r+iml14a >= max_cov):
                    df_B_im14-=1
                    df_B_im-=1
                    df_B_14-=1
                else:
                    for j,xx in enumerate(xim14):
                        ll_Y_im14+=(-(xx-uim14)**2)/(2*vim14[j])
                    ll_Y_im+=ll_Y_im14
                    ll_Y_14+=ll_Y_im14
                        
                if ((xqe14 < p_min and xql14 < p_min) or (xqe14 > p_max and xql14 > p_max) or qe14r+qe14a <= min_cov or ql14r+ql14a <= min_cov or qe14r+qe14a >= max_cov or ql14r+ql14a >= max_cov):
                    df_B_q14-=1
                    df_B_q-=1
                    df_B_14-=1
                else:
                    for j,xx in enumerate(xq14):                    
                        ll_Y_q14+=(-(xx-uq14)**2)/(2*vq14[j])   
                    ll_Y_q+=ll_Y_q14
                    ll_Y_14+=ll_Y_q14                                     
                
                    
                
                """Calculate Fsts"""                

                if df_B_im14 > 0 and df_B_q14>0:
                    pime14=math.sin(xime14/2)**2
                    piml14=math.sin(ximl14/2)**2
                    pim14=(pime14+piml14)/2
                    pqe14=math.sin(xqe14/2)**2
                    pql14=math.sin(xql14/2)**2
                    pq14=(pqe14+pql14)/2
                    pe14=(pime14+pqe14)/2
                    pl14=(piml14+pql14)/2
                    HsIM14=((2*(pime14)*(1-pime14))+(2*(piml14)*(1-piml14)))/2
                    HtIM14=(2*pim14*(1-pim14))            
                    HsQ14=((2*(pqe14)*(1-pqe14))+(2*(pql14)*(1-pql14)))/2
                    HtQ14=(2*pq14*(1-pq14))
                    HsE14=((2*(pqe14)*(1-pqe14))+(2*(pime14)*(1-pime14)))/2
                    HtE14=(2*pe14*(1-pe14))
                    HsL14=((2*(pql14)*(1-pql14))+(2*(piml14)*(1-piml14)))/2
                    HtL14=(2*pl14*(1-pl14))
                    try:
                        fstL14=(HtL14-HsL14)/HtL14
                    except ZeroDivisionError:
                        fstL14="-"
                    try:
                        fstE14=(HtE14-HsE14)/HtE14
                    except ZeroDivisionError:
                        fstE14="-"
                    try:
                        fstIM14=(HtIM14-HsIM14)/HtIM14
                    except ZeroDivisionError:
                        fstIM14="-"
                    try:
                        fstQ14=(HtQ14-HsQ14)/HtQ14
                    except ZeroDivisionError:
                        fstQ14="-"                                      
                    out4x.writerow([pime14,piml14,pqe14,pql14,fstIM14,fstQ14,fstE14,fstL14])   

                else:
                    p_B_im14="-"
                    lrt_B_im14="-"
                    pime14="-"
                    piml14="-"
                if df_B_im13 > 0 and df_B_q13>0:
                    pime13=math.sin(xime13/2)**2
                    piml13=math.sin(ximl13/2)**2
                    pim13=(pime13+piml13)/2
                    pqe13=math.sin(xqe13/2)**2
                    pql13=math.sin(xql13/2)**2
                    pq13=(pqe13+pql13)/2
                    pe13=(pime13+pqe13)/2
                    pl13=(piml13+pql13)/2
                    HsIM13=((2*(pime13)*(1-pime13))+(2*(piml13)*(1-piml13)))/2
                    HtIM13=(2*pim13*(1-pim13))                  
                    HsQ13=((2*(pqe13)*(1-pqe13))+(2*(pql13)*(1-pql13)))/2
                    HtQ13=(2*pq13*(1-pq13))
                    HsE13=((2*(pqe13)*(1-pqe13))+(2*(pime13)*(1-pime13)))/2
                    HtE13=(2*pe13*(1-pe13))
                    HsL13=((2*(pql13)*(1-pql13))+(2*(piml13)*(1-piml13)))/2
                    HtL13=(2*pl13*(1-pl13))
                    try:
                        fstL13=(HtL13-HsL13)/HtL13
                    except ZeroDivisionError:
                        fstL13="-"
                    try:
                        fstE13=(HtE13-HsE13)/HtE13
                    except ZeroDivisionError:
                        fstE13="-"
                    try:
                        fstIM13=(HtIM13-HsIM13)/HtIM13
                    except ZeroDivisionError:
                        fstIM13="-"
                    try:
                        fstQ13=(HtQ13-HsQ13)/HtQ13
                    except ZeroDivisionError:
                        fstQ13="-"
                    out3x.writerow([pime13,piml13,pqe13,pql13,fstIM13,fstQ13,fstE13,fstL13])   
                else:
                    p_B_im13="-"
                    lrt_B_im13="-"
                    pime13="-"
                    piml13="-"
