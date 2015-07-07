# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:37:06 2015

@author: patrick
"""

import csv
import math

EL_Likelihoods=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_Likelihoods.csv","wb")
like=csv.writer(EL_Likelihoods,delimiter=",",dialect='excel')
like.writerow(["scaff","pos","bre13r","bre13a","brl13r","brl13a","bre14r","bre14a","ime13r","ime13a","iml13r","iml13a","ime14r","ime14a","iml14r","iml14a","qe13r","qe13a","ql13r","ql13a","qe14r","qe14a","ql14r","ql14a"])

snp_null=[10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99]
snp_Y=[10.0,10.0,10.0,10.0,20.0,1.0,10.0,10.0,10.0,10.0,20.0,1.0,20.0,1.0,10.0,10.0,10.0,10.0,20.0,1.0,20.0,1.0,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99]
snp_P=[20.0,1.0,20.0,1.0,20.0,1.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,1.0,20.0,1.0,20.0,1.0,20.0,1.0,20.0,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99]
snp_B=[20.0,1.0,1.0,20.0,20.0,1.0,20.0,1.0,1.0,20.0,20.0,1.0,1.0,20.0,20.0,1.0,1.0,20.0,20.0,1.0,1.0,20.0,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99]
snp_YP=[10.0,10.0,10.0,10.0,20.0,1.0,1.0,20.0,1.0,20.0,20.0,1.0,20.0,1.0,20.0,1.0,20.0,1.0,1.0,20.0,1.0,20.0,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99]
snp_YB=[10.0,10.0,1.0,20.0,20.0,1.0,10.0,10.0,1.0,20.0,20.0,1.0,10.0,10.0,10.0,10.0,1.0,20.0,20.0,1.0,10.0,10.0,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99]
snp_PB=[20.0,1.0,17.0,3.0,20.0,1.0,10.0,10.0,10.0,5.0,10.0,10.0,10.0,5.0,1.0,20.0,5.0,20.0,1.0,20.0,5.0,20.0,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99]

with open("/Users/patrick/Google Drive/Research/EarlyLate/counts_filtered_6-13-15.txt","rb") as sites_file:   
    LL_S=0.0
    LL_Y=0.0
    LL_P=0.0
    LL_B=0.0
    LL_YP=0.0
    LL_YB=0.0
    LL_PB=0.0
    LL_C=0.0  
    for i,site in enumerate(sites_file):
        if i==0:
#            site=site.strip("\n")
#            site=site.split(" ")
#            site=site[1:]
#            scaff=site[0]
#            pos=float(site[1])
#            site=[float(a) for a in site[4:]]
            bre13r,bre13a,brl13r,brl13a,bre14r,bre14a,ime13r,ime13a,iml13r,iml13a,ime14r,ime14a,iml14r,iml14a,qe13r,qe13a,ql13r,ql13a,qe14r,qe14a,ql14r,ql14a,brer,brea,brlr,brla,imer,imea,imlr,imla,qer,qea,qlr,qla,im14t,im13t,br13t,br14t,q13t,q14t=snp_B
            xbre13=math.asin((bre13r/(bre13r+bre13a))**0.5)
            xbrl13=math.asin((brl13r/(brl13r+brl13a))**0.5)
            xbre14=math.asin((bre14r/(bre14r+bre14a))**0.5)
            xime13=math.asin((ime13r/(ime13r+ime13a))**0.5)
            ximl13=math.asin((iml13r/(iml13r+iml13a))**0.5)
            xime14=math.asin((ime14r/(ime14r+ime14a))**0.5)
            ximl14=math.asin((iml14r/(iml14r+iml14a))**0.5)
            xqe13=math.asin((qe13r/(qe13r+qe13a))**0.5)
            xql13=math.asin((ql13r/(ql13r+ql13a))**0.5)
            xqe14=math.asin((qe14r/(qe14r+qe14a))**0.5)
            xql14=math.asin((ql14r/(ql14r+ql14a))**0.5)
            
            vbre13=(0.02+(1/(bre13r+bre13a)))
            vbrl13=(0.02+(1/(brl13r+brl13a)))
            vbre14=(0.02+(1/(bre14r+bre14a)))
            vime13=(0.02+(1/(ime13r+ime13a)))
            viml13=(0.02+(1/(iml13r+iml13a)))
            vime14=(0.02+(1/(ime14r+ime14a)))
            viml14=(0.02+(1/(iml14r+iml14a)))
            vqe13=(0.02+(1/(qe13r+qe13a)))
            vql13=(0.02+(1/(ql13r+ql13a)))
            vqe14=(0.02+(1/(qe14r+qe14a)))
            vql14=(0.02 +(1/(ql14r+ql14a)))
            
            """Model 1 = Simplest Model.  All bulks share one mean."""                        
            x=[xbre13,xbrl13,xbre14,xime13,ximl13,xime14,ximl14,xqe13,xql13,xqe14,xql14]
            v=[vbre13,vbrl13,vbre14,vime13,viml13,vime14,viml14,vqe13,vql13,vqe14,vql14]
            w=[1/vbre13,1/vbrl13,1/vbre14,1/vime13,1/viml13,1/vime14,1/viml14,1/vqe13,1/vql13,1/vqe14,1/vql14]            
            W=sum(w)            
            uS=sum([xx*(ww/W) for xx,ww in zip(x,w)])
            ll_S=0.0
            for j,xx in enumerate(x):
                ll_S+=(-(xx-uS)**2)/(2*v[j])
            LL_S+=ll_S
            
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
            LL_Y+=ll_Y
            
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
            LL_P+=ll_P
            
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
            LL_B+=ll_B
            
            """Model 5 = Year and Pop effect.  No differences among bulks.""" 
            ll_YP=0.0
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
            
            for j,xx in enumerate(xbr13):
                ll_YP+=(-(xx-ubr13)**2)/(2*vbr13[j])
            for j,xx in enumerate(xim13):
                ll_YP+=(-(xx-uim13)**2)/(2*vim13[j])
            for j,xx in enumerate(xq13):
                ll_YP+=(-(xx-uq13)**2)/(2*vq13[j])
            for j,xx in enumerate(xbr14):
                ll_YP+=(-(xx-ubr14)**2)/(2*vbr14[j])
            for j,xx in enumerate(xim14):
                ll_YP+=(-(xx-uim14)**2)/(2*vim14[j])
            for j,xx in enumerate(xq14):
                ll_YP+=(-(xx-uq14)**2)/(2*vq14[j])
            LL_YP+=ll_YP
            
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
            
            for j,xx in enumerate(xe14):
                ll_YB+=(-(xx-ue14)**2)/(2*ve14[j])
            for j,xx in enumerate(xl14):
                ll_YB+=(-(xx-ul14)**2)/(2*vl14[j])
            for j,xx in enumerate(xe13):
                ll_YB+=(-(xx-ue13)**2)/(2*ve13[j])
            for j,xx in enumerate(xl13):
                ll_YB+=(-(xx-ul13)**2)/(2*vl13[j])
            LL_YB+=ll_YB

            """Model 7 = Pop and bulk effect.  No effect of year."""
            ll_PB=0.0
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
            
            for j,xx in enumerate(xbre):
                ll_PB+=(-(xx-ubre)**2)/(2*vbre[j])
            for j,xx in enumerate(xbrl):
                ll_PB+=(-(xx-ubrl)**2)/(2*vbrl[j])
            for j,xx in enumerate(xime):
                ll_PB+=(-(xx-uime)**2)/(2*vime[j])
            for j,xx in enumerate(ximl):
                ll_PB+=(-(xx-uiml)**2)/(2*viml[j])
            for j,xx in enumerate(xqe):
                ll_PB+=(-(xx-uqe)**2)/(2*vqe[j])
            for j,xx in enumerate(xql):
                ll_PB+=(-(xx-uql)**2)/(2*vql[j])
            LL_PB+=ll_PB  
            
            """Model 8 = Most complex model.  All bulks have distinct mean."""
            ll_C=0.0
            
            lrt_Y=chisqprob(-2*(ll_PB-ll_C),5)
            lrt_P=chisqprob(-2*(ll_YB-ll_C),7)
            lrt_B=chisqprob(-2*(ll_YP-ll_C),5)
            lrt_YP=chisqprob(-2*(ll_B-ll_C),9)
            lrt_YB=chisqprob(-2*(ll_P-ll_C),8)
            lrt_PB=chisqprob(-2*(ll_Y-ll_C),9)
            lrt_Yb=chisqprob(-2*(ll_B-ll_YB),2)
            lrt_Yp=chisqprob(-2*(ll_P-ll_YP),2)
            lrt_Py=chisqprob(-2*(ll_Y-ll_YP),2)
            like.writerow([bre13r,bre13a,brl13r,brl13a,bre14r,bre14a,ime13r,ime13a,iml13r,iml13a,ime14r,ime14a,iml14r,iml14a,qe13r,qe13a,ql13r,ql13a,qe14r,qe14a,ql14r,ql14a,ll_S,ll_Y,ll_P,ll_B,ll_YP,ll_YB,ll_PB,ll_C,lrt_Y,lrt_P,lrt_B,lrt_YP,lrt_YB,lrt_PB])   
            
    print "LL_S=",LL_S
    print "LL_Y=",LL_Y,lrt_Y
    print "LL_P=",LL_P,lrt_P
    print "LL_B=",LL_B,lrt_B
    print "LL_YP=",LL_YP,lrt_YP,lrt_Py,lrt_Yp
    print "LL_YB=",LL_YB,lrt_Yb
    print "LL_PB=",LL_PB,lrt_PB
    