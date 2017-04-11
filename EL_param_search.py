# -*- coding: utf - 8 -*-
"""
Created on Wed May  4 17:33:01 2016

@author: patrick monnahan
"""

# -*- coding: utf - 8 -*-
"""
Created on Mon Mar 28 10:50:52 2016

Purpose: Simulate data using a range of values for the fraction of sites under selection (f0) as well as effect sizes (ES). 
Record the proportion of sites that exceed certain lrt values (10 and 15).  Use this to find a set of parameters that
most closely match the observed distributions of lrts for the test of model 1 versus model 2. (consistent effects).

Use this to guide input of f0 and ES for EL_parametric_sim.py
"""

from scipy.stats import *
import csv
import math
import time
import timeit
from random import *
import numpy as np
import sys

def simulate(effect_size, fraction, InFile, Outdir):
    start = timeit.default_timer()
    
    INPUT_FILE = InFile
    OUTDIR = Outdir
    
    timestr = time.strftime("%Y%m%d-%H%M")
    
    # filters
    p_min1 = 0.05
    p_max1 = 0.95
    p_min = 2 * math.asin(p_min1**0.5)
    p_max = 2 * math.asin(p_max1**0.5)
    min_cov = 25
    max_cov = 100
    
    # Simulation parameters

    intensity = 1.755  # intensity of selection. Calculated from truncated normal in which 10% of individuals exceed truncation point
    
    timestr = time.strftime("%Y%m%d-%H%M")
        
    # Run Description
    CustomMessage = ""
       
    FDR = 0.1
    
    paramfile1 = open(OUTDIR + "EL_Likelihoods_ParamSim_m1f0Search_"+ timestr + ".txt", "w")

    # BS variance factor : Paired + Single
    vb1 = 0.03225929
    vb2 = 0.03552918
    vb3 = 0.01646089

    vi1 = 0.01414067
    vi2 = 0.02000715
    vi3 = 0.00663494
    vi4 = 0.01101396

    vq1 = 0.01524366
    vq2 = 0.02004908
    vq3 = 0.00833159
    vq4 = 0.01379415

    fo = [fraction]
    es = [effect_size]
    for if0, f0 in enumerate(fo):
        for ES in np.arange(es[if0], es[if0] + 0.11, 0.01):
            im_15_avg = 0.0
            im_10_avg = 0.0
            q_15_avg = 0.0
            q_10_avg = 0.0
            i13_10_avg = 0.0
            i13_15_avg = 0.0
            i14_10_avg = 0.0
            i14_15_avg = 0.0
            for rep in range(0, 5):
                with open(INPUT_FILE, "rb") as sites_file:   
                    sites = 0
                    fixed = 0
                    filtered = 0
                    dist = 0
                    pos = 0
                    wind_sites = 0
                    wind_df = 0
                    pvals_Y_im = []
                    pvals_Y_q = []
                    pvals_B_br13 = []
                    pvals_B_im13 = []
                    pvals_B_q13 = []
                    pvals_B_im14 = []
                    pvals_B_q14 = []
                    pvals_B_q_M12 = []
                    pvals_B_q_M13 = []
                    pvals_B_q_M23 = []
                    pvals_B_im_M12 = []
                    pvals_B_im_M13 = []
                    pvals_B_im_M23 = []
                    pvals_B_13_M12 = []
                    pvals_B_13_M13 = []
                    pvals_B_13_M23 = []
                    pvals_B_14_M12 = []
                    pvals_B_14_M13 = []
                    pvals_B_14_M23 = []
                    pvals_P = []
                    pvals_Y_windows = []
                    pvals_P_windows = []
                    pvals_B_windows = []
                    exceptions = 0
                    im_xs = 0.0
                    im_x = 0.0
                    q_xs = 0.0
                    q_x = 0.0
                    xs13 = 0.0
                    x13 = 0.0
                    xs14 = 0.0
                    x14 = 0.0
                    im_12_count = 0.0
                    q_12_count = 0.0
                    i13_12_count = 0.0
                    i14_12_count = 0.0
                    im_12_count1 = 0.0
                    q_12_count1 = 0.0
                    i13_12_count1 = 0.0
                    i14_12_count1 = 0.0
                    im_13_count = 0.0
                    q_13_count = 0.0
                    i13_13_count = 0.0
                    i14_13_count = 0.0
                    im_23_count = 0.0
                    q_23_count = 0.0
                    i13_23_count = 0.0
                    i14_23_count = 0.0
                    
                    for i, site in enumerate(sites_file):
                        if i>0 and i % 100 == 0:  # SAMPLE ONLY 1 OUT OF EVERY 100 SITES
                            rx = np.random.uniform()
                            if rx>f0:
                                effect_size = 0.0
                            else:
                                rx1 = np.random.uniform()
                                if rx1>0.5:
                                    effect_size = ES
                                else:
                                    effect_size = -ES
                                    
                            site = site.strip("\n")
                            site = site.split("\t")
                            scaff = site[0]
                            scaff = scaff
                            pos = float(site[1])
                            site = [float(a) for a in site[4:]]
                            bre13r, bre13a, brl13r, brl13a, ime14r, ime14a, iml14r, iml14a, ime13r, ime13a, iml13r, iml13a, qe14r, qe14a, ql14r, ql14a, qe13r, qe13a, ql13r, ql13a = site
                            
                            bre14r = 0.0
                            bre14a = 0.0
                            xbre14 = -99.0
                            vbre14 = -99.0
                            
                            if bre13r + bre13a <= min_cov or bre13r + bre13a >= max_cov or brl13r + brl13a <= min_cov or brl13r + brl13a >= max_cov:
                                pbr13 = -99.0
                                pbre13 = -99.0
                                pbrl13 = -99.0
                                xbre13 = -99.0
                                vbre13 = -99.0
                                xbrl13 = -99.0
                                vbrl13 = -99.0
                                                
                            else:                
                                pbr13 = ((bre13r / (bre13r + bre13a)) + (brl13r / (brl13r + brl13a))) / 2
                                pbre13 = pbr13 - (0.5 * intensity * effect_size * pbr13 * (1 - pbr13))
                                pbrl13 = pbr13 + (0.5 * intensity * effect_size * pbr13 * (1 - pbr13))
                                if pbre13>1.0:
                                    pbre13 = 1.0
                                elif pbre13<0.0:
                                    pbre13 = 0.0
                                if pbrl13>1.0:
                                    pbrl13 = 1.0
                                elif pbrl13<0.0:
                                    pbrl13 = 0.0
                                vbre13 = (vb1 + (1 / (bre13r + bre13a)))
                                xbre13 = 2 * math.asin(pbre13**0.5) + normalvariate(0, vbre13**0.5)
                                vbrl13 = (vb2 + (1 / (brl13r + brl13a)))
                                xbrl13 = 2 * math.asin(pbrl13**0.5) + normalvariate(0, vbrl13**0.5)
                                if xbre13>math.pi:
                                    xbre13 = math.pi
                                elif xbre13<0.0:
                                    xbre13 = 0.0
                                if xbrl13>math.pi:
                                    xbrl13 = math.pi
                                elif xbrl13<0.0:
                                    xbrl13 = 0.0
                                    
                            if ime13r + ime13a <= min_cov or ime13r + ime13a >= max_cov or iml13r + iml13a <= min_cov or iml13r + iml13a >= max_cov:
                                pim13 = -99.0
                                pime13 = -99.0
                                piml13 = -99.0
                                xime13 = -99.0
                                vime13 = -99.0
                                ximl13 = -99.0
                                viml13 = -99.0                                           
                            else:                
                                pim13 = ((ime13r / (ime13r + ime13a)) + (iml13r / (iml13r + iml13a))) / 2
                                pime13 = pim13 - (0.5 * intensity * effect_size * pim13 * (1 - pim13))
                                piml13 = pim13 + (0.5 * intensity * effect_size * pim13 * (1 - pim13))
                                if pime13>1.0:
                                    pime13 = 1.0
                                elif pime13<0.0:
                                    pime13 = 0.0
                                if piml13>1.0:
                                    piml13 = 1.0
                                elif piml13<0.0:
                                    piml13 = 0.0

                                vime13 = (vi1 + (1 / (ime13r + ime13a)))
                                dev1 = normalvariate(0, vime13**0.5)
                                xime13 = 2 * math.asin(pime13**0.5) + dev1
                                viml13 = (vi2 + (1 / (iml13r + iml13a)))
                                dev2 = normalvariate(0, viml13**0.5)
                                ximl13 = 2 * math.asin(piml13**0.5) + dev2
                                
                                if xime13>math.pi:
                                    xime13 = math.pi
                                elif xime13<0.0:
                                    xime13 = 0.0
                                if ximl13>math.pi:
                                    ximl13 = math.pi
                                elif ximl13<0.0:
                                    ximl13 = 0.0
                                
                            if ime14r + ime14a <= min_cov or ime14r + ime14a >= max_cov or iml14r + iml14a <= min_cov or iml14r + iml14a >= max_cov:
                                pim14 = -99.0
                                pime14 = -99.0
                                piml14 = -99.0
                                xime14 = -99.0
                                vime14 = -99.0
                                ximl14 = -99.0
                                viml14 = -99.0
                            else:               
                                pim14 = ((ime14r / (ime14r + ime14a)) + (iml14r / (iml14r + iml14a))) / 2
                                pime14 = pim14 - (0.5 * intensity * effect_size * pim14 * (1 - pim14))
                                piml14 = pim14 + (0.5 * intensity * effect_size * pim14 * (1 - pim14))
                                if pime14>1.0:
                                    pime14 = 1.0
                                elif pime14<0.0:
                                    pime14 = 0.0
                                if piml14>1.0:
                                    piml14 = 1.0
                                elif piml14<0.0:
                                    piml14 = 0.0
                                vime14 = (vi3 + (1 / (ime14r + ime14a)))
                                xime14 = 2 * math.asin(pime14**0.5) + normalvariate(0, vime14**0.5)
                                viml14 = (vi4 + (1 / (iml14r + iml14a)))
                                ximl14 = 2 * math.asin(piml14**0.5) + normalvariate(0, viml14**0.5)
                                
                                if xime14>math.pi:
                                    xime14 = math.pi
                                elif xime14<0.0:
                                    xime14 = 0.0
                                if ximl14>math.pi:
                                    ximl14 = math.pi
                                elif ximl14<0.0:
                                    ximl14 = 0.0
                                
                            if qe13r + qe13a <= min_cov or qe13r + qe13a >= max_cov or ql13r + ql13a <= min_cov or ql13r + ql13a >= max_cov:
                                pq13 = -99.0
                                pqe13 = -99.0
                                pql13 = -99.0
                                xqe13 = -99.0
                                vqe13 = -99.0
                                xql13 = -99.0
                                vql13 = -99.0
                            else:               
                                pq13 = ((qe13r / (qe13r + qe13a)) + (ql13r / (ql13r + ql13a))) / 2
                                pqe13 = pq13 - (0.5 * intensity * effect_size * pq13 * (1 - pq13))
                                pql13 = pq13 + (0.5 * intensity * effect_size * pq13 * (1 - pq13))
                                if pqe13>1.0:
                                    pqe13 = 1.0
                                elif pqe13<0.0:
                                    pqe13 = 0.0
                                if pql13>1.0:
                                    pql13 = 1.0
                                elif pql13<0.0:
                                    pql13 = 0.0
                                vqe13 = (vq1 + (1 / (qe13r + qe13a)))
                                xqe13 = 2 * math.asin(pqe13**0.5) + normalvariate(0, vqe13**0.5)
                                vql13 = (vq2 + (1 / (ql13r + ql13a)))
                                xql13 = 2 * math.asin(pql13**0.5) + normalvariate(0, vql13**0.5)
                                
                                if xqe13>math.pi:
                                    xqe13 = math.pi
                                elif xqe13<0.0:
                                    xqe13 = 0.0
                                if xql13>math.pi:
                                    xql13 = math.pi
                                elif xql13<0.0:
                                    xql13 = 0.0
                                
                            if qe14r + qe14a <= min_cov or qe14r + qe14a >= max_cov or ql14r + ql14a <= min_cov or ql14r + ql14a >= max_cov:
                                pq14 = -99.0
                                pqe14 = -99.0
                                pql14 = -99.0
                                xqe14 = -99.0
                                vqe14 = -99.0
                                xql14 = -99.0
                                vql14 = -99.0
                            else:         
                                pq14 = ((qe14r / (qe14r + qe14a)) + (ql14r / (ql14r + ql14a))) / 2
                                pqe14 = pq14 - (0.5 * intensity * effect_size * pq14 * (1 - pq14))
                                pql14 = pq14 + (0.5 * intensity * effect_size * pq14 * (1 - pq14))
                                if pqe14>1.0:
                                    pqe14 = 1.0
                                elif pqe14<0.0:
                                    pqe14 = 0.0
                                if pql14>1.0:
                                    pql14 = 1.0
                                elif pql14<0.0:
                                    pql14 = 0.0
                                vqe14 = (vq3 + (1 / (qe14r + qe14a)))
                                xqe14 = 2 * math.asin(pqe14**0.5) + normalvariate(0, vqe14**0.5)
                                vql14 = (vq4 + (1 / (ql14r + ql14a)))
                                xql14 = 2 * math.asin(pql14**0.5) + normalvariate(0, vql14**0.5)
                                
                                if xqe14>math.pi:
                                    xqe14 = math.pi
                                elif xqe14<0.0:
                                    xqe14 = 0.0
                                if xql14>math.pi:
                                    xql14 = math.pi
                                elif xql14<0.0:
                                    xql14 = 0.0
                
                            if (all(k == 0.0 for k in [xime13, ximl13, xime14, ximl14, xbre13, xbrl13, xqe13, xql13, xqe14, xql14]) or all(k == math.pi for k in [xime13, ximl13, xime14, ximl14, xbre13, xbrl13, xqe13, xql13, xqe14, xql14])):
                                fixed += 1
                                filtered += 1
                            elif (all(k==-99.0 for k in [xime13, ximl13, xime14, ximl14, xbre13, xbrl13, xqe13, xql13, xqe14, xql14])):                
                                filtered += 1
                            
                            else:
                                sites += 1
                                    
                                """Model 5 = Year and Pop effect.  No differences among bulks.""" 
                                    
                                ll_Y_br13 = 0.0
                                ll_Y_im13 = 0.0
                                ll_Y_q13 = 0.0
                                ll_Y_im14 = 0.0
                                ll_Y_q14 = 0.0
                                ll_Y_im = 0.0
                                ll_Y_q = 0.0
                                ll_Y_13 = 0.0
                                ll_Y_14 = 0.0
                                xbr13 = [xbre13, xbrl13]
                                vbr13 = [vbre13, vbrl13]
                                wbr13 = [1 / vbre13, 1 / vbrl13]
                                xim13 = [xime13, ximl13]
                                xim14 = [xime14, ximl14]
                                vim13 = [vime13, viml13]
                                vim14 = [vime14, viml14]
                                wim13 = [1 / vime13, 1 / viml13]
                                wim14 = [1 / vime14, 1 / viml14]
                                xq13 = [xqe13, xql13]
                                xq14 = [xqe14, xql14]
                                vq13 = [vqe13, vql13]
                                vq14 = [vqe14, vql14]
                                wq13 = [1 / vqe13, 1 / vql13]
                                wq14 = [1 / vqe14, 1 / vql14]
                                ubr13 = sum([xx * (ww / sum(wbr13)) for xx, ww in zip(xbr13, wbr13)])
                                uim13 = sum([xx * (ww / sum(wim13)) for xx, ww in zip(xim13, wim13)])
                                uim14 = sum([xx * (ww / sum(wim14)) for xx, ww in zip(xim14, wim14)])
                                uq13 = sum([xx * (ww / sum(wq13)) for xx, ww in zip(xq13, wq13)])            
                                uq14 = sum([xx * (ww / sum(wq14)) for xx, ww in zip(xq14, wq14)])
                                
                                uim13_yy = uim13
                                uim14_yy = uim14
                                uq13_yy = uq13
                                
                                df_B_br13 = 1
                                df_B_im13 = 1
                                df_B_q13 = 1
                                df_B_im14 = 1
                                df_B_q14 = 1  
                                df_B_im = 2
                                df_B_q = 2 
                                df_B_13 = 2
                                df_B_14 = 2                 
                
                                if ((xbre13 < p_min and xbrl13 < p_min) or (xbre13 > p_max and xbrl13 > p_max) or bre13r + bre13a <= min_cov or brl13r + brl13a <= min_cov or bre13r + bre13a >= max_cov or brl13r + brl13a >= max_cov):
                                    df_B_br13-=1
                                else:
                                    for j, xx in enumerate(xbr13):                   
                                        ll_Y_br13+=(-(xx - ubr13)**2) / (2 * vbr13[j])
                                        
                                if ((xime13 < p_min and ximl13 < p_min) or (xime13 > p_max and ximl13 > p_max) or ime13r + ime13a <= min_cov or iml13r + iml13a <= min_cov or ime13r + ime13a >= max_cov or iml13r + iml13a >= max_cov):
                                    df_B_im13-=1
                                    df_B_im-=1
                                    df_B_13-=1
                                else:
                                    for j, xx in enumerate(xim13):
                                        ll_Y_im13+=(-(xx - uim13)**2) / (2 * vim13[j])
                                    ll_Y_im += ll_Y_im13
                                    ll_Y_13 += ll_Y_im13
                                if ((xqe13 < p_min and xql13 < p_min) or (xqe13 > p_max and xql13 > p_max) or qe13r + qe13a <= min_cov or ql13r + ql13a <= min_cov or qe13r + qe13a >= max_cov or ql13r + ql13a >= max_cov):
                                    df_B_q13-=1 
                                    df_B_q-=1
                                    df_B_13-=1
                                else:
                                    for j, xx in enumerate(xq13):
                                        ll_Y_q13+=(-(xx - uq13)**2) / (2 * vq13[j])
                                    ll_Y_q += ll_Y_q13
                                    ll_Y_13 += ll_Y_q13
                                    
                                if ((xime14 < p_min and ximl14 < p_min) or (xime14 > p_max and ximl14 > p_max) or ime14r + ime14a <= min_cov or iml14r + iml14a <= min_cov or ime14r + ime14a >= max_cov or iml14r + iml14a >= max_cov):
                                    df_B_im14-=1
                                    df_B_im-=1
                                    df_B_14-=1
                                else:
                                    for j, xx in enumerate(xim14):
                                        ll_Y_im14+=(-(xx - uim14)**2) / (2 * vim14[j])
                                    ll_Y_im += ll_Y_im14
                                    ll_Y_14 += ll_Y_im14
                                        
                                if ((xqe14 < p_min and xql14 < p_min) or (xqe14 > p_max and xql14 > p_max) or qe14r + qe14a <= min_cov or ql14r + ql14a <= min_cov or qe14r + qe14a >= max_cov or ql14r + ql14a >= max_cov):
                                    df_B_q14-=1
                                    df_B_q-=1
                                    df_B_14-=1
                                else:
                                    for j, xx in enumerate(xq14):                    
                                        ll_Y_q14+=(-(xx - uq14)**2) / (2 * vq14[j])   
                                    ll_Y_q += ll_Y_q14
                                    ll_Y_14 += ll_Y_q14
                                
                                """Model 6 = Year and bulk effect.  No differences across pops."""
                                ll_YB = 0.0
                                xe13 = [xbre13, xime13, xqe13]
                                xe14 = [xime14, xqe14]
                                xl13 = [xbrl13, ximl13, xql13]
                                xl14 = [ximl14, xql14]
                                ve13 = [vbre13, vime13, vqe13]
                                ve14 = [vime14, vqe14]
                                vl13 = [vbrl13, viml13, vql13]
                                vl14 = [viml14, vql14]
                                we13 = [1 / vbre13, 1 / vime13, 1 / vqe13]
                                we14 = [1 / vime14, 1 / vqe14]
                                wl13 = [1 / vbrl13, 1 / viml13, 1 / vql13]
                                wl14 = [1 / viml14, 1 / vql14]
                                ue13 = sum([xx * (ww / sum(we13)) for xx, ww in zip(xe13, we13)])
                                ue14 = sum([xx * (ww / sum(we14)) for xx, ww in zip(xe14, we14)])
                                ul13 = sum([xx * (ww / sum(wl13)) for xx, ww in zip(xl13, wl13)])
                                ul14 = sum([xx * (ww / sum(wl14)) for xx, ww in zip(xl14, wl14)])
                                
                                df_P = 6
                                
                                if (any(k ==-99.0 for k in [xime14, ximl14, xqe14, xql14])):
                                    df_P -= 2
                                else:
                                    for j, xx in enumerate(xe14):
                                        ll_YB+=(-(xx - ue14)**2) / (2 * ve14[j])
                                    for j, xx in enumerate(xl14):
                                        ll_YB+=(-(xx - ul14)**2) / (2 * vl14[j])
                                if (any(k ==-99.0 for k in [xime13, ximl13, xqe13, xql13, xbre13, xbrl13])):
                                    df_P -= 4
                                else:                                    
                                    for j, xx in enumerate(xe13):
                                        ll_YB+=(-(xx - ue13)**2) / (2 * ve13[j])
                                    for j, xx in enumerate(xl13):
                                        ll_YB+=(-(xx - ul13)**2) / (2 * vl13[j])
                                
                                """Model 6b = Year effect and consistent effect of bulk"""
                                ll_YCB_im = 0.0
                                ll_YCB_q = 0.0
                                uim13 = -((-viml13 * xime13 - vime13 * ximl13 - vime14 * ximl13 - viml14 * ximl13 + viml13 * xime14 - viml13 * ximl14) / (vime13 + viml13 + vime14 + viml14))
                                uim14 = -((viml14 * xime13 - viml14 * ximl13 - viml14 * xime14 - vime13 * ximl14 - viml13 * ximl14 - vime14 * ximl14) / (vime13 + viml13 + vime14 + viml14))                
                                uq13 = -((-vql13 * xqe13 - vqe13 * xql13 - vqe14 * xql13 - vql14 * xql13 + vql13 * xqe14 - vql13 * xql14) / (vqe13 + vql13 + vqe14 + vql14))
                                uq14 = -((vql14 * xqe13 - vql14 * xql13 - vql14 * xqe14 - vqe13 * xql14 - vql13 * xql14 - vqe14 * xql14) / (vqe13 + vql13 + vqe14 + vql14))                                
                                aim = -((-vime14 * xime13 - viml14 * xime13 + vime14 * ximl13 + viml14 * ximl13 - vime13 * xime14 - viml13 * xime14 + vime13 * ximl14 + viml13 * ximl14) / (vime13 + viml13 + vime14 + viml14))
                                aq = -((-vqe14 * xqe13 - vql14 * xqe13 + vqe14 * xql13 + vql14 * xql13 - vqe13 * xqe14 - vql13 * xqe14 + vqe13 * xql14 + vql13 * xql14) / (vqe13 + vql13 + vqe14 + vql14))
                                aim_xx = aim
                                uim13_xx = uim13
                                uim14_xx = uim14
                                df_YCB_im = 1
                                df_YCB_q = 1                    
                                        
                                if (xime13 < p_min and ximl13 < p_min) or (xime13 > p_max and ximl13 > p_max) or (xime14 < p_min and ximl14 < p_min) or (xime14 > p_max and ximl14 > p_max) or (any(k <= min_cov for k in [ime13r + ime13a, iml13r + iml13a, ime14r + ime14a, iml14r + iml14a])) or (any(k >= max_cov for k in [ime13r + ime13a, iml13r + iml13a, ime14r + ime14a, iml14r + iml14a])):
                                    df_YCB_im-=1
                                else:
                                    ll_YCB_im+=(-(xime13 - (uim13 + aim))**2) / (2 * vime13)
                                    ll_YCB_im+=(-(xime14 - (uim14 + aim))**2) / (2 * vime14)
                                    ll_YCB_im+=(-(ximl13 - (uim13))**2) / (2 * viml13)
                                    ll_YCB_im+=(-(ximl14 - (uim14))**2) / (2 * viml14)
                                        
                                if (xqe13 < p_min and xql13 < p_min) or (xqe13 > p_max and xql13 > p_max) or (xqe14 < p_min and xql14 < p_min) or (xqe14 > p_max and xql14 > p_max) or (any(k <= min_cov for k in [qe13r + qe13a, ql13r + ql13a, qe14r + qe14a, ql14r + ql14a])) or (any(k >= max_cov for k in [qe13r + qe13a, ql13r + ql13a, qe14r + qe14a, ql14r + ql14a])):
                                    df_YCB_q+=(-1) 
                                else:
                                    ll_YCB_q+=(-(xqe13 - (uq13 + aq))**2) / (2 * vqe13)
                                    ll_YCB_q+=(-(xqe14 - (uq14 + aq))**2) / (2 * vqe14)
                                    ll_YCB_q+=(-(xql13 - (uq13))**2) / (2 * vql13)
                                    ll_YCB_q+=(-(xql14 - (uq14))**2) / (2 * vql14)
                                    
                                """Model 6c = Consistent effect of bulk across populations within a year."""
                                ll_PCB_13 = 0.0
                                ll_PCB_14 = 0.0
                                uim13 = -((-viml13 * xime13 - vime13 * ximl13 - vqe13 * ximl13 - vql13 * ximl13 + viml13 * xqe13 - viml13 * xql13) / (vime13 + viml13 + vqe13 + vql13))
                                uq13 = -((vql13 * xime13 - vql13 * ximl13 - vql13 * xqe13 - vime13 * xql13 - viml13 * xql13 - vqe13 * xql13) / (vime13 + viml13 + vqe13 + vql13))                
                                uim14 = -((-viml14 * xime14 - vime14 * ximl14 - vqe14 * ximl14 - vql14 * ximl14 + viml14 * xqe14 - viml14 * xql14) / (vime14 + viml14 + vqe14 + vql14))
                                uq14 = -((vql14 * xime14 - vql14 * ximl14 - vql14 * xqe14 - vime14 * xql14 - viml14 * xql14 - vqe14 * xql14) / (vime14 + viml14 + vqe14 + vql14))                                
                                a13 = -((-vqe13 * xime13 - vql13 * xime13 + vqe13 * ximl13 + vql13 * ximl13 - vime13 * xqe13 - viml13 * xqe13 + vime13 * xql13 + viml13 * xql13) / (vime13 + viml13 + vqe13 + vql13))
                                a14 = -((-vqe14 * xime14 - vql14 * xime14 + vqe14 * ximl14 + vql14 * ximl14 - vime14 * xqe14 - viml14 * xqe14 + vime14 * xql14 + viml14 * xql14) / (vime14 + viml14 + vqe14 + vql14))
                                
                                df_PCB_13 = 1
                                df_PCB_14 = 1                    
                                        
                                if (xime13 < p_min and ximl13 < p_min) or (xime13 > p_max and ximl13 > p_max) or (xqe13 < p_min and xql13 < p_min) or (xqe13 > p_max and xql13 > p_max) or (any(k <= min_cov for k in [ime13r + ime13a, iml13r + iml13a, qe13r + qe13a, ql13r + ql13a])) or (any(k >= max_cov for k in [ime13r + ime13a, iml13r + iml13a, qe13r + qe13a, ql13r + ql13a])):
                                    df_PCB_13-=1
                                else:
                                    ll_PCB_13+=(-(xime13 - (uim13 + a13))**2) / (2 * vime13)
                                    ll_PCB_13+=(-(xqe13 - (uq13 + a13))**2) / (2 * vqe13)
                                    ll_PCB_13+=(-(ximl13 - (uim13))**2) / (2 * viml13)
                                    ll_PCB_13+=(-(xql13 - (uq13))**2) / (2 * vql13)
                                        
                                if (xime14 < p_min and ximl14 < p_min) or (xime14 > p_max and ximl14 > p_max) or (xqe14 < p_min and xql14 < p_min) or (xqe14 > p_max and xql14 > p_max) or (any(k <= min_cov for k in [ime14r + ime14a, iml14r + iml14a, qe14r + qe14a, ql14r + ql14a])) or (any(k >= max_cov for k in [ime14r + ime14a, iml14r + iml14a, qe14r + qe14a, ql14r + ql14a])):
                                    df_PCB_14+=(-1) 
                                else:
                                    ll_PCB_14+=(-(xime14 - (uim14 + a14))**2) / (2 * vime14)
                                    ll_PCB_14+=(-(xqe14 - (uq14 + a14))**2) / (2 * vqe14)
                                    ll_PCB_14+=(-(ximl14 - (uim14))**2) / (2 * viml14)
                                    ll_PCB_14+=(-(xql14 - (uq14))**2) / (2 * vql14)
                
                    
                                """Model 7 = Pop and bulk effect.  No effect of year."""
                                ll_B_im = 0.0
                                ll_B_q = 0.0
                                xime = [xime13, xime14]
                                ximl = [ximl13, ximl14]
                                vime = [vime13, vime14]
                                viml = [viml13, viml14]
                                wime = [1 / vime13, 1 / vime14]
                                wiml = [1 / viml13, 1 / viml14]
                                uime = sum([xx * (ww / sum(wime)) for xx, ww in zip(xime, wime)])
                                uiml = sum([xx * (ww / sum(wiml)) for xx, ww in zip(ximl, wiml)])
                                xqe = [xqe13, xqe14]
                                xql = [xql13, xql14]
                                vqe = [vqe13, vqe14]
                                vql = [vql13, vql14]
                                wqe = [1 / vqe13, 1 / vqe14]
                                wql = [1 / vql13, 1 / vql14]
                                uqe = sum([xx * (ww / sum(wqe)) for xx, ww in zip(xqe, wqe)])
                                uql = sum([xx * (ww / sum(wql)) for xx, ww in zip(xql, wql)])
                                                               
                                df_Y_im = 2
                                df_Y_q = 2
                                                           
                                if ((xime13 < p_min and ximl13 < p_min and xime14 < p_min and ximl14 < p_min) or (xime13 > p_max and ximl13 > p_max and xime14 > p_max and ximl14 > p_max) or (any(k <= min_cov for k in [ime13r + ime13a, iml13r + iml13a, ime14r + ime14a, iml14r + iml14a])) or (any(k >= max_cov for k in [ime13r + ime13a, iml13r + iml13a, ime14r + ime14a, iml14r + iml14a]))):
                                    df_Y_im-=2
                                else:
                                    for j, xx in enumerate(xime):
                                        ll_B_im+=(-(xx - uime)**2) / (2 * vime[j])
                                    for j, xx in enumerate(ximl):
                                        ll_B_im+=(-(xx - uiml)**2) / (2 * viml[j])
                                
                                if ((xqe13 < p_min and xql13 < p_min and xqe14 < p_min and xql14 < p_min) or (xqe13 > p_max and xql13 > p_max and xqe14 > p_max and xql14 > p_max) or (any(k <= min_cov for k in [qe13r + qe13a, ql13r + ql13a, qe14r + qe14a, ql14r + ql14a])) or (any(k >= max_cov for k in [qe13r + qe13a, ql13r + ql13a, qe14r + qe14a, ql14r + ql14a]))):
                                    df_Y_q-=2
                                else:
                                    for j, xx in enumerate(xqe):
                                        ll_B_q+=(-(xx - uqe)**2) / (2 * vqe[j])
                                    for j, xx in enumerate(xql):
                                        ll_B_q+=(-(xx - uql)**2) / (2 * vql[j])
                                
                                """Model 8 = Most complex model.  All bulks have distinct mean."""
                                ll_C = 0.0
                                
                                """tests for each population individually"""                
                                if df_B_br13 > 0:
                                    lrt_B_br13 = -2 * ll_Y_br13
                                    p_B_br13 = chisqprob(lrt_B_br13, df_B_br13)
                                    pvals_B_br13.append(p_B_br13)
                                else:
                                    p_B_br13 = "-"
                                    lrt_B_br13 = "-"
                                if df_Y_im > 0:
                                    lrt_Y_im = -2 * ll_B_im
                                    p_Y_im = chisqprob(lrt_Y_im, df_Y_im)
                                    pvals_Y_im.append(p_Y_im)
                                else:
                                    p_Y_im = "-"
                                    lrt_Y_im = "-"
                                if df_B_im14 > 0:
                                    lrt_B_im14 = -2 * ll_Y_im14
                                    p_B_im14 = chisqprob(lrt_B_im14, df_B_im14)
                                    pvals_B_im14.append(p_B_im14)
                                else:
                                    p_B_im14 = "-"
                                    lrt_B_im14 = "-"
                                if df_B_im13 > 0:
                                    lrt_B_im13 = -2 * ll_Y_im13
                                    p_B_im13 = chisqprob(lrt_B_im13, df_B_im13)
                                    pvals_B_im13.append(p_B_im13)          
                                else:
                                    p_B_im13 = "-"
                                    lrt_B_im13 = "-"
                                if df_Y_q > 0:
                                    lrt_Y_q = -2 * ll_B_q
                                    p_Y_q = chisqprob(lrt_Y_q, df_Y_q)
                                    pvals_Y_q.append(p_Y_q)
                                else:
                                    p_Y_q = "-"
                                    lrt_Y_q = "-"
                                if df_B_q14 > 0:
                                    lrt_B_q14 = -2 * ll_Y_q14
                                    p_B_q14 = chisqprob(lrt_B_q14, df_B_q14)
                                    pvals_B_q14.append(p_B_q14)
                                else:
                                    p_B_q14 = "-"
                                    lrt_B_q14 = "-"
                                if df_B_q13 > 0:
                                    lrt_B_q13 = -2 * ll_Y_q13
                                    p_B_q13 = chisqprob(lrt_B_q13, df_B_q13)
                                    pvals_B_q13.append(p_B_q13)
                                else:
                                    p_B_q13 = "-"
                                    lrt_B_q13 = "-"
                                if df_B_im > 1:  
                                    lrt_B_im_M12 = -2 * (ll_Y_im - ll_YCB_im)  # Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)
                                    lrt_B_im_M13 = -2 * (ll_Y_im)
                                    lrt_B_im_M23 = -2 * (ll_YCB_im)
                                    p_B_im_M12 = chisqprob(lrt_B_im_M12, 1)
                                    p_B_im_M13 = chisqprob(lrt_B_im_M13, 2)
                                    p_B_im_M23 = chisqprob(lrt_B_im_M23, 1)
                                    pvals_B_im_M12.append(p_B_im_M12)
                                    pvals_B_im_M13.append(p_B_im_M13)
                                    pvals_B_im_M23.append(p_B_im_M23)
                                    im_xs += lrt_B_im_M12**2.0
                                    im_x += lrt_B_im_M12
                                else:
                                    lrt_B_im_M12 = "-"  # Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)
                                    lrt_B_im_M13 = "-"
                                    lrt_B_im_M23 = "-"
                                    p_B_im_M12 = "-"
                                    p_B_im_M13 = "-"
                                    p_B_im_M23 = "-"
                                if df_B_q > 1:  
                                    lrt_B_q_M12 = -2 * (ll_Y_q - ll_YCB_q)  # Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)
                                    lrt_B_q_M13 = -2 * (ll_Y_q)
                                    lrt_B_q_M23 = -2 * (ll_YCB_q)
                                    p_B_q_M12 = chisqprob(lrt_B_q_M12, 1)
                                    p_B_q_M13 = chisqprob(lrt_B_q_M13, 2)
                                    p_B_q_M23 = chisqprob(lrt_B_q_M23, 1)
                                    pvals_B_q_M12.append(p_B_q_M12)
                                    pvals_B_q_M13.append(p_B_q_M13)
                                    pvals_B_q_M23.append(p_B_q_M23)
                                    q_xs += lrt_B_q_M12**2.0
                                    q_x += lrt_B_q_M12
                                else:
                                    lrt_B_q_M12 = "-"  # Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)
                                    lrt_B_q_M13 = "-"
                                    lrt_B_q_M23 = "-"
                                    p_B_q_M12 = "-"
                                    p_B_q_M13 = "-"
                                    p_B_q_M23 = "-"
                                if df_B_13 > 1:  
                                    lrt_B_13_M12 = -2 * (ll_Y_13 - ll_PCB_13)  # Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)
                                    lrt_B_13_M13 = -2 * (ll_Y_13)
                                    lrt_B_13_M23 = -2 * (ll_PCB_13)
                                    p_B_13_M12 = chisqprob(lrt_B_13_M12, 1)
                                    p_B_13_M13 = chisqprob(lrt_B_13_M13, 2)
                                    p_B_13_M23 = chisqprob(lrt_B_13_M23, 1)
                                    pvals_B_13_M12.append(p_B_13_M12)
                                    pvals_B_13_M13.append(p_B_13_M13)
                                    pvals_B_13_M23.append(p_B_13_M23)
                                    xs13 += lrt_B_13_M12**2.0
                                    x13 += lrt_B_13_M12
                                else:
                                    lrt_B_13_M12 = "-"  # Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)
                                    lrt_B_13_M13 = "-"
                                    lrt_B_13_M23 = "-"
                                    p_B_13_M12 = "-"
                                    p_B_13_M13 = "-"
                                    p_B_13_M23 = "-"
                                if df_B_14 > 1:  
                                    lrt_B_14_M12 = -2 * (ll_Y_14 - ll_PCB_14)  # Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)                 
                                    lrt_B_14_M13 = -2 * (ll_Y_14)
                                    lrt_B_14_M23 = -2 * (ll_PCB_14)
                                    p_B_14_M12 = chisqprob(lrt_B_14_M12, 1)
                                    p_B_14_M13 = chisqprob(lrt_B_14_M13, 2)
                                    p_B_14_M23 = chisqprob(lrt_B_14_M23, 1)
                                    pvals_B_14_M12.append(p_B_14_M12)
                                    pvals_B_14_M13.append(p_B_14_M13)
                                    pvals_B_14_M23.append(p_B_14_M23)
                                    xs14 += lrt_B_14_M12**2.0
                                    x14 += lrt_B_14_M12
                                else:
                                    lrt_B_14_M12 = "-"  # Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)
                                    lrt_B_14_M13 = "-"
                                    lrt_B_14_M23 = "-"
                                    p_B_14_M12 = "-"
                                    p_B_14_M13 = "-"
                                    p_B_14_M23 = "-"
                                    
                                
                                
                                """test for all pops together"""                
                                if df_P > 0:    
                                    lrt_P = -2 * (ll_YB - ll_C)
                                    p_P = chisqprob(lrt_P, df_P)
                                    pvals_P.append(p_P)
                                else:
                                    lrt_P = "-"
                                    p_P = "-"
                                    
                                if lrt_B_im_M12 == "-":
                                    pass
                                elif lrt_B_im_M12>=10.0:
                                    im_12_count1 += 1.0
                                    if lrt_B_im_M12>=15.0:
                                        im_12_count += 1.0
                                if lrt_B_q_M12 == "-":
                                    pass
                                elif lrt_B_q_M12>=10.0:
                                    q_12_count1 += 1.0
                                    if lrt_B_q_M12>=15.0:
                                        q_12_count += 1.0
                                if lrt_B_13_M12 == "-":
                                    pass
                                elif lrt_B_13_M12>=10.0:
                                    i13_12_count1 += 1.0
                                    if lrt_B_13_M12>=15.0:
                                        i13_12_count += 1.0
                                if lrt_B_14_M12 == "-":
                                    pass
                                elif lrt_B_14_M12>=10.0:
                                    i14_12_count1 += 1.0
                                    if lrt_B_14_M12>=15.0:
                                        i14_12_count += 1.0
                                
                                if lrt_B_im_M13 == "-":
                                    pass
                                elif lrt_B_im_M13>=15.0:
                                    im_13_count += 1.0
                                if lrt_B_q_M13 == "-":
                                    pass
                                elif lrt_B_q_M13>=15.0:
                                    q_13_count += 1.0
                                if lrt_B_13_M13 == "-":
                                    pass
                                elif lrt_B_13_M13>=15.0:
                                    i13_13_count += 1.0
                                if lrt_B_14_M13 == "-":
                                    pass
                                elif lrt_B_14_M13>=15.0:
                                    i14_13_count += 1.0
                                
                                if lrt_B_im_M23 == "-":
                                    pass
                                elif lrt_B_im_M23>=15.0:
                                    im_23_count += 1.0
                                if lrt_B_q_M23 == "-":
                                    pass
                                elif lrt_B_q_M23>=15.0:
                                    q_23_count += 1.0
                                if lrt_B_13_M23 == "-":
                                    pass
                                elif lrt_B_13_M23>=15.0:
                                    i13_23_count += 1.0
                                if lrt_B_14_M23 == "-":
                                    pass
                                elif lrt_B_14_M23>=15.0:
                                    i14_23_count += 1.0
                                    
                    """Determine p - value tresholds given a specified FDR"""
                    pvals_Y_im.sort()
                    pvals_Y_q.sort()
                    pvals_B_br13.sort()
                    pvals_B_im13.sort()
                    pvals_B_q13.sort()
                    pvals_B_im14.sort()
                    pvals_B_q14.sort()
                    pvals_B_q_M12.sort()
                    pvals_B_q_M13.sort()
                    pvals_B_q_M23.sort()
                    pvals_B_im_M12.sort()
                    pvals_B_im_M13.sort()
                    pvals_B_im_M23.sort()
                    pvals_B_13_M12.sort()
                    pvals_B_13_M13.sort()
                    pvals_B_13_M23.sort()
                    pvals_B_14_M12.sort()
                    pvals_B_14_M13.sort()
                    pvals_B_14_M23.sort()
                    pvals_P.sort()
                    pvals_Y_windows.sort()
                    pvals_P_windows.sort()
                    pvals_B_windows.sort()
                    
                    p_P_cutoff = 0.0
                    p_Y_im_cutoff = 0.0
                    p_Y_q_cutoff = 0.0
                    p_B_br13_cutoff = 0.0
                    p_B_im13_cutoff = 0.0
                    p_B_q13_cutoff = 0.0
                    p_B_im14_cutoff = 0.0
                    p_B_q14_cutoff = 0.0
                    p_B_q_M12_cutoff = 0.0
                    p_B_q_M13_cutoff = 0.0
                    p_B_q_M23_cutoff = 0.0
                    p_B_im_M12_cutoff = 0.0
                    p_B_im_M13_cutoff = 0.0
                    p_B_im_M23_cutoff = 0.0
                    p_B_13_M12_cutoff = 0.0
                    p_B_13_M13_cutoff = 0.0
                    p_B_13_M23_cutoff = 0.0
                    p_B_14_M12_cutoff = 0.0
                    p_B_14_M13_cutoff = 0.0
                    p_B_14_M23_cutoff = 0.0
                            
                    o = 0
                    for m, p in enumerate(pvals_P):
                        if p < ((float(m) + 1.0) / float(len(pvals_P))) * FDR:
                            p_P_cutoff = p
                            o = m + 1
                    x = 0
                    for m, p in enumerate(pvals_Y_im):
                        if p < ((float(m) + 1.0) / float(len(pvals_Y_im))) * FDR:
                            p_Y_im_cutoff = p
                            x = m + 1
                    y = 0   
                    for m, p in enumerate(pvals_Y_q):
                        if p < ((float(m) + 1.0) / float(len(pvals_Y_q))) * FDR:
                            p_Y_q_cutoff = p
                            y = m + 1
                    t = 0    
                    for m, p in enumerate(pvals_B_im13):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_im13))) * FDR:
                            p_B_im13_cutoff = p
                            t = m + 1
                    qq = 0    
                    for m, p in enumerate(pvals_B_q13):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_q13))) * FDR:
                            p_B_q13_cutoff = p
                            qq = m + 1
                    z = 0      
                    for m, p in enumerate(pvals_B_br13):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_br13))) * FDR:
                            p_B_br13_cutoff = p
                            z = m + 1 
                    
                    tt = 0    
                    for m, p in enumerate(pvals_B_im14):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_im14))) * FDR:
                            p_B_im14_cutoff = p
                            tt = m + 1
                    www = 0    
                    for m, p in enumerate(pvals_B_q14):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_q14))) * FDR:
                            p_B_q14_cutoff = p
                            www = m + 1
                                
                    h = 0           
                    for m, p in enumerate(pvals_B_im_M12):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_im_M12))) * FDR:
                            p_B_im_M12_cutoff = p
                            h = m + 1
                            
                    hh = 0           
                    for m, p in enumerate(pvals_B_im_M13):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_im_M13))) * FDR:
                            p_B_im_M13_cutoff = p
                            hh = m + 1
                            
                    hhh = 0           
                    for m, p in enumerate(pvals_B_im_M23):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_im_M23))) * FDR:
                            p_B_im_M23_cutoff = p
                            hhh = m + 1
                            
                    j = 0 
                    for m, p in enumerate(pvals_B_q_M12):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_q_M12))) * FDR:
                            p_B_q_M12_cutoff = p
                            j = m + 1
                            
                    jj = 0 
                    for m, p in enumerate(pvals_B_q_M13):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_q_M13))) * FDR:
                            p_B_q_M13_cutoff = p
                            jj = m + 1
                            
                    jjj = 0 
                    for m, p in enumerate(pvals_B_q_M23):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_q_M23))) * FDR:
                            p_B_q_M23_cutoff = p
                            jjj = m + 1
                            
                    k = 0           
                    for m, p in enumerate(pvals_B_13_M12):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_13_M12))) * FDR:
                            p_B_13_M12_cutoff = p
                            k = m + 1
                            
                    kk = 0           
                    for m, p in enumerate(pvals_B_13_M13):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_13_M13))) * FDR:
                            p_B_13_M13_cutoff = p
                            kk = m + 1
                            
                    kkk = 0           
                    for m, p in enumerate(pvals_B_13_M23):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_13_M23))) * FDR:
                            p_B_13_M23_cutoff = p
                            kkk = m + 1
                            
                    l = 0 
                    for m, p in enumerate(pvals_B_14_M12):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_14_M12))) * FDR:
                            p_B_14_M12_cutoff = p
                            l = m + 1
                            
                    ll = 0 
                    for m, p in enumerate(pvals_B_14_M13):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_14_M13))) * FDR:
                            p_B_14_M13_cutoff = p
                            ll = m + 1
                            
                    lll = 0 
                    for m, p in enumerate(pvals_B_14_M23):
                        if p < ((float(m) + 1.0) / float(len(pvals_B_14_M23))) * FDR:
                            p_B_14_M23_cutoff = p
                            lll = m + 1
                
                        
                    stop = timeit.default_timer()
                    runtime = stop - start
                    im_15_avg += im_12_count / len(pvals_B_im_M12)
                    im_10_avg += im_12_count1 / len(pvals_B_im_M12)
                    q_15_avg += q_12_count / len(pvals_B_q_M12)
                    q_10_avg += q_12_count1 / len(pvals_B_q_M12)
                    i13_15_avg += i13_12_count / len(pvals_B_13_M12)
                    i13_10_avg += i13_12_count1 / len(pvals_B_13_M12)
                    i14_15_avg += i14_12_count / len(pvals_B_14_M12)
                    i14_10_avg += i14_12_count1 / len(pvals_B_14_M12)
                    
                    print "FDR = ", FDR
                    print "p_P_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_P_cutoff, x, len(pvals_P))     
                    print "p_Y_im_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_Y_im_cutoff, x, len(pvals_Y_im))     
                    print "p_Y_q_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_Y_q_cutoff, y, len(pvals_Y_q))     
                    print "p_B_im13_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_im13_cutoff, t, len(pvals_B_im13))
                    print "p_B_q13_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_q13_cutoff, qq, len(pvals_B_q13))
                    print "p_B_br13_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_br13_cutoff, z, len(pvals_B_br13))
                    print "p_B_im14_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_im14_cutoff, tt, len(pvals_B_im14))
                    print "p_B_q14_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_q14_cutoff, www, len(pvals_B_q14))
                    print "p_B_im_M12_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_im_M12_cutoff, h, len(pvals_B_im_M12))
                    print "p_B_im_M13_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_im_M13_cutoff, hh, len(pvals_B_im_M13))
                    print "p_B_im_M23_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_im_M23_cutoff, hhh, len(pvals_B_im_M23))
                    print "p_B_q_M12_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_q_M12_cutoff, j, len(pvals_B_q_M12))
                    print "p_B_q_M13_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_q_M13_cutoff, jj, len(pvals_B_q_M13))
                    print "p_B_q_M23_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_q_M23_cutoff, jjj, len(pvals_B_q_M23))     
                    print "p_B_13_M12_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_13_M12_cutoff, k, len(pvals_B_13_M12))
                    print "p_B_13_M13_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_13_M13_cutoff, kk, len(pvals_B_13_M13))
                    print "p_B_13_M23_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_13_M23_cutoff, kkk, len(pvals_B_13_M23))     
                    print "p_B_14_M12_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_14_M12_cutoff, l, len(pvals_B_14_M12))
                    print "p_B_14_M13_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_14_M13_cutoff, ll, len(pvals_B_14_M13))
                    print "p_B_14_M23_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_B_14_M23_cutoff, lll, len(pvals_B_14_M23))     
                    print "number of sites fixed for alt = ", fixed
                    print "filtered sites = ", filtered
                    print "number of sites =", sites
                    print "program run time = ", runtime / 60
                    print "num_lines = ", i
            paramfile1.write(str(ES) + "\t" + str(f0) + "\t" + str(im_10_avg / 5.0) + "\t" + str(im_15_avg / 5.0) + "\t" + str(q_10_avg / 5.0) + "\t" + str(q_15_avg / 5.0) + "\t" + str(i13_10_avg / 5.0) + "\t" + str(i13_15_avg / 5.0) + "\t" + str(i14_10_avg / 5.0) + "\t" + str(i14_15_avg / 5.0) + "\n")
        
            # All values above are divided by 5.0 because we are performing 5 replicates for each value of f0 and c.

    paramfile1.close()      

if __name__ == '__main__':
    import sys
    effect_size = float(sys.argv[1])
    fraction = float(sys.argv[2])
    infile = sys.argv[3]
    outdir = sys.argv[4]
    print sys.argv[1], sys.argv[2]
    simulate(effect_size, fraction, infile, outdir)
