# - * - coding: utf - 8 - *  -
"""
Created on Wed Jun 17 10:37:06 2015

Important Note:  Variable names may seem inconsistent or confusing.  For the likelihood, ll_Y_im is the likelihood in which bulk is constrained (i.e. year is allowed to vary, which will be compared to a model in which both yeah AND bulk car vary...hence the Y in lrt_Y)
                I changed the letters between the likelihoods and the likelihood ratios because if ll_Y_br (year allowed to vary) is compared to a model in which year and bulk are unconstrained (ll_C = 0.0), this provides a test of the effect of the bulk.

@author: patrick
  """
from scipy.stats import *
import csv
import math
import time
import timeit

start = timeit.default_timer()

INPUT_FILE = "/Users/patrick/Documents/Research/EarlyLate/Current_Files/counts_massiveJan16_AF250.txt"
# OUTDIR = " / Volumes / TOSHIBA EXT / EarlyLate / "
OUTDIR = "/Users/patrick/Documents/Research/EarlyLate/"
# OUTDIR = " / Volumes / avery / Research / EarlyLate / "

timestr = time.strftime("%Y%m%d-%H%M")

# Indicator variables used to specify which output files to produce
file1 = 0
file2 = 0
file3 = 0
file6 = 0

# filters
p_min1 = 0.05
p_max1 = 0.95
p_min = 2 * math.asin(p_min1 ** 0.5)
p_max = 2 * math.asin(p_max1 ** 0.5)
min_cov = 25
max_cov = 100

# Run Description
CustomMessage = ""

FDR = 0.1

# Open files and print headers
if file1 != 0:
    EL_Likelihoods = open(OUTDIR + "EL_Likelihoods_IndTests_" + timestr + ".csv", "wb")
    like = csv.writer(EL_Likelihoods, delimiter=",", dialect='excel')
    like.writerow(["scaff", "pos", "lrt_B", "p_B", "pop"])
if file2 != 0:
    out2 = open(OUTDIR + "EL_Likelihoods_Interactions_" + timestr + ".csv", "wb")
    out2x = csv.writer(out2, delimiter=",", dialect='excel')
    out2x.writerow(["scaff", "pos", "lrt_B_M12", "p_B_M12", "lrt_B_M13", "p_B_M13", "lrt_B_M23", "p_B_M23", "pop"])
if file3 != 0:
    out3 = open(OUTDIR + "EL_Likelihoods_Slim13_RealData_" + timestr + ".csv", "wb")
    out3x = csv.writer(out3, delimiter=",", dialect='excel')
    out3x.writerow(["lrt_B_im13", "p_B_im13", "lrt_B_q13", "p_B_q13", "lrt_B_13_M12", "p_B_13_M12", "lrt_B_13_M13", "p_B_13_M13", "lrt_B_13_M23", "p_B_13_M23", "pime13", "piml13", "pqe13", "pql13"])
    out4 = open(OUTDIR + "EL_Likelihoods_Slim14_RealData_" + timestr + ".csv", "wb")
    out4x = csv.writer(out4, delimiter=",", dialect='excel')
    out4x.writerow(["lrt_B_im14", "p_B_im14", "lrt_B_q14", "p_B_q14", "lrt_B_14_M12", "p_B_14_M12", "lrt_B_14_M13", "p_B_14_M13", "lrt_B_14_M23", "p_B_14_M23", "pime14", "piml14", "pqe14", "pql14"])
    out5 = open(OUTDIR + "EL_Likelihoods_SlimIM_RealData_" + timestr + ".csv", "wb")
    out5x = csv.writer(out5, delimiter=",", dialect='excel')
    out5x.writerow(["scaff", "pos", "lrt_B_im14", "p_B_im14", "lrt_B_im13", "p_B_im13", "lrt_B_im_M12", "p_B_im_M12", "lrt_B_im_M13", "p_B_im_M13", "lrt_B_im_M23", "p_B_im_M23", "pime13", "piml13", "pime14", "piml14"])
    out6 = open(OUTDIR + "EL_Likelihoods_SlimQ_RealData_" + timestr + ".csv", "wb")
    out6x = csv.writer(out6, delimiter=",", dialect='excel')
    out6x.writerow(["scaff", "pos", "lrt_B_q14", "p_B_q14", "lrt_B_q13", "p_B_q13", "lrt_B_q_M12", "p_B_q_M12", "lrt_B_q_M13", "p_B_q_M13", "lrt_B_q_M23", "p_B_q_M23", "pqe13", "pql13", "pqe14", "pql14"])

# bulk variance (v) terms
vb1 = 0.03225929  # All v - terms re - estimated prior to resubmission to Genetics
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

# Begin loop over input file
with open(INPUT_FILE, "rb") as sites_file:
    sites = 0
    fixed = 0
    filtered = 0
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
    polim = 0
    polq = 0
    polbr = 0
    polimq = 0
    polimbr = 0
    polqbr = 0
    polimqbr = 0
    polim1 = 0
    polq1 = 0
    polimq1 = 0
    for i, site in enumerate(sites_file):
        if i > 0:
            # Parse information from line of input file
            site = site.strip("\n")
            site = site.split("\t")
            scaff = site[0]
            scaff = scaff
            pos = float(site[1])
            site = [float(a) for a in site[4:]]
            bre13r, bre13a, brl13r, brl13a, ime14r, ime14a, iml14r, iml14a, ime13r, ime13a, iml13r, iml13a, qe14r, qe14a, ql14r, ql14a, qe13r, qe13a, ql13r, ql13a = site

            # Transform allele frequency and set to missing value (-99.0) if read depth does not pass filter
            if bre13r + bre13a <= min_cov or bre13r + bre13a >= max_cov:
                xbre13 = - 99.0
                vbre13 = - 99.0
            else:
                xbre13 = 2 * math.asin((bre13r / (bre13r + bre13a))**0.5)
                vbre13 = (vb1 + (1 / (bre13r + bre13a)))

            if brl13r + brl13a <= min_cov or brl13r + brl13a >= max_cov:
                xbrl13 = -99.0
                vbrl13 = -99.0
            else:
                xbrl13 = 2 * math.asin((brl13r / (brl13r + brl13a))**0.5)
                vbrl13 = (vb2 + (1 / (brl13r + brl13a)))

            if ime13r + ime13a <= min_cov or ime13r + ime13a >= max_cov:
                xime13 = -99.0
                vime13 = -99.0
            else:
                xime13 = 2 * math.asin((ime13r / (ime13r + ime13a))**0.5)
                vime13 = (vi1 + (1 / (ime13r + ime13a)))

            if iml13r + iml13a <= min_cov or iml13r + iml13a >= max_cov:
                ximl13 = -99.0
                viml13 = -99.0
            else:
                ximl13 = 2 * math.asin((iml13r / (iml13r + iml13a))**0.5)
                viml13 = (vi2 + (1 / (iml13r + iml13a)))

            if ime14r + ime14a <= min_cov or ime14r + ime14a >= max_cov:
                xime14 = -99.0
                vime14 = -99.0
            else:
                xime14 = 2 * math.asin((ime14r / (ime14r + ime14a))**0.5)
                vime14 = (vi3 + (1 / (ime14r + ime14a)))

            if iml14r + iml14a <= min_cov or iml14r + iml14a >= max_cov:
                ximl14 = -99.0
                viml14 = -99.0
            else:
                ximl14 = 2 * math.asin((iml14r / (iml14r + iml14a))**0.5)
                viml14 = (vi4 + (1 / (iml14r + iml14a)))

            if qe13r + qe13a <= min_cov or qe13r + qe13a >= max_cov:
                xqe13 = -99.0
                vqe13 = -99.0
            else:
                xqe13 = 2 * math.asin((qe13r / (qe13r + qe13a))**0.5)
                vqe13 = (vq1 + (1 / (qe13r + qe13a)))

            if ql13r + ql13a <= min_cov or ql13r + ql13a >= max_cov:
                xql13 = -99.0
                vql13 = -99.0
            else:
                xql13 = 2 * math.asin((ql13r / (ql13r + ql13a))**0.5)
                vql13 = (vq2 + (1 / (ql13r + ql13a)))

            if qe14r + qe14a <= min_cov or qe14r + qe14a >= max_cov:
                xqe14 = -99.0
                vqe14 = -99.0
            else:
                xqe14 = 2 * math.asin((qe14r / (qe14r + qe14a))**0.5)
                vqe14 = (vq3 + (1 / (qe14r + qe14a)))

            if ql14r + ql14a <= min_cov or ql14r + ql14a >= max_cov:
                xql14 = -99.0
                vql14 = -99.0
            else:
                xql14 = 2 * math.asin((ql14r / (ql14r + ql14a))**0.5)
                vql14 = (vq4 + (1 / (ql14r + ql14a)))

            # Filter out sites that are fixed for reference or alternative allele across all samples
            if (all(k == 0.0 for k in [xime13, ximl13, xime14, ximl14, xbre13, xbrl13, xqe13, xql13, xqe14, xql14]) or all(k == math.pi for k in [xime13, ximl13, xime14, ximl14, xbre13, xbrl13, xqe13, xql13, xqe14, xql14])):
                fixed += 1
                filtered += 1
            # Filter out sites that have insufficient read depth for all samples
            elif all(k == -99.0 for k in [xime13, ximl13, xime14, ximl14, xbre13, xbrl13, xqe13, xql13, xqe14, xql14]):
                filtered += 1
            # Perform analysis for site
            else:
                sites += 1
                if sites % 10000 == 0:
                    print(sites)

                """Calculate likelihood for each population / year assuming no differences among bulks"""
                # EXPLICIT COMMENTING PERFORMED ONLY FOR THIS TEST.  REMAINDER OF TESTS FOLLOW SIMILAR FORMAT
                ll_Y_br13 = 0.0  # Likelihood for no difference between Early and Late flowering plants in Browder Ridge in 2013
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

                # Calculate null mean value as weighted average across samples
                ubr13 = sum([xx * (ww / sum(wbr13)) for xx, ww in zip(xbr13, wbr13)])
                uim13 = sum([xx * (ww / sum(wim13)) for xx, ww in zip(xim13, wim13)])
                uim14 = sum([xx * (ww / sum(wim14)) for xx, ww in zip(xim14, wim14)])
                uq13 = sum([xx * (ww / sum(wq13)) for xx, ww in zip(xq13, wq13)])
                uq14 = sum([xx * (ww / sum(wq14)) for xx, ww in zip(xq14, wq14)])

                # Set degrees of freedom
                df_B_br13 = 1
                df_B_im13 = 1
                df_B_q13 = 1
                df_B_im14 = 1
                df_B_q14 = 1
                df_B_im = 2
                df_B_q = 2
                df_B_13 = 2
                df_B_14 = 2

                # These filters subtract degrees of freedom if samples do not pass filter.  P-values and LRT's are not calculated for sites with 0 degrees of freedom
                if all(k < p_min for k in [xbre13, xbrl13]) or all(k > p_max for k in [xbre13, xbrl13]) or any(k == -99.0 for k in [xbre13, xbrl13]):
                    df_B_br13 -= 1
                # Calculate likelihood if site passes filter
                else:
                    for j, xx in enumerate(xbr13):
                        ll_Y_br13 += (-(xx - ubr13)**2) / (2 * vbr13[j])

                if all(k < p_min for k in [xime13, ximl13]) or all(k > p_max for k in [xime13, ximl13]) or any(k == -99.0 for k in [xime13, ximl13]):
                    df_B_im13 -= 1
                    df_B_im -= 1
                    df_B_13 -= 1
                else:
                    for j, xx in enumerate(xim13):
                        ll_Y_im13 += (-(xx - uim13)**2) / (2 * vim13[j])
                    ll_Y_im += ll_Y_im13
                    ll_Y_13 += ll_Y_im13

                if all(k < p_min for k in [xqe13, xql13]) or all(k > p_max for k in [xqe13, xql13]) or any(k == -99.0 for k in [xqe13, xql13]):
                    df_B_q13 -= 1
                    df_B_q -= 1
                    df_B_13 -= 1
                else:
                    for j, xx in enumerate(xq13):
                        ll_Y_q13 += (-(xx - uq13)**2) / (2 * vq13[j])
                    ll_Y_q += ll_Y_q13
                    ll_Y_13 += ll_Y_q13

                if all(k < p_min for k in [xime14, ximl14]) or all(k > p_max for k in [xime14, ximl14]) or any(k == -99.0 for k in [xime14, ximl14]):
                    df_B_im14 -= 1
                    df_B_im -= 1
                    df_B_14 -= 1
                else:
                    for j, xx in enumerate(xim14):
                        ll_Y_im14 += (-(xx - uim14)**2) / (2 * vim14[j])
                    ll_Y_im += ll_Y_im14
                    ll_Y_14 += ll_Y_im14

                if all(k < p_min for k in [xqe14, xql14]) or all(k > p_max for k in [xqe14, xql14]) or any(k == -99.0 for k in [xqe14, xql14]):
                    df_B_q14 -= 1
                    df_B_q -= 1
                    df_B_14 -= 1
                else:
                    for j, xx in enumerate(xq14):
                        ll_Y_q14 += (-(xx - uq14)**2) / (2 * vq14[j])
                    ll_Y_q += ll_Y_q14
                    ll_Y_14 += ll_Y_q14

                """Year and bulk effect.  No differences across pops."""
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

                if (any(k == - 99.0 for k in [xime14, ximl14, xqe14, xql14])):
                    df_P -= 2
                else:
                    for j, xx in enumerate(xe14):
                        ll_YB += (-(xx - ue14)**2) / (2 * ve14[j])
                    for j, xx in enumerate(xl14):
                        ll_YB += (-(xx - ul14)**2) / (2 * vl14[j])
                if (any(k == - 99.0 for k in [xime13, ximl13, xqe13, xql13, xbre13, xbrl13])):
                    df_P -= 4
                else:
                    for j, xx in enumerate(xe13):
                        ll_YB += (-(xx - ue13)**2) / (2 * ve13[j])
                    for j, xx in enumerate(xl13):
                        ll_YB += (-(xx - ul13)**2) / (2 * vl13[j])

                """Year effect and consistent effect of bulk"""
                ll_YCB_im = 0.0
                ll_YCB_q = 0.0
                uim13 = -((-viml13 * xime13 - vime13 * ximl13 - vime14 * ximl13 - viml14 * ximl13 + viml13 * xime14 - viml13 * ximl14) / (vime13 + viml13 + vime14 + viml14))
                uim14 = -((viml14 * xime13 - viml14 * ximl13 - viml14 * xime14 - vime13 * ximl14 - viml13 * ximl14 - vime14 * ximl14) / (vime13 + viml13 + vime14 + viml14))
                uq13 = -((- vql13 * xqe13 - vqe13 * xql13 - vqe14 * xql13 - vql14 * xql13 + vql13 * xqe14 - vql13 * xql14) / (vqe13 + vql13 + vqe14 + vql14))
                uq14 = -((vql14 * xqe13 - vql14 * xql13 - vql14 * xqe14 - vqe13 * xql14 - vql13 * xql14 - vqe14 * xql14) / (vqe13 + vql13 + vqe14 + vql14))
                # MLE for Consistent effect
                aim = -((- vime14 * xime13 - viml14 * xime13 + vime14 * ximl13 + viml14 * ximl13 - vime13 * xime14 - viml13 * xime14 + vime13 * ximl14 + viml13 * ximl14) / (vime13 + viml13 + vime14 + viml14))
                aq = -((- vqe14 * xqe13 - vql14 * xqe13 + vqe14 * xql13 + vql14 * xql13 - vqe13 * xqe14 - vql13 * xqe14 + vqe13 * xql14 + vql13 * xql14) / (vqe13 + vql13 + vqe14 + vql14))

                df_YCB_im = 1
                df_YCB_q = 1

                if all(k < p_min for k in [xime13, ximl13, xime14, ximl14]) or all(k > p_max for k in [xime13, ximl13, xime14, ximl14]) or (any(k == -99.0 for k in [xime13, xime14, ximl13, ximl14])):
                    df_YCB_im -= 1
                else:
                    ll_YCB_im += (-(xime13 - (uim13 + aim))**2) / (2 * vime13)
                    ll_YCB_im += (-(xime14 - (uim14 + aim))**2) / (2 * vime14)
                    ll_YCB_im += (-(ximl13 - (uim13))**2) / (2 * viml13)
                    ll_YCB_im += (-(ximl14 - (uim14))**2) / (2 * viml14)

                if all(k < p_min for k in [xqe13, xql13, xqe14, xql14]) or all(k > p_max for k in [xqe13, xql13, xqe14, xql14]) or (any(k == -99.0 for k in [xqe13, xqe14, xql13, xql14])):
                    df_YCB_q -= 1
                else:
                    ll_YCB_q += (-(xqe13 - (uq13 + aq))**2) / (2 * vqe13)
                    ll_YCB_q += (-(xqe14 - (uq14 + aq))**2) / (2 * vqe14)
                    ll_YCB_q += (-(xql13 - (uq13))**2) / (2 * vql13)
                    ll_YCB_q += (-(xql14 - (uq14))**2) / (2 * vql14)

                """Consistent effect of bulk across populations within a year."""
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
                        
                if all(k < p_min for k in [xime13, ximl13]) or all(k > p_max for k in [xime13, ximl13]) or all(k < p_min for k in [xqe13, xql13]) or all(k > p_max for k in [xqe13, xql13]) or any(k == -99.0 for k in [xime13, ximl13, xqe13, xql13]):
                    df_PCB_13 -= 1
                else:
                    ll_PCB_13 += (-(xime13 - (uim13 + a13))**2) / (2 * vime13)
                    ll_PCB_13 += (-(xqe13 - (uq13 +a13))**2) / (2 * vqe13)
                    ll_PCB_13 += (-(ximl13 - (uim13))**2) / (2 * viml13)
                    ll_PCB_13 += (-(xql13 - (uq13))**2) / (2 * vql13)
                
                if all(k < p_min for k in [xime14, ximl14]) or all(k > p_max for k in [xime14, ximl14]) or all(k < p_min for k in [xqe14, xql14]) or all(k > p_max for k in [xqe14, xql14]) or any(k == -99.0 for k in [xime14, ximl14, xqe14, xql14]):
                    df_PCB_14 -= 1
                else:
                    ll_PCB_14 += (-(xime14 - (uim14 + a14))**2) / (2 * vime14)
                    ll_PCB_14 += (-(xqe14 - (uq14 + a14))**2) / (2 * vqe14)
                    ll_PCB_14 += (-(ximl14 - (uim14))**2) / (2 * viml14)
                    ll_PCB_14 += (-(xql14 - (uq14))**2) / (2 * vql14)

    
                """Pop and bulk effect.  No effect of year."""
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
                
                if all(k < p_min for k in [xime13, ximl13]) or all(k > p_max for k in [xime13, ximl13]) or all(k < p_min for k in [xime14, ximl14]) or all(k > p_max for k in [xime14, ximl14]) or any(k == -99.0 for k in [xime13, ximl13, xime14, ximl14]):              
                    df_Y_im -= 2
                else:
                    for j, xx in enumerate(xime):
                        ll_B_im += (-(xx - uime)**2) / (2 * vime[j])
                    for j, xx in enumerate(ximl):
                        ll_B_im += (-(xx - uiml)**2) / (2 * viml[j])
                
                if all(k < p_min for k in [xqe13, xql13]) or all(k > p_max for k in [xqe13, xql13]) or all(k < p_min for k in [xqe14, xql14]) or all(k > p_max for k in [xqe14, xql14]) or any(k == -99.0 for k in [xqe13, xql13, xqe14, xql14]):
                    df_Y_q -= 2
                else:
                    for j, xx in enumerate(xqe):
                        ll_B_q += (-(xx - uqe)**2) / (2 * vqe[j])
                    for j, xx in enumerate(xql):
                        ll_B_q += (-(xx - uql)**2) / (2 * vql[j])
                
                """Most complex model.  All bulks have distinct mean which equals observed value.  Likelihoods reduce to 0"""
                ll_C = 0.0
                
                """Calculate Likelihood ratio test and p-values."""                
                if df_B_br13 > 0:
                    lrt_B_br13 = -2 * ll_Y_br13
                    p_B_br13 = chisqprob(lrt_B_br13, df_B_br13)
                    pvals_B_br13.append(p_B_br13)
                else:
                    p_B_br13 = "-"  # Set p-values and lrt's to missing data if site did not pass filter for this test
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
                    pime14 = math.sin(xime14 / 2)**2
                    piml14 = math.sin(ximl14 / 2)**2
                else:
                    p_B_im14 = "-"
                    lrt_B_im14 = "-"
                    pime14 = "-"
                    piml14 = "-"
                if df_B_im13 > 0:
                    lrt_B_im13 = -2 * ll_Y_im13
                    p_B_im13 = chisqprob(lrt_B_im13, df_B_im13)
                    pvals_B_im13.append(p_B_im13)
                    pime13 = math.sin(xime13 / 2)**2
                    piml13 = math.sin(ximl13 / 2)**2
                else:
                    p_B_im13 = "-"
                    lrt_B_im13 = "-"
                    pime13 = "-"
                    piml13 = "-"
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
                    pqe14 = math.sin(xqe14 / 2)**2
                    pql14 = math.sin(xql14 / 2)**2
                else:
                    p_B_q14 = "-"
                    lrt_B_q14 = "-"
                    pqe14 = "-"
                    pql14 = "-"
                if df_B_q13 > 0:
                    lrt_B_q13 = -2 * ll_Y_q13
                    p_B_q13 = chisqprob(lrt_B_q13, df_B_q13)
                    pvals_B_q13.append(p_B_q13)
                    pqe13 = math.sin(xqe13 / 2)**2
                    pql13 = math.sin(xql13 / 2)**2
                else:
                    p_B_q13 = "-"
                    lrt_B_q13 = "-"
                    pqe13 = "-"
                    pql13 = "-"
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
                else:
                    lrt_B_14_M12 = "-"  # Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)
                    lrt_B_14_M13 = "-"
                    lrt_B_14_M23 = "-"
                    p_B_14_M12 = "-"
                    p_B_14_M13 = "-"
                    p_B_14_M23 = "-"
                if df_P > 0:    
                    lrt_P = -2 * (ll_YB - ll_C)
                    p_P = chisqprob(lrt_P, df_P)
                    pvals_P.append(p_P)
                else:
                    lrt_P = "-"
                    p_P = "-"
                    
                # Determine pattern of polymorphism across samples for this site
                if p_B_im13 != "-" and p_B_q13 == "-" and p_B_br13 == "-":
                    polim += 1
                elif p_B_im13 == "-" and p_B_q13 != "-" and p_B_br13 == "-":
                    polq += 1
                elif p_B_im13 == "-" and p_B_q13 == "-" and p_B_br13 != "-":
                    polbr += 1
                elif p_B_im13 != "-" and p_B_q13 == "-" and p_B_br13 != "-":
                    polimbr += 1
                elif p_B_im13 == "-" and p_B_q13 != "-" and p_B_br13 != "-":
                    polqbr += 1
                elif p_B_im13 != "-" and p_B_q13 != "-" and p_B_br13 == "-":
                    polimq += 1
                elif p_B_im13 != "-" and p_B_q13 != "-" and p_B_br13 != "-":
                    polimqbr += 1
                    
                if p_B_im14 != "-" and p_B_q14 == "-":
                    polim1 += 1
                elif p_B_im14 == "-" and p_B_q14 != "-":
                    polq1 += 1
                elif p_B_im13 != "-" and p_B_q13 != "-":
                    polimq1 += 1
                    
                """Print to file"""
                if file1 != 0:
                    like.writerow([scaff, pos, lrt_B_br13, p_B_br13, "BR13"])
                    like.writerow([scaff, pos, lrt_B_im13, p_B_im13, "IM13"])
                    like.writerow([scaff, pos, lrt_B_q13, p_B_q13, "Q13"])
                    like.writerow([scaff, pos, lrt_B_im14, p_B_im14, "IM14"])
                    like.writerow([scaff, pos, lrt_B_q14, p_B_q14, "Q14"])

                if file2 != 0:
                    out2x.writerow([scaff, pos, lrt_B_im_M12, p_B_im_M12, lrt_B_im_M13, p_B_im_M13, lrt_B_im_M23, p_B_im_M23, "IM"])
                    out2x.writerow([scaff, pos, lrt_B_q_M12, p_B_q_M12, lrt_B_q_M13, p_B_q_M13, lrt_B_q_M23, p_B_q_M23, "Q"])
                    out2x.writerow([scaff, pos, lrt_B_13_M12, p_B_13_M12, lrt_B_13_M13, p_B_13_M13, lrt_B_13_M23, p_B_13_M23, "2013"])
                    out2x.writerow([scaff, pos, lrt_B_14_M12, p_B_14_M12, lrt_B_14_M13, p_B_14_M13, lrt_B_14_M23, p_B_14_M23, "2014"])
                    
                if file3 != 0:
                    out3x.writerow([lrt_B_im13, p_B_im13, lrt_B_q13, p_B_q13, lrt_B_13_M12, p_B_13_M12, lrt_B_13_M13, p_B_13_M13, lrt_B_13_M23, p_B_13_M23, pime13, piml13, pqe13, pql13])   
                    out4x.writerow([lrt_B_im14, p_B_im14, lrt_B_q14, p_B_q14, lrt_B_14_M12, p_B_14_M12, lrt_B_14_M13, p_B_14_M13, lrt_B_14_M23, p_B_14_M23, pime14, piml14, pqe14, pql14])   
                    out5x.writerow([scaff, pos, lrt_B_im14, p_B_im14, lrt_B_im13, p_B_im13, lrt_B_im_M12, p_B_im_M12, lrt_B_im_M13, p_B_im_M13, lrt_B_im_M23, p_B_im_M23, pime13, piml13, pime14, piml14])   
                    out6x.writerow([scaff, pos, lrt_B_q14, p_B_q14, lrt_B_q13, p_B_q13, lrt_B_q_M12, p_B_q_M12, lrt_B_q_M13, p_B_q_M13, lrt_B_q_M23, p_B_q_M23, pqe13, pql13, pqe14, pql14])   
               
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
    kk = 0    
    for m, p in enumerate(pvals_B_q14):
        if p < ((float(m) + 1.0) / float(len(pvals_B_q14))) * FDR:
            p_B_q14_cutoff = p
            kk = m + 1
                
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
    
    paramfile = open(OUTDIR + "EL_Likelihoods_" + timestr + ".params.txt", "w")
    paramfile.write("File ID: " + timestr + "\n")
    paramfile.write("Input: " + INPUT_FILE + "\n")
    paramfile.write("p_min: " + str(p_min1) + "\n")
    paramfile.write("p_max: " + str(p_max1) + "\n")
    paramfile.write("min_cov: " + str(min_cov) + "\n")
    paramfile.write("max_cov: " + str(max_cov) + "\n")
    paramfile.write("FDR: " + str(FDR) + "\n")
    paramfile.write("Number sites filtered: " + str(filtered) + "\n")
    paramfile.write("Total number of sites: " + str(sites) + "\n")
    paramfile.write("Program run time: " + str(runtime / 60) + " minutes \n")
    paramfile.write("p_P_cutoff = " + str(p_P_cutoff) + " ; num_sig_tests = " + str(o) + " ; num_tests = " +str(len(pvals_P)) + "\n")
    paramfile.write("p_Y_im_cutoff = " + str(p_Y_im_cutoff) + " ; num_sig_tests = " + str(x) + " ; num_tests = " +str(len(pvals_Y_im)) + "\n")
    paramfile.write("p_Y_q_cutoff = " + str(p_Y_q_cutoff) + " ; num_sig_tests = " + str(y) + " ; num_tests = " +str(len(pvals_Y_q)) + "\n")
    paramfile.write("p_B_q13_cutoff = " + str(p_B_q13_cutoff) + " ; num_sig_tests = " + str(qq) + " ; num_tests = " +str(len(pvals_B_q13)) + "\n")
    paramfile.write("p_B_im13_cutoff = " + str(p_B_im13_cutoff) + " ; num_sig_tests = " + str(t) + " ; num_tests = " +str(len(pvals_B_im13)) + "\n")
    paramfile.write("p_B_br13_cutoff = " + str(p_B_br13_cutoff) + " ; num_sig_tests = " + str(z) + " ; num_tests = " +str(len(pvals_B_br13)) + "\n")
    paramfile.write("p_B_q14_cutoff = " + str(p_B_q14_cutoff) + " ; num_sig_tests = " + str(kk) + " ; num_tests = " +str(len(pvals_B_q14)) + "\n")
    paramfile.write("p_B_im14_cutoff = " + str(p_B_im14_cutoff) + " ; num_sig_tests = " + str(tt) + " ; num_tests = " +str(len(pvals_B_im14)) + "\n")
    paramfile.write("p_B_im_M12_cutoff = " + str(p_B_im_M12_cutoff) + " ; num_sig_tests = " + str(h) + " ; num_tests = " +str(len(pvals_B_im_M12)) + "\n")
    paramfile.write("p_B_im_M13_cutoff = " + str(p_B_im_M13_cutoff) + " ; num_sig_tests = " + str(hh) + " ; num_tests = " +str(len(pvals_B_im_M13)) + "\n")
    paramfile.write("p_B_im_M23_cutoff = " + str(p_B_im_M23_cutoff) + " ; num_sig_tests = " + str(hhh) + " ; num_tests = " +str(len(pvals_B_im_M23)) + "\n")
    paramfile.write("p_B_q_M12_cutoff = " + str(p_B_q_M12_cutoff) + " ; num_sig_tests = " + str(j) + " ; num_tests = " +str(len(pvals_B_q_M12)) + "\n")
    paramfile.write("p_B_q_M13_cutoff = " + str(p_B_q_M13_cutoff) + " ; num_sig_tests = " + str(jj) + " ; num_tests = " +str(len(pvals_B_q_M13)) + "\n")
    paramfile.write("p_B_q_M23_cutoff = " + str(p_B_q_M23_cutoff) + " ; num_sig_tests = " + str(jjj) + " ; num_tests = " +str(len(pvals_B_q_M23)) + "\n")
    paramfile.write("p_B_13_M12_cutoff = " + str(p_B_13_M12_cutoff) + " ; num_sig_tests = " + str(k) + " ; num_tests = " +str(len(pvals_B_13_M12)) + "\n")
    paramfile.write("p_B_13_M13_cutoff = " + str(p_B_13_M13_cutoff) + " ; num_sig_tests = " + str(kk) + " ; num_tests = " +str(len(pvals_B_13_M13)) + "\n")
    paramfile.write("p_B_13_M23_cutoff = " + str(p_B_13_M23_cutoff) + " ; num_sig_tests = " + str(kkk) + " ; num_tests = " +str(len(pvals_B_13_M23)) + "\n")
    paramfile.write("p_B_14_M12_cutoff = " + str(p_B_14_M12_cutoff) + " ; num_sig_tests = " + str(l) + " ; num_tests = " +str(len(pvals_B_14_M12)) + "\n")
    paramfile.write("p_B_14_M13_cutoff = " + str(p_B_14_M13_cutoff) + " ; num_sig_tests = " + str(ll) + " ; num_tests = " +str(len(pvals_B_14_M13)) + "\n")
    paramfile.write("p_B_14_M23_cutoff = " + str(p_B_14_M23_cutoff) + " ; num_sig_tests = " + str(lll) + " ; num_tests = " +str(len(pvals_B_14_M23)) + "\n")
    paramfile.write(CustomMessage + "\n")    
    paramfile.close()
    
    print "FDR = ", FDR
    print "p_P_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_P_cutoff, x, len(pvals_P))     
    print "p_Y_im_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_Y_im_cutoff, x, len(pvals_Y_im))     
    print "p_Y_q_cutoff = %.19f ; num_sig_tests = %d ; num_tests = %d" % (p_Y_q_cutoff, y, len(pvals_Y_q))     
    print "p_B_im13_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_im13_cutoff, t, len(pvals_B_im13))
    print "p_B_q13_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_q13_cutoff, qq, len(pvals_B_q13))
    print "p_B_br13_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_br13_cutoff, z, len(pvals_B_br13))
    print "p_B_im14_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_im14_cutoff, tt, len(pvals_B_im14))
    print "p_B_q14_cutoff = %.19f ; num sig tests = %d ; num_tests = %d" % (p_B_q14_cutoff, kk, len(pvals_B_q14))
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
    print "number polymorphic 2013 IM Only = ", polim
    print "number polymorphic 2013 Q Only = ", polq
    print "number polymorphic 2013 BR Only = ", polbr
    print "number polymorphic 2013 IM Q = ", polimq
    print "number polymorphic 2013 IM BR = ", polimbr
    print "number polymorphic 2013 Q BR Only = ", polqbr
    print "number polymorphic 2013 IM Q BR = ", polimqbr
    print "number polymorphic 2013 total = ", polim + polq + polbr + polimq + polimbr + polqbr + polimqbr   
    print "number polymorphic 2014 IM Only = ", polim1 
    print "number polymorphic 2014 Q Only = ", polq1
    print "number polymorphic 2014 IM Q = ", polimq1
    print "number of sites fixed for alt = ", fixed
    print "filtered sites = ", filtered
    print "number of sites = ", sites
    print "program run time = ", runtime / 60
    print "num_lines = ", i
