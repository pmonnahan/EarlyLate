# Calculates the core contrast variance (i.e. total var in contrast minus read depth variance).  
#This becomes the z vector for the general least squares calculation of library specific variance.

import math
import sys


DIR="/Users/patrick/Documents/Research/EarlyLate/"
#comparison = 5 # 0 to 5
src  =open(DIR+"Q.dz.txt", "rU")
#out1 =open(comparison+".IM.txt","w")


# evaluate contents of each line of input file
z=[]
varm=[]
zmeans=[0.0,0.0]
vmeans=[0.0,0.0]

for comparison in range (0,6):
    src  =open(DIR+"IM.dz.txt", "rU")
    for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t') 
    		
        if cols[2*int(comparison)+2] != "-9":
    
            z.append(float(cols[2*int(comparison)+2]))
            varm.append(float(cols[2*int(comparison)+3]))	
            zmeans[0]+=float(cols[2*int(comparison)+2])
            zmeans[1]+=float(cols[2*int(comparison)+2])**2.0
            vmeans[0]+=float(cols[2*int(comparison)+3])
            vmeans[1]+=float(cols[2*int(comparison)+3])**2.0
    
    SnpCC=len(z)
    
    print "Contrast ",comparison ,"snps = ",SnpCC
    print "Mean, Var z = ",zmeans[0]/float(SnpCC),zmeans[1]/float(SnpCC)-(zmeans[0]/float(SnpCC))**2.0
    RDV= vmeans[0]/float(SnpCC) # average read depth variance
    print "Read depth Var = ", RDV
    
    ranked_z=sorted(z)
    n25=ranked_z[int(SnpCC/4)]
    n50=ranked_z[int(SnpCC/2)]
    n75=ranked_z[int(3*SnpCC/4)]
    print "Z percentiles ",n25,n50,n75
    print "total variance in Z (based on IQR) ",((n75-n25)/1.349)**2
    print "estimated core var ",((n75-n25)/1.349)**2 - RDV

#ranked_varm=sorted(varm)

#out1.write(outstring+'\n')
