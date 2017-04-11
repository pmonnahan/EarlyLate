# matching Patrick for EL data 12/7/15
# Takes a file containing all contrasts
# 2a    take output from 2: IM.snps.txt and perform contrasts
# 2a    take output from 2a: IM.z_diffs.txt and calculate variance stats

import math
import sys
import csv


DIR="/Users/patrick/Documents/Research/EarlyLate/Current_Files/"
INPUT_FILE="Q.dz.minDP25.noInversions.txt"
OUTPUT_FILE="JackVars10-Q.minDP25.noInversions.csv"
out1=open(DIR+OUTPUT_FILE,"wb")
out1x=csv.writer(out1,delimiter=",",dialect='excel')
Num_Reps=10 #MUST BE A POWER OF TEN

out1x.writerow(["Replicate","num_lines","Con0","Con1","Con2","Con3","Con4","Con5"])
for rep in range(Num_Reps):
    # evaluate contents of each line of input file
    src  =open(DIR+INPUT_FILE, "rU")
    z0=[]
    z1=[]
    z2=[]
    z3=[]
    z4=[]
    z5=[]
    varm=[]
    zmeans=[0.0,0.0,0.0,0.0,0.0,0.0]
    zmeanSq=[0.0,0.0,0.0,0.0,0.0,0.0]
    vmeans=[0.0,0.0,0.0,0.0,0.0,0.0]
    lines=0    
    excluded=0
    for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t')
        try:
            if int(cols[1][-(str(Num_Reps).count(""))+2:]) != rep: #This throws out a portion of the data log10-proportional to Num_Reps.  Checks the last 1, 2 or 3 digist of the site position (correpsonding to 10, 100, or 1,000 reps).  Basically, throws out a tenth, hunredth, etc, for jacknife replicates.  
                if cols[2] != "-9":
                    z0.append(float(cols[2]))
                    zmeans[0]+=float(cols[2])
                    zmeanSq[0]+=float(cols[2])**2.0
                    vmeans[0]+=float(cols[3])
                if cols[4] != "-9":
                    z1.append(float(cols[4]))
                    zmeans[1]+=float(cols[4])
                    zmeanSq[1]+=float(cols[4])**2.0
                    vmeans[1]+=float(cols[5])
                if cols[6] != "-9":
                    z2.append(float(cols[6]))
                    zmeans[2]+=float(cols[6])
                    zmeanSq[2]+=float(cols[6])**2.0
                    vmeans[2]+=float(cols[7])    
                if cols[8] != "-9":
                    z3.append(float(cols[8]))    
                    zmeans[3]+=float(cols[8])
                    zmeanSq[3]+=float(cols[8])**2.0
                    vmeans[3]+=float(cols[9])
                if cols[10] != "-9":
                    z4.append(float(cols[10]))
                    zmeans[4]+=float(cols[10])
                    zmeanSq[4]+=float(cols[10])**2.0
                    vmeans[4]+=float(cols[11])
                if cols[12] != "-9":
                    z5.append(float(cols[12]))
                    zmeans[5]+=float(cols[12])
                    zmeanSq[5]+=float(cols[12])**2.0
                    vmeans[5]+=float(cols[13])
                lines+=1
            else:
                excluded+=1
        except ValueError:
            pass

    SnpNums=[len(z0),len(z1),len(z2),len(z3),len(z4),len(z5)]
    
    RDV = [vms/float(snps) for vms,snps in zip(vmeans,SnpNums)] # average read depth variance

    rz0=sorted(z0)
    rz1=sorted(z1)
    rz2=sorted(z2)
    rz3=sorted(z3)
    rz4=sorted(z4)
    rz5=sorted(z5)
    
    n250=rz0[int(SnpNums[0]/4)] 
    n500=rz0[int(SnpNums[0]/2)]
    n750=rz0[int(3*SnpNums[0]/4)]
    
    n251=rz1[int(SnpNums[1]/4)]
    n501=rz1[int(SnpNums[1]/2)]
    n751=rz1[int(3*SnpNums[1]/4)]
    
    n252=rz2[int(SnpNums[2]/4)]
    n502=rz2[int(SnpNums[2]/2)]
    n752=rz2[int(3*SnpNums[2]/4)]
    
    n253=rz3[int(SnpNums[3]/4)]
    n503=rz3[int(SnpNums[3]/2)]
    n753=rz3[int(3*SnpNums[3]/4)]
    
    n254=rz4[int(SnpNums[4]/4)]
    n504=rz4[int(SnpNums[4]/2)]
    n754=rz4[int(3*SnpNums[4]/4)]
    
    n255=rz5[int(SnpNums[5]/4)]
    n505=rz5[int(SnpNums[5]/2)]
    n755=rz5[int(3*SnpNums[5]/4)]
    
    v0=((n750-n250)/1.349)**2 - RDV[0]
    v1=((n751-n251)/1.349)**2 - RDV[1]
    v2=((n752-n252)/1.349)**2 - RDV[2]
    v3=((n753-n253)/1.349)**2 - RDV[3]
    v4=((n754-n254)/1.349)**2 - RDV[4]
    v5=((n755-n255)/1.349)**2 - RDV[5]
    
    src.close()
    out1x.writerow([rep,lines,v0,v1,v2,v3,v4,v5])
    print rep, lines, excluded

#print "Contrast ",comparison ,"snps = ",SnpCC
#print "Mean, Var z = ",zmeans[0]/float(SnpCC),zmeans[1]/float(SnpCC)-(zmeans[0]/float(SnpCC))**2.0
#print "Read depth Var = ", RDV

#print "Z percentiles ",n25,n50,n75
#print "total variance in Z (based on IQR) ",((n75-n25)/1.349)**2
#print "estimated core var ",((n75-n25)/1.349)**2 - RDV

#ranked_varm=sorted(varm)

#out1.write(outstring+'\n')
