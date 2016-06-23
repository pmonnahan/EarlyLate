# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 19:24:52 2015

@author: patrick
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 19:44:29 2014

@author: patrick
"""
##NEEDS TO BE FIXED.  EVERY HIT WITHIN A GENE IS BEING REPORTED...RESULTING SOMETIMES IN HUNDREDS OF HITS PER GENE.  THIS PRODUCES A MASSIVE FILE...>25 GIGS.  ONLY NEED TO REPORT ONE HIT (PERHAPS MOST SIGNIFICANT ONE) PER GENE.

def main(input_file,sigcol):
    import csv
    
    snp_dict={}
    
    #new_csv=open("/Users/patrick/Documents/Research/EarlyLate/EL_queried_all_snps_7-19-15.csv","wb")
    #writer=csv.writer(new_csv,dialect="excel")
    
    OUTDIR="/Volumes/avery/Research/EarlyLate/GeneSearch/Results/"
    INDIR="/Volumes/avery/Research/EarlyLate/GeneSearch/"
    
    INPUT_FILE=input_file
    
    #OUTPUT_FILE1=DIR+INPUT_FILE[:-4]+"_gffqueriedB.csv"
    
    #out1=open(OUTPUT_FILE1,"wb")
    
    #writer1=csv.writer(out1,dialect="excel",delimiter=",")
    
    numsigB=0
    
    numhits=0
    
    with open(INDIR+INPUT_FILE,"rU") as sites_file:  
        scaffs=[]
        snps=[]
        SigB=[]
        for i,site in enumerate(sites_file):
            if i == 0:  
                site=site.strip("\n")
                site=site.strip("\r")
                site=site.split(",")
                for k,j in enumerate(site):
                    print k,j
            if i>1:
                site=site.strip('\n')
                site=site.split(",")
                scaff=site[0].split("_")
                scaff=int(scaff[1])
                if scaff==999:
                    pass
                else:
                    try:
                        sigB=site[sigcol]
                    except IndexError:
                        pass
    
                snps.append(float(site[1]))
                scaffs.append(scaff)
                SigBim.append(sigBim)
                SigBq.append(sigBq)
                SigBbr.append(sigBbr)
    
        total = i+1
    
        try:
            x=[[] for i in range(max(scaffs)+1)]
        except ValueError:
            print site
    
        z=[[] for i in range(max(scaffs)+1)]
        zz=[[] for i in range(max(scaffs)+1)]
        zzz=[[] for i in range(max(scaffs)+1)]
        for i,scaff in enumerate(scaffs):
            x[scaff].append(snps[i])
            z[scaff].append(SigBim[i])
            zz[scaff].append(SigBq[i])
            zzz[scaff].append(SigBbr[i])
     
    with open("/Users/patrick/Documents/Research/Mimulus/Mguttatus_256_v2.0.gene_exons.txt","rU") as gff:
    #    writer1.writerow(["scaff","pos","lower","upper","Gene id","description","GO.biological.process","GO.cellular.component","GO.molecular.function","Sequence Length","Hit.ACC","Evalue","sigBim","sigBq","sigBbr"])        
        for i,line in enumerate(gff):
            if i%100==0:
                print i
            if i>1:
                line=line.split("\t")
                scaff=int(line[0].split("_")[1])
                upper_bound=float(line[4])
                lower_bound=float(line[3])
                try:
                    snps=x[scaff]
                    sigBim=z[scaff]
                    sigBq=zz[scaff]
                    sigBbr=zzz[scaff]
                except IndexError:
                    snps=[]
                    counts=[]
                boo=True
                for i,snp in enumerate(snps):
                    if boo==True:
                        if (snp > lower_bound and snp < upper_bound and line[2]=="gene") and (sigBim[i]=="1" or sigBq[i]=="1" or sigBbr[i]=="1"):
                            #writer.writerow([scaff,snp,lower_bound,upper_bound,counts[i][0],counts[i][1],counts[i][2],counts[i][3],counts[i][4],counts[i][5],counts[i][6],counts[i][7],counts[i][8],counts[i][9],counts[i][10],counts[i][11],counts[i][12],counts[i][13],counts[i][14],counts[i][15],counts[i][16],counts[i][17],counts[i][18],counts[i][19],counts[i][20],counts[i][21],counts[i][22],counts[i][23],line[0],line[4],line[5],line[6],line[7],line[8],line[9],line[10]])
                            boo=False                    
    #                    writer1.writerow([scaff,snp,lower_bound,upper_bound]+line+[sigBim[i],sigBq[i],sigBbr[i]])
                            numhits+=1
                        
    
                
                    
    print "Number of hits to query = ", numhits
    print "Number Total = ", total                   


if __name__ == '__main__':
    import sys
    input_file=sys.argv[1] 
    sigcol=sys.argv[2] #Column of file that contains 0/1 indicator of significance.
    main(input_file)