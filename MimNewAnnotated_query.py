# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 12:12:58 2015

@author: patrick
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 19:44:29 2014

@author: patrick
"""
def main(input_file):    
    import csv
    
    snp_dict={}
    
    #new_csv=open("/Users/patrick/Documents/Research/EarlyLate/EL_queried_all_snps_7-19-15.csv","wb")
    #writer=csv.writer(new_csv,dialect="excel")
    
    #DIR="/Users/patrick/Documents/Research/EarlyLate/"
    OUTDIR="/Volumes/avery/Research/EarlyLate/GeneSearch/Results/"
    INDIR="/Volumes/avery/Research/EarlyLate/GeneSearch/Data/"
    INPUT_FILE=input_file
    
    OUTPUT_FILE1=OUTDIR+INPUT_FILE[:-4]+"_queried.csv"
    
    out1=open(OUTPUT_FILE1,"wb")
    
    writer1=csv.writer(out1,dialect="excel",delimiter=",")
    
    #a=open(DIR+INPUT_FILE[:-4]+"_queried_GeneIDs_Y.txt","wb")
    #b=open(DIR+INPUT_FILE[:-4]+"_queried_GeneIDs_B.txt","wb")
    #c=open(DIR+INPUT_FILE[:-4]+"_queried_GeneIDs_P.txt","wb")
    d=open(OUTDIR+INPUT_FILE[:-4]+"_queried_GeneIDs.txt","wb")
    
    numsigB=0
    
    numhits=0
    
    with open("/Users/patrick/Documents/Research/Mimulus/MimNewAnnotated.csv","rU") as gff:
        writer1.writerow(["scaff","pos","lower","upper","pval","Gene id","description","GO.biological.process","GO.cellular.component","GO.molecular.function","Sequence Length","Hit.ACC","Evalue"])    
        scaffs=[]
        upper=[]
        lower=[]
        geneids=[]
        geneinfo=[]
        for i,line in enumerate(gff):
#            if i == 0 :
#                print line
            if i>0:
                line=line.split(",")
                scaff=int(line[1].split("_")[1])
                scaffs.append(scaff)
                upper.append(float(line[3]))
                lower.append(float(line[2]))
                geneids.append(line[0])
                geneinfo.append(line)
        
        upper_list=[[] for i in range(max(scaffs)+1)]        
        lower_list=[[] for i in range(max(scaffs)+1)]
        geneid_list=[[] for i in range(max(scaffs)+1)]
        gene_info=[[] for i in range(max(scaffs)+1)]
        for i,scaff in enumerate(scaffs):
            upper_list[scaff].append(upper[i])
            lower_list[scaff].append(lower[i])
            geneid_list[scaff].append(geneids[i])
            gene_info[scaff].append(geneinfo[i])
            
    
    numhits=0
    with open(INDIR+INPUT_FILE,"rU") as sites_file:  
        previd=None
        for i,site in enumerate(sites_file):
#            if i%100000==0:
#                print i
            if i == 0:  
                site=site.strip("\n")
                site=site.strip("\r")
                site=site.split(",")
#                for k,j in enumerate(site):
#                    print k,j
            if i>1:
                site=site.strip('\n')
                site=site.split(",")
#                scaff=site[0].split('"')
#                scaff=int(scaff[1].split("_")[1])
                scaff=int(site[0])
                
                snp=float(site[1])
    #            if scaff==999:
    #                pass
    #            else:
    #                sigP=site[40]
    #                sigB=site[39]
    #                sigY=site[38]
    
                try:
                    upper=upper_list[scaff]
                    lower=lower_list[scaff]
                    geneids=geneid_list[scaff]
                    info=gene_info[scaff]
                except IndexError:
                    upper=[]
                    lower=[]
                    geneids=[]
                    info=[]
                for j,lower_bound in enumerate(lower):
                    if snp > lower_bound-2000 and snp <= upper[j]+2000:
#                        print snp,lower_bound,upper[j]
                        #writer.writerow([scaff,snp,lower_bound,upper_bound,counts[i][0],counts[i][1],counts[i][2],counts[i][3],counts[i][4],counts[i][5],counts[i][6],counts[i][7],counts[i][8],counts[i][9],counts[i][10],counts[i][11],counts[i][12],counts[i][13],counts[i][14],counts[i][15],counts[i][16],counts[i][17],counts[i][18],counts[i][19],counts[i][20],counts[i][21],counts[i][22],counts[i][23],line[0],line[4],line[5],line[6],line[7],line[8],line[9],line[10]])
    #                    print line                    
                        writer1.writerow([scaff,snp,lower_bound,upper[j]]+info)
                        if geneids[j]!=previd:
                            d.writelines(str(geneids[j])+"\n")
                            numhits+=1
                        previd=geneids[j]
    #                    if sigY[i] == '0':
    #                        pass
    #                    elif sigY[i] == '1':
    #                        a.writelines(str(geneid)+"\n")
    #                        numsigY+=1
    #                    else:
    #                        print "WTF1", sigY
    #                    if sigB[i] == '0':
    #                        pass
    #                    elif sigB[i] == '1':
    #                        b.writelines(str(geneid)+"\n")
    #                        numsigB+=1
    #                    else:
    #                        print "WTF2", sigB
    #                    if sigP[i] == '0':
    #                        pass
    #                    elif sigP[i] == '1':
    #                        c.writelines(str(geneid)+"\n") 
    #                        numsigP+=1
    #                    else:
    #                        print "WTF3", sigP
                
                    
    print "Number of hits to query = ", numhits
    print "Number of B hits = ", numsigB
    #print "Number Total = ", total                   
    #a.close()
    #b.close()
    #c.close()
    d.close()

if __name__ == '__main__':
    import sys
    input_file=sys.argv[1].split("/")
    input_file=input_file[-1]    
    main(input_file)