# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 19:16:16 2015

@author: patrick
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 19:44:29 2014

@author: patrick
"""

import csv

snp_dict={}

#new_csv=open("/Users/patrick/Documents/Research/EarlyLate/EL_queried_all_snps_7-19-15.csv","wb")
#writer=csv.writer(new_csv,dialect="excel")

DIR="/Users/patrick/Documents/Research/Mimulus/"

INPUT_FILE="Friedman_Mguttatus_140_Flowering_genes_location.csv"
synonyms_file="Mguttatus_256_v2.0.synonym.txt"

OUTPUT_FILE1=DIR+INPUT_FILE[:-4]+"_queried_simple.csv"

out1=open(OUTPUT_FILE1,"wb")

writer1=csv.writer(out1,dialect="excel",delimiter=",")

#a=open(DIR+INPUT_FILE[:-4]+"_queried_GeneIDs_Y.txt","wb")
#b=open(DIR+INPUT_FILE[:-4]+"_queried_GeneIDs_B.txt","wb")
#c=open(DIR+INPUT_FILE[:-4]+"_queried_GeneIDs_P.txt","wb")

#numsigY=0
#numsigB=0
#numsigP=0

numhits=0

with open(DIR+INPUT_FILE,"rU") as sites_file:  
    names=[]
    sites=[]
    for i,site in enumerate(sites_file):
        site=site.strip('\n')
        site=site.split(",")
        if i==0:
            head=site
        else:
            names.append(site[0])
            sites.append(site)
            
with open(DIR+synonyms_file,"rU") as syns:
    newnames=[[] for _ in range(len(names))]
    for k,syn in enumerate(syns):
        syn=syn.strip("\n")
        syn=syn.split("\t")
#        index=names.index(syn[0])
        try:
            index=names.index(str(syn[1][:-2]))
            newnames[index].append(syn[0][:-2])
        except ValueError:
            pass
        
ii=0
for jj,nich in enumerate(newnames):
    for kk in nich:
        if kk ==None:
            print names[jj]
            ii+=1
        
 
with open(DIR+"Mguttatus_256_v2.0.gene_exons.txt","rU") as gff: 
    for i,line in enumerate(gff):
        if i == 0:
            writer1.writerow(["gffscaff","phytozome?","gene","start","stop","?","??","???","INFO","v2Name"]+head)
        if i>1:
            line=line.strip("\n")
            line=line.split("\t")
            info=line[8].split("=")
            name=info[1]
            name=name.split(".")
            name=name[0]+"."+name[1]
            for j,kkk in enumerate(newnames):
                for l,geneid in enumerate(kkk):
                    if geneid==name and line[2]=="gene":
                        writer1.writerow(line+[geneid]+sites[j])
                        numhits+=1
                        newnames[j].pop(l)
            if i%10000==0:
                print i
                
print "Number of hits to gff = ", numhits
print "Number of ids searched for = ", len(names)                   