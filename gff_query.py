# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 19:44:29 2014

@author: patrick
"""

import csv

snp_dict={}
new_csv=open("/Users/patrick/Google Drive/Research/EarlyLate/EL_fdr05_queried.csv","wb")
writer=csv.writer(new_csv,dialect="excel")
with open("/Users/patrick/Google Drive/Research/EarlyLate/EL_fdr05.csv","rU") as sites_file:  
    scaffs=[]
    snps=[]
    counts=[]
    for i,site in enumerate(sites_file):
        if i>0:
            site=site.split(",")
            scaff=site[1].split("_")
            scaff=int(scaff[1])
            scaffs.append(scaff)
            snps.append(int(site[2]))
            counts.append(site[3:27])
#    for i in range(max(scaffs)):
    x=[[] for i in range(max(scaffs))]
    y=[[] for i in range(max(scaffs))]
    for i,scaff in enumerate(scaffs):
        x[scaff-1].append(snps[i])
        y[scaff-1].append(counts[i])
    
with open("/Users/patrick/Google Drive/Research/Mimulus/MimNewAnnotated.csv","rU") as gff:
    writer.writerow(["scaff","pos","lower","upper","bre_ref","bre_alt","brl_ref","brl_alt","ime_ref","ime_alt","iml_ref","iml_alt","qe_ref","qe_alt","ql_ref","ql_alt","Z_br","p_br","Z_im","p_im","Z_q","p_q","num_pops","Z_cum","p_cum","Fst_IMBR","Fst_IMQ","Fst_QBR","Gene id","description","GO.biological.process","GO.cellular.component","GO.molecular.function","Sequence Length","Hit.ACC","Evalue"])    
    for i,line in enumerate(gff):
        if i>0:
            line=line.split(",")
            scaff=int(line[1].split("_")[1])
            upper_bound=int(line[3])
            lower_bound=int(line[2])
            try:
                snps=x[scaff]
                counts=y[scaff]
            except IndexError:
                snps=[]
                counts=[]
            for i,snp in enumerate(snps):
                if snp > lower_bound and snp < upper_bound:
                    print counts
                    writer.writerow([scaff,snp,lower_bound,upper_bound,counts[i][0],counts[i][1],counts[i][2],counts[i][3],counts[i][4],counts[i][5],counts[i][6],counts[i][7],counts[i][8],counts[i][9],counts[i][10],counts[i][11],counts[i][12],counts[i][13],counts[i][14],counts[i][15],counts[i][16],counts[i][17],counts[i][18],counts[i][19],counts[i][20],counts[i][21],counts[i][22],counts[i][23],line[0],line[4],line[5],line[6],line[7],line[8],line[9],line[10]])
        
