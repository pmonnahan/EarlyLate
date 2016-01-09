# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 17:37:06 2015

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

a=open("/Users/patrick/Documents/Research/Mimulus/FTGenes_AthalianaHomologs_PJM_intervals1.csv","w")
writer=csv.writer(a,dialect="excel",delimiter=",")

with open("/Users/patrick/Documents/Research/Mimulus/Friedman_Mguttatus_140_Flowering_genes_location.txt","rU") as sites_file:  
    scaffs=[]
    upper=[]
    mid=[]
    lower=[]
    sites=[]
    for i,site in enumerate(sites_file):
        if i==0:
            site=site.strip("\n")
            site=site.split("\t")
            first_line=site
        elif i>0:
            site=site.strip("\n")
            site=site.split("\t")
            scaff=site[7]
            middle=int(site[8])
            length=abs(int(site[5])-int(site[6]))
            upper_bound=middle+length
            lower_bound=middle-length
            try:
                upper.append(upper_bound)
                lower.append(lower_bound)
                mid.append(middle)
                scaffs.append(int(scaff))
                sites.append(site)
            except ValueError:
                pass
    
#    for i in range(max(scaffs)):
    upper_list=[[] for i in range(max(scaffs)+1)]
    lower_list=[[] for i in range(max(scaffs)+1)]
    mid_list=[[] for i in range(max(scaffs)+1)]
    site_list=[[] for i in range(max(scaffs)+1)]
    for i,scaff in enumerate(scaffs):
        upper_list[scaff].append(upper[i])
        lower_list[scaff].append(lower[i])
        mid_list[scaff].append(mid[i])
        site_list[scaff].append(sites[i])
print upper_list[2][1]
print lower_list[2][1]
print site_list[2][1]       

with open("/Users/patrick/Documents/Research/Mimulus/FTGenes_AthalianaHomologs_PJM.csv","rU") as gff:
    for i,line in enumerate(gff):
        if i==0:
            line=line.strip("\n")
            line=line.split(",")       
            writer.writerow(line+first_line+["v2Start","v2Stop"])    
    
        elif i>0:
            line=line.strip("\n")
            line=line.split(",")
            scaff=line[2].split("_")
            scaff=int(scaff[1])
            snp=int(line[9])
            
            
            try:
                upper_bound=upper_list[scaff]
                middle=mid_list[scaff]
                lower_bound=lower_list[scaff]
                site=site_list[scaff]
            except IndexError:
                upper_bound=[]
                middle=[]
                lower_bound=[]
                site=[]
#            print snp
            boo=True
            
            for i,midd in enumerate(middle):
#                if snp == 22825054:
#                    print upper_bound[i],snp,lower_bound[i],midd
#                    print scaff
                if (snp > lower_bound[i] and snp <= midd):                    
                    #writer.writerow([scaff,snp,lower_bound,upper_bound,counts[i][0],counts[i][1],counts[i][2],counts[i][3],counts[i][4],counts[i][5],counts[i][6],counts[i][7],counts[i][8],counts[i][9],counts[i][10],counts[i][11],counts[i][12],counts[i][13],counts[i][14],counts[i][15],counts[i][16],counts[i][17],counts[i][18],counts[i][19],counts[i][20],counts[i][21],counts[i][22],counts[i][23],line[0],line[4],line[5],line[6],line[7],line[8],line[9],line[10]])
                    try:                    
                        line+=site[i]+[lower_bound[i],midd]
                        writer.writerow(line)
                    except IndexError:
                        pass
                    boo=False
                elif (snp > midd and snp <= upper_bound[i]):
                    try:
                        line+=site[i]+[midd,upper_bound[i]] 
                        writer.writerow(line)
                    except IndexError:
                        pass
                    boo=False
                    
            if boo==True:
                boo2=True
                for i,ss in enumerate(site_list):
                    if len(ss)!=0:
                        for k,jj in enumerate(ss):
                            if jj[1][0:-1]==line[0][0:-1] and boo2==True:
                                line1=line+jj
                                boo2=False
                                writer.writerow(line1)

                if boo2==True:
                    writer.writerow(line)

                    
                    