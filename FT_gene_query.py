# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:52:28 2016

@author: patrick
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 16:41:54 2015

@author: patrick
"""

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
def main(input_file):
    import csv
    import scipy.stats as stats
    
    snp_dict={}
    
    #new_csv=open("/Users/patrick/Documents/Research/EarlyLate/EL_queried_all_snps_7-19-15.csv","wb")
    #writer=csv.writer(new_csv,dialect="excel")
    
    #p_Y_cutoff = 0.000148
    p_B_cutoff = 0.00015940722695
    #p_P_cutoff = 0.048506
    #p_Y_im_cutoff = 0.000037
    #p_YCB_im_cutoff = 0.000007
    #p_B_im_cutoff = 0.000028
    #p_Y_br_cutoff = 0.000010
    #p_B_br_cutoff = 0.000002
    #p_Y_q_cutoff = 0.000106
    #p_B_q_cutoff = 0.000180
    #p_YCB_q_cutoff = 0.000017
    
    sigFTG=0
    nsigFTG=0
    sigFTUS=0
    nsigFTUS=0
    sigNonFTG=0
    nsigNonFTG=0
    sigNonFTUS=0
    nsigNonFTUS=0
    sigFTBoth=0
    sigNonFTBoth=0
    
    INDIR="/Volumes/avery/Research/EarlyLate/GeneSearch/Data/"
    OUTDIR="/Volumes/avery/Research/EarlyLate/GeneSearch/Results/"
    
    INPUT_FILE="/Users/patrick/Documents/Research/Mimulus/FriedBlackman_GeneList.csv"
    Sig_QUERY_FILE=input_file
#    NS_QUERY_FILE=nonsig_input_file
    
    OUTPUT_FILE1=OUTDIR+Sig_QUERY_FILE[:-4]+"_FTqueried_sigB_inGenes.csv"
    OUTPUT_FILE2=OUTDIR+Sig_QUERY_FILE[:-4]+"_FTqueried_sigB_upstream.csv"
    OUTPUT_FILE3=OUTDIR+Sig_QUERY_FILE[:-4]+"_FTqueried_sigB_both.csv"
    #OUTPUT_FILE4=DIR+"FriedBlackman_GeneList_sigCounts.csv"
    #
    out1=open(OUTPUT_FILE1,"wb")
    out2=open(OUTPUT_FILE2,"wb")
    out3=open(OUTPUT_FILE3,"wb")
    #out4=open(OUTPUT_FILE4,"wb")
    #
    writer1=csv.writer(out1,dialect="excel",delimiter=",")
    writer2=csv.writer(out2,dialect="excel",delimiter=",")
    writer3=csv.writer(out3,dialect="excel",delimiter=",")
    #writer4=csv.writer(out4,dialect="excel",delimiter=",")
    
    with open("/Users/patrick/Documents/Research/Mimulus/Mguttatus_256_v2.0.gene_exons.txt","rU") as gff: 
        scaffs=[]
        upper=[]
        lower=[]
        sc_list=[]
        names=[]    
        numgenes=0
        for i,line in enumerate(gff):

            if i>1:
                line=line.strip("\n")
                line=line.split("\t")
                if line[2] == "gene":
                    info=line[8].split("=")
                    name=info[1]
                    name=name.split(".")
                    name=name[0]+"."+name[1]
                    upper_bound=float(line[4])
                    lower_bound=float(line[3])
                    scaff=line[0].split("_")[1]
                    numgenes+=1
                    try:
                        upper.append(upper_bound)
                        lower.append(lower_bound)
                        scaffs.append(int(scaff))
                        sc_list.append([[0,0],[0,0]])
                        names.append(name)
                    except ValueError:
                        pass
                    
        upper_list=[[] for i in range(max(scaffs)+1)]
        lower_list=[[] for i in range(max(scaffs)+1)]
        SC_list=[[] for i in range(max(scaffs)+1)]
        names_list=[[] for i in range(max(scaffs)+1)]
        
        for i,scaff in enumerate(scaffs):
            upper_list[scaff].append(upper[i])
            lower_list[scaff].append(lower[i])
            SC_list[scaff].append(sc_list[i])
            names_list[scaff].append(names[i])
    
    with open(INPUT_FILE,"rU") as FTG:  
        sites=[]
        ft_list=[]
        numsites=0
        for i,site in enumerate(FTG):
            if i==0:
                site=site.strip("\n")
                site=site.split(",")
                first_line=site
            elif i>0:
                site=site.strip("\n")
                site=site.split(",")
                ft_list.append(site[3])
                sites.append(site)
    
                  
    FT_list=[[] for i in range(max(scaffs)+1)]
    sites_list=[[] for i in range(max(scaffs)+1)]
    numFTGenes=0
    for ii,scaff in enumerate(names_list):
        ft_scaff=[]
        sites_scaff=[]
        for jj,name in enumerate(scaff):
            boo=True            
            for kk,ft_name in enumerate(ft_list):
                if ft_name == name and boo==True:
                    sites_scaff.append(sites[kk])
                    ft_scaff.append(1)
#                    print sites[kk],jj
#                    print lower_list[ii][jj]
                    numFTGenes+=1
                    boo=False                
            if boo==True:
                ft_scaff.append(0)
                sites_scaff.append(" ")
        FT_list[ii]=ft_scaff
        sites_list[ii]=sites_scaff

    with open(INDIR+Sig_QUERY_FILE,"rU") as sites_file:
        geneid=99.0
        usid=99.0
        for i,line in enumerate(sites_file):
#            if i%10000==0:
#                    print i
            if i==0:
                line=line.strip("\n")
                line=line.split(",")       
                writer1.writerow(line+first_line+["sigB"])
                writer2.writerow(line+first_line+["sigB"])
                writer3.writerow(line+first_line+["sigB"])
#                for k,j in enumerate(line):
#                    print k,j
        
            elif i>0:
                line1=line.strip("\n")
                line1=line1.split(",")
                line2=line.strip("\n")
                line2=line2.split(",")
                scaff=line1[0].split('"')
                scaff=scaff[1].split("_")
                scaff=int(scaff[1])
                snp=float(line1[1])
                
                try:
                    upper_bound=upper_list[scaff]
                    lower_bound=lower_list[scaff]
                    Site=sites_list[scaff]
                    ft=FT_list[scaff]
    
                except IndexError:
                    upper_bound=[]
                    lower_bound=[]
                    Site=[]
                    ft=[]

                for i,lower in enumerate(lower_bound):
                    if (snp > lower and snp <= upper_bound[i]):                      
                        #writer.writerow([scaff,snp,lower_bound,upper_bound,counts[i][0],counts[i][1],counts[i][2],counts[i][3],counts[i][4],counts[i][5],counts[i][6],counts[i][7],counts[i][8],counts[i][9],counts[i][10],counts[i][11],counts[i][12],counts[i][13],counts[i][14],counts[i][15],counts[i][16],counts[i][17],counts[i][18],counts[i][19],counts[i][20],counts[i][21],counts[i][22],counts[i][23],line[0],line[4],line[5],line[6],line[7],line[8],line[9],line[10]])                   
                        SC_list[scaff][i][0][0]+=1                        
                        if ft[i]==1 and lower!=geneid:                                                                                           
                            try:
                                line1+=Site[i]
                                line1+=Sig_QUERY_FILE
                                writer1.writerow(line1)
                                line1+=["inGene"]
                                writer3.writerow(line1)
                                geneid=lower
                            except IndexError:
                                pass                    
    
                        
                    elif (snp > lower-2000.0 and snp <= lower) or (snp >=upper_bound[i] and snp < upper_bound[i]+2000.0):
                        SC_list[scaff][i][1][0]+=1
                        if ft[i]==1 and lower!=usid:                            
                            try:
                                line2+=Site[i]
                                line2+=Sig_QUERY_FILE
                                writer2.writerow(line2)
                                line2+=["upstream"]
                                writer3.writerow(line2)
                                usid=lower
                            except IndexError:
                                pass  
                            
#    with open(INDIR+Sig_QUERY_FILE,"rU") as sites_file:
#        for i,line in enumerate(sites_file):
#            if i==0:
#                line=line.strip("\n")
#                line=line.split(",")       
#        
#            elif i>0:
#                line=line.strip("\n")
#                line=line.split(",")
#                scaff=line[0].split('"')
#                scaff=scaff[1].split("_")
#                scaff=int(scaff[1])
#                snp=float(line[1])
#                
#                try:
#                    upper_bound=upper_list[scaff]
#                    lower_bound=lower_list[scaff]
#                    site=sites_list[scaff]
#                    ft=FT_list[scaff]
#    
#                except IndexError:
#                    upper_bound=[]
#                    lower_bound=[]
#                    site=[]
#                    ft=[]
#                boo=True        
#                for i,lower in enumerate(lower_bound):
#                    if (snp > lower and snp <= upper_bound[i]) and boo==True:                                        
#                        SC_list[scaff][i][0][0]+=1                                           
#                            
#                    elif (snp > lower-2000 and snp <= lower) and boo==True:
#                        SC_list[scaff][i][1][0]+=1

                            
#    with open(INDIR+NS_QUERY_FILE,"rU") as sites_file:
#        for i,line in enumerate(sites_file):
##            if i%10000==0:
##                    print i
#            if i==0:
#                line=line.strip("\n")
#                line=line.split(",")       
##                for k,j in enumerate(line):
##                    print k,j
#        
#            elif i>0:
#                line=line.strip("\n")
#                line=line.split(",")
#                scaff=line[0].split('"')
#                scaff=scaff[1].split("_")
#                scaff=int(scaff[1])
#                snp=float(line[1])
#                
#                try:
#                    upper_bound=upper_list[scaff]
#                    lower_bound=lower_list[scaff]
#                    site=sites_list[scaff]
#                    ft=FT_list[scaff]
#    
#                except IndexError:
#                    upper_bound=[]
#                    lower_bound=[]
#                    site=[]
#                    ft=[]
#    #            print snp
#                boo=True        
#                for i,lower in enumerate(lower_bound):
#                    if (snp > lower and snp <= upper_bound[i]) and boo==True:                                        
#                        #writer.writerow([scaff,snp,lower_bound,upper_bound,counts[i][0],counts[i][1],counts[i][2],counts[i][3],counts[i][4],counts[i][5],counts[i][6],counts[i][7],counts[i][8],counts[i][9],counts[i][10],counts[i][11],counts[i][12],counts[i][13],counts[i][14],counts[i][15],counts[i][16],counts[i][17],counts[i][18],counts[i][19],counts[i][20],counts[i][21],counts[i][22],counts[i][23],line[0],line[4],line[5],line[6],line[7],line[8],line[9],line[10]])                   
#                        SC_list[scaff][i][0][1]+=1
#                        boo=False                                                                                                               
#                        
#                    elif (snp > lower-2000 and snp <= lower) and boo==True:
#                        SC_list[scaff][i][1][1]+=1
#                        boo=False

#    if other_file==1:
#        Other_File="/Users/patrick/Documents/Research/EarlyLate/EL_Likelihoods_Interactions_20160602-0925_nonsig"
#    else:
#        Other_File="/Users/patrick/Documents/Research/EarlyLate/EL_Likelihoods_IndTests_20160602-0925_nonsig.csv"
#
#    with open(Other_File,"rU") as sites_file:
#        for i,line in enumerate(sites_file):
##            if i%10000==0:
##                    print i
#            if i==0:
#                line=line.strip("\n")
#                line=line.split(",")       
##                for k,j in enumerate(line):
##                    print k,j
#        
#            elif i>0:
#                line=line.strip("\n")
#                line=line.split(",")
#                scaff=line[0].split("_")
#                scaff=int(scaff[1])
#                snp=float(line[1])
#                
#                try:
#                    upper_bound=upper_list[scaff]
#                    lower_bound=lower_list[scaff]
#                    site=sites_list[scaff]
#                    ft=FT_list[scaff]
#    
#                except IndexError:
#                    upper_bound=[]
#                    lower_bound=[]
#                    site=[]
#                    ft=[]
#    #            print snp
#                boo=True        
#                for i,lower in enumerate(lower_bound):
#                    if (snp > lower and snp <= upper_bound[i]) and boo==True:                                        
#                        #writer.writerow([scaff,snp,lower_bound,upper_bound,counts[i][0],counts[i][1],counts[i][2],counts[i][3],counts[i][4],counts[i][5],counts[i][6],counts[i][7],counts[i][8],counts[i][9],counts[i][10],counts[i][11],counts[i][12],counts[i][13],counts[i][14],counts[i][15],counts[i][16],counts[i][17],counts[i][18],counts[i][19],counts[i][20],counts[i][21],counts[i][22],counts[i][23],line[0],line[4],line[5],line[6],line[7],line[8],line[9],line[10]])                   
#                        SC_list[scaff][i][0][1]+=1
#                        boo=False                                                                                                               
#                        
#                    elif (snp > lower-2000 and snp <= lower) and boo==True:
#                        SC_list[scaff][i][1][1]+=1
#                        boo=False
        
                    
                        
    
    #linex=first_line+["IGsig","IGnotsig","USsig","USnotsig"]
    #writer4.writerow(linex)
    for j,scaff1 in enumerate(SC_list):
        for k,gene in enumerate(scaff1):
    #        gene+=SC_list[j][k]
    #        writer4.writerow(gene)
            if FT_list[j][k]==1:  
                if SC_list[j][k][0][0] > 0:
                    sigFTG+=1
                elif SC_list[j][k][0][1] > 0:
                    nsigFTG+=1
                if SC_list[j][k][1][0] > 0:
                    sigFTUS+=1
                elif SC_list[j][k][1][1] > 0:
                    nsigFTUS+=1
                if SC_list[j][k][0][0] > 0 or SC_list[j][k][1][0] > 0:
                    sigFTBoth+=1
            elif FT_list[j][k]==0:
                if SC_list[j][k][0][0] > 0:
                    sigNonFTG+=1
                elif SC_list[j][k][0][1] > 0:
                    nsigNonFTG+=1
                if SC_list[j][k][1][0] > 0:
                    sigNonFTUS+=1
                elif SC_list[j][k][1][1] > 0:
                    nsigNonFTUS+=1
                if SC_list[j][k][0][0] > 0 or SC_list[j][k][1][0] > 0:
                    sigNonFTBoth+=1
            
    Gene_oddsratio, Gene_pvalue = stats.fisher_exact([[sigFTG, numFTGenes-sigFTG], [sigNonFTG, numgenes-sigNonFTG-(numFTGenes-sigFTG)]])
    Upstream_oddsratio, Upstream_pvalue = stats.fisher_exact([[sigFTUS, numFTGenes-sigFTUS], [sigNonFTUS, numgenes-sigNonFTUS-(numFTGenes-sigFTUS)]])
    Both_oddsratio, Both_pvalue = stats.fisher_exact([[sigFTBoth, numFTGenes-sigFTBoth], [sigNonFTBoth, numgenes-sigNonFTBoth-(numFTGenes-sigFTBoth)]])
    
    
    print "Sig FT Genes = ",sigFTG
    print "NonSig FT Genes = ", numFTGenes-sigFTG
    print "Sig NonFTGenes = ", sigNonFTG
    print "NonSig NonFTGenes = ", numgenes-sigNonFTG-numFTGenes
    print "InGenes OR = ", Gene_oddsratio
    print "InGenes Pval = ", Gene_pvalue
    print "Sig FTUpstream = ", sigFTUS
    print "NonSig FTUpstream = ",numFTGenes-sigFTUS
    print "Sig NonFTUpstream = ", sigNonFTUS
    print "NonSig NonFTUpstream = ", numgenes-sigNonFTUS-numFTGenes
    print "Upstream OR = ", Upstream_oddsratio
    print "Upstream Pval = ", Upstream_pvalue
    print "Sig Either = ", sigFTBoth
    print "NonSig Either = ",numFTGenes-sigFTBoth
    print "Sig NonHits = ", sigNonFTBoth
    print "NonSig NonHits = ", numgenes-sigNonFTBoth-numFTGenes 
    print "Both OR = ", Both_oddsratio
    print "Both Pval = ", Both_pvalue
#    print "In Genes: Odds Ratio = %f P-value = %f" %(Gene_oddsratio, Gene_pvalue)
#    print "Upstream: Odds Ratio = %f P-value = %f" %(Upstream_oddsratio, Upstream_pvalue)
#    print "Both: Odds Ratio = %f P-value = %f" %(Both_oddsratio, Both_pvalue)
                        


                    
if __name__ == '__main__':
    import sys
    input_file=sys.argv[1].split("/")
    input_file=input_file[-1]
#    nonsig_input_file=sys.argv[2].split("/")
#    nonsig_input_file=nonsig_input_file[-1]
    main(input_file)