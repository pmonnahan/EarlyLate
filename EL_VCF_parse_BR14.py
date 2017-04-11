# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 16:37:04 2016

@author: patrick
"""

#-------------------------------------------------------------------------------
# Name:        VCF parse
# Purpose:     Takes VCF and returns allele counts for each sample in the file outfilepath0.
#Notes: Run the program one time with permissive settings for Min_read_depth and Max_read_depth.  Look at the distribution of coverage per site in the file 'read_depth.txt'.  
#       Use the distribution of coverage, particularly the median, to set the min and max read depth.  
#-------------------------------------------------------------------------------

import sys

#AF = allele frequency analysis; RN = read number analysis

def main(inpath1="black0.txt", infilepath="/Volumes/avery/Early_Late_all.vcf", outfilepath0="/Users/patrick/Documents/Research/EarlyLate/BR14_counts.txt"):
    src  =open(infilepath, "rU")
    out0 =open(outfilepath0, "w")


#       Counters
    depth_cc=0
    raw_read_count=0
    line_num=0
    cutsites = 0
    noncutsites = 0
#    out1.write("Scaff"+"\t"+"pos"+"\t"+"DP"+"\t"+"\t"+"BRE_13"+"\t"+"BRL_13"+"\t"+"BR_count"+"\t"+"IME_14"+"\t"+"IML_14"+"\t"+"IME_13"+"\t"+"IML_13"+"\t"+"IM_count"+"\t"+"QE_14"+"\t"+"QL_14"+"\t"+"QE_13"+"\t"+"QL_13"+"\t"+"Q_count"+"\n")
    # evaluate contents of each line of input file
    for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t') 
        # print line_idx
        if len(cols) < 2:               ## This should be header
            continue
        elif cols[0] == "#CHROM":
            for i in range(len(cols)):
                print i,cols[i]
            out0.write("Scaff"+"\t"+"pos"+"\t"+"ref"+"\t"+"alt"+"\t"+"BRE14_R"+"\t"+"BRE14_A"+"\n")
        else:
            scaff=cols[0]                    
            position = int(cols[1])
            ref_base = cols[3]
            alt_base = cols[4]
            if len(alt_base) > 1:    # multiple bases at site
                pass #
            else:
                line_num+=1
                # AC=2;AF=0.500;AN=4;BaseQRankSum=-0.278;DP=11;Dels=0.00;...                                                    
                infoSNP = cols[7].split(";")
                ccREADS=infoSNP[3].split("=")
                if ccREADS[0]=="DP":
                    num_reads=int(ccREADS[1])
                else:
                    cxREADS=infoSNP[4].split("=")
                    if cxREADS[0]=="DP":
                        num_reads=int(cxREADS[1])
                    else:
                        print "wtf"
                depth_cc+=1
                raw_read_count+=num_reads   
                
                #Browder Ridge 2014
                if cols[9] != "./." : #Early Bulk
                    bases = cols[9].split(":")  
                    BRE14_ref,BRE14_alt = bases[1].split(",")
                else:
                    BRE14_ref=BRE14_alt=0
                
                out0.write(scaff+"\t"+str(position)+"\t"+ref_base+"\t"+alt_base+"\t"+str(BRE14_ref)+"\t"+str(BRE14_alt)+"\n")
            if line_idx % 100000 == 0:
                print line_idx
    
    out0.close()  
                     

if __name__ == "__main__":
    main()

