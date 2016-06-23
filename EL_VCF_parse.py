#-------------------------------------------------------------------------------
# Name:        VCF parse
# Purpose:     Takes VCF and returns allele counts for each sample in the file outfilepath0.
#Notes: Run the program one time with permissive settings for Min_read_depth and Max_read_depth.  Look at the distribution of coverage per site in the file 'read_depth.txt'.  
#       Use the distribution of coverage, particularly the median, to set the min and max read depth.  
#-------------------------------------------------------------------------------

import sys

#AF = allele frequency analysis; RN = read number analysis

def main(inpath1="black0.txt", infilepath="/Users/patrick/Documents/Research/EarlyLate/EL.massive.Jan16.vcf", outfilepath0="/Users/patrick/Documents/Research/EarlyLate/counts_massiveJan16_AF250_22-6-16.txt", outfilepath1="/Users/patrick/Documents/Research/EarlyLate/read_depth_massiveJan16_22-6-16.txt", outfilepath6="/Users/patrick/Documents/Research/EarlyLate/counts_parsed_8-11-15.txt",outfilepath10="/Users/patrick/Documents/Research/EarlyLate/counts_dt_SandP_RN250_22-6-16.txt"):
    src  =open(infilepath, "rU")
#    in0  =open(inpath1, "rU")
    out0 =open(outfilepath0, "w")
    out1 =open(outfilepath1, "w")
    out10=open(outfilepath10,"w")
         
    

#       USER DEFINED
    
#   These are the 99th percentile for read depth for each bulk
    BRE14_cutoff = 69
    BRE13_cutoff = 91
    BRL13_cutoff = 68
    IME13_cutoff = 63
    IML13_cutoff = 104
    IME14_cutoff = 78
    IML14_cutoff = 77
    QE13_cutoff = 168
    QL13_cutoff = 96
    QE14_cutoff = 90
    QL14_cutoff = 89


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
            out1.write("scaff"+"\t"+"position"+"\t"+"num_reads"+"\t"+"BRE13"+"\t"+"BRL13"+"\t"+"BR13"+"\t"+"IME14"+"\t"+"IML14"+"\t"+"IME13"+"\t"+"IML13"+"\t"+"IM"+"\t"+"QE14"+"\t"+"QL14"+"\t"+"QE13"+"\t"+"QL13"+"\t"+"Q"+"\n")
            out0.write("Scaff"+"\t"+"pos"+"\t"+"ref"+"\t"+"alt"+"\t"+"BRE13_R"+"\t"+"BRE13_A"+"\t"+"BRL13_R"+"\t"+"BRL13_A"+"\t"+"IME14_R"+"\t"+"IME14_A"+"\t"+"IML14_R"+"\t"+"IML14_A"+"\t"+"IME13_R"+"\t"+"IME13_A"+"\t"+"IML13_R"+"\t"+"IML13_A"+"\t"+"QE14_R"+"\t"+"QE14_A"+"\t"+"QL14_R"+"\t"+"QL14_A"+"\t"+"QE13_R"+"\t"+"QE13_A"+"\t"+"QL13_R"+"\t"+"QL13_A"+"\n")
            out10.write("Scaff"+"\t"+"pos"+"\t"+"ref"+"\t"+"alt"+"\t"+"BRE13_R"+"\t"+"BRE13_A"+"\t"+"BRL13_R"+"\t"+"BRL13_A"+"\t"+"IME14_R"+"\t"+"IME14_A"+"\t"+"IML14_R"+"\t"+"IML14_A"+"\t"+"IME13_R"+"\t"+"IME13_A"+"\t"+"IML13_R"+"\t"+"IML13_A"+"\t"+"QE14_R"+"\t"+"QE14_A"+"\t"+"QL14_R"+"\t"+"QL14_A"+"\t"+"QE13_R"+"\t"+"QE13_A"+"\t"+"QL13_R"+"\t"+"QL13_A"+"\n")
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
                
                #Browder Ridge 2013
                if cols[9] != "./." : #Early Bulk
                    bases = cols[9].split(":")  
                    BRE13_ref,BRE13_alt = bases[1].split(",")
                else:
                    BRE13_ref=BRE13_alt=0

                if cols[10] != "./.": #Late Bulk
                    bases = cols[10].split(":") 
                    BRL13_ref,BRL13_alt = bases[1].split(",")
                else:
                    BRL13_ref=BRL13_alt=0
                
                #Iron Mountain 2013
                if cols[11] != "./.":#Early Bulk
                    bases = cols[11].split(":")  
                    IME13_ref,IME13_alt = bases[1].split(",")  
                else:
                    IME13_ref=IME13_alt=0
                    
                if cols[12] !="./.": #Late Bulk
                    bases = cols[12].split(":") 
                    IML13_ref,IML13_alt = bases[1].split(",")
                else:
                    IML13_ref=IML13_alt=0
                
                #Iron Mountain 2014 
                if cols[13] != "./.": #Early Bulk
                    bases = cols[13].split(":")  
                    IME14_ref,IME14_alt = bases[1].split(",")
                else:
                    IME14_ref=IME14_alt=0

                if cols[14] != "./.":
                    bases = cols[14].split(":") #Late Bulk
                    IML14_ref,IML14_alt = bases[1].split(",")
                else:
                    IML14_ref=IML14_alt=0
                    
                #Quarry 2013
                if cols[18] != "./.": #Early Bulk
                    bases = cols[18].split(":")  
                    QE13_ref,QE13_alt = bases[1].split(",")
                else:
                    QE13_ref=QE13_alt=0

                if cols[15] != "./.": #Late Bulk
                    bases = cols[15].split(":") 
                    QL13_ref,QL13_alt = bases[1].split(",")
                else:
                    QL13_ref=QL13_alt=0
                    
                #Quarry 2014
                if cols[16] != "./.":#Early Bulk
                    bases = cols[16].split(":") 
                    QE14_ref,QE14_alt = bases[1].split(",")
                else:
                    QE14_ref=QE14_alt=0   
                     
                if cols[17] != "./.":
                    bases = cols[17].split(":") #Late
                    QL14_ref,QL14_alt = bases[1].split(",")
                else:
                    QL14_ref=QL14_alt=0
                
                BRE13_count=int(BRE13_ref)+int(BRE13_alt)
                BRL13_count=int(BRL13_ref)+int(BRL13_alt)
                IME14_count=int(IME14_ref)+int(IME14_alt)
                IME13_count=int(IME13_ref)+int(IME13_alt)
                IML13_count=int(IML13_ref)+int(IML13_alt)
                IML14_count=int(IML14_ref)+int(IML14_alt)
                QE14_count=int(QE14_ref)+int(QE14_alt)
                QE13_count=int(QE13_ref)+int(QE13_alt)
                QL13_count=int(QL13_ref)+int(QL13_alt)
                QL14_count=int(QL14_ref)+int(QL14_alt)
                
                if (any(k > 250 for k in [BRE13_count, BRL13_count, IME13_count, IML13_count, IME14_count, IML14_count, QE13_count, QL13_count, QE14_count, QL14_count])):
                    out10.write(scaff+"\t"+str(position)+"\t"+ref_base+"\t"+alt_base+"\t"+str(BRE13_ref)+"\t"+str(BRE13_alt)+"\t"+str(BRL13_ref)+"\t"+str(BRL13_alt)+"\t"+str(IME14_ref)+"\t"+str(IME14_alt)+"\t"+str(IML14_ref)+"\t"+str(IML14_alt)+"\t"+str(IME13_ref)+"\t"+str(IME13_alt)+"\t"+str(IML13_ref)+"\t"+str(IML13_alt)+"\t"+str(QE14_ref)+"\t"+str(QE14_alt)+"\t"+str(QL14_ref)+"\t"+str(QL14_alt)+"\t"+str(QE13_ref)+"\t"+str(QE13_alt)+"\t"+str(QL13_ref)+"\t"+str(QL13_alt)+"\n")
                    cutsites+=1
                else:
                    out0.write(scaff+"\t"+str(position)+"\t"+ref_base+"\t"+alt_base+"\t"+str(BRE13_ref)+"\t"+str(BRE13_alt)+"\t"+str(BRL13_ref)+"\t"+str(BRL13_alt)+"\t"+str(IME14_ref)+"\t"+str(IME14_alt)+"\t"+str(IML14_ref)+"\t"+str(IML14_alt)+"\t"+str(IME13_ref)+"\t"+str(IME13_alt)+"\t"+str(IML13_ref)+"\t"+str(IML13_alt)+"\t"+str(QE14_ref)+"\t"+str(QE14_alt)+"\t"+str(QL14_ref)+"\t"+str(QL14_alt)+"\t"+str(QE13_ref)+"\t"+str(QE13_alt)+"\t"+str(QL13_ref)+"\t"+str(QL13_alt)+"\n")        
                    noncutsites+=1
                out1.write(scaff+"\t"+str(position)+"\t"+str(num_reads)+"\t"+str(BRE13_count)+"\t"+str(BRL13_count)+"\t"+str((BRE13_count+BRL13_count))+"\t"+str(IME14_count)+"\t"+str(IML14_count)+"\t"+str(IME13_count)+"\t"+str(IML13_count)+"\t"+str((IME13_count+IML13_count+IME14_count+IML14_count))+"\t"+str(QE14_count)+"\t"+str(QL14_count)+"\t"+str(QE13_count)+"\t"+str(QL13_count)+"\t"+str((QE13_count+QL13_count+QE14_count+QL14_count))+"\n")
            if line_idx % 100000 == 0:
                print line_idx
    
    out1.close()
    out0.close()
    out10.close()    
    
    print "Mean read depth across samples ",raw_read_count/depth_cc      
    print "number of sites for RN analysis = ", cutsites
    print "number of sites for AF analysis = ", noncutsites                        

if __name__ == "__main__":
    main()

