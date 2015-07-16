#-------------------------------------------------------------------------------
# Name:        VCF parse
# Purpose:     Takes VCF and returns allele counts for each sample in the file outfilepath0.
#Notes: Run the program one time with permissive settings for Min_read_depth and Max_read_depth.  Look at the distribution of coverage per site in the file 'read_depth.txt'.  
#       Use the distribution of coverage, particularly the median, to set the min and max read depth.  
#-------------------------------------------------------------------------------

import sys

def main(inpath1="black0.txt", infilepath="Early_Late_all.vcf", outfilepath5="black1.txt", outfilepath4="out4.txt", outfilepath3="Ordered_scaff_list.txt", outfilepath2="genotypes_X.txt", outfilepath0="counts_bothYears_5-26-15.txt", outfilepath1="read_depth.txt"):
    src  =open(infilepath, "rU")
#    in0  =open(inpath1, "rU")
    out0 =open(outfilepath0, "w")
    out1 =open(outfilepath1, "w")
    out2 =open(outfilepath2, "w")
    out3 =open(outfilepath3, "w")
    out4 =open(outfilepath4, "w")
    out5 =open(outfilepath5, "w")

#       USER DEFINED
    Min_read_depth = 5 # per sample
    Max_read_depth = 200
    Number_RILs = 11
    RareQ = 0.00

    Scaffold = []
    BPpos = []
    LineID= []
    genotypes = [[] for i in range(Number_RILs)]
    BLMarkers = {}

        # read black listed markers
#        for line_idx, line in enumerate(in0):
#                cols = line.replace('\n', '').split('\t') 
#        key1 = cols[0]+"_"+cols[1]
#        BLMarkers[key1]=1

#       Counters
    depth_cc=0
    raw_read_count=0
    valid_snps=0
    newlines=[]
    line_num=0
    out1.write("DP"+"\t"+"BRE_14"+"\t"+"BRE_13"+"\t"+"BRL_13"+"\t"+"BR_count"+"\t"+"IME_14"+"\t"+"IML_14"+"\t"+"IME_13"+"\t"+"IML_13"+"\t"+"IM_count"+"\t"+"QE_14"+"\t"+"QL_14"+"\t"+"QE_13"+"\t"+"QL_13"+"\t"+"Q_count"+"\n")
        # evaluate contents of each line of input file
    for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t') 
        # print line_idx
        if len(cols) < 2:               ## This should be header
            continue
        elif cols[0] == "#CHROM":
            for i in range(len(cols)):
                print i,cols[i]
                subcols = cols[i].split('.')
                if len(subcols) > 1 and i>8: 
                    LineID.append(subcols[0])
            out0.write("Scaff"+"\t"+"pos"+"\t"+"ref"+"\t"+"alt"+"\t")
#            for i in range(len(LineID)):
#                out0.write(LineID[i]+"R"+"\t"+LineID[i]+"A"+"\t")
            out0.write("BRE13_R"+"\t"+"BRE13_A"+"\t"+"BRL13_R"+"\t"+"BRL13_A"+"\t"+"BRE14_R"+"\t"+"BRE14_A"+"\t"+"IME13_R"+"\t"+"IME13_A"+"\t"+"IML13_R"+"\t"+"IML13_A"+"\t"+"IME14_R"+"\t"+"IME14_A"+"\t"+"IML14_R"+"\t"+"IML14_A"+"\t"+"QE13_R"+"\t"+"QE13_A"+"\t"+"QL13_R"+"\t"+"QL13_A"+"\t"+"QE14_R"+"\t"+"QE14_A"+"\t"+"QL14_R"+"\t"+"QL14_A"+"\t"+"BRE_R"+"\t"+"BRE_A"+"\t"+"BRL_R"+"\t"+"BRL_A"+"\t"+"IME_R"+"\t"+"IME_A"+"\t"+"IML_R"+"\t"+"IML_A"+"\t"+"QE_R"+"\t"+"QE_A"+"\t"+"QL_R"+"\t"+"QL_A"+"\n")

        else:
#            if line_idx<2000:
            key1 = cols[0]+"_"+cols[1]
            try:
                y=BLMarkers[key1]
                #print "black listed"
            except (KeyError):                        
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
                    
                    pop_count=0
                    if cols[10] != "./." and cols[11] != "./.":
                        bases = cols[10].split(":")
                        #print bases[1]
                        BRE13_ref,BRE13_alt = bases[1].split(",")
                        bases = cols[11].split(":")
                        #print bases
                        BRL13_ref,BRL13_alt = bases[1].split(",")
                        BRE13_count= int(BRE13_ref)+int(BRE13_alt)
                        BRL13_count= int(BRL13_ref)+int(BRL13_alt)
                        qtot = (float(BRE13_ref)+float(BRL13_ref))/float(float(BRE13_count)+float(BRL13_count))
                        if BRE13_count >= Min_read_depth and BRL13_count >= Min_read_depth and BRE13_count <= Max_read_depth and BRL13_count <= Max_read_depth and qtot>=RareQ and qtot<=(1.0-RareQ):
                            pop_count+=1
                            BR13_ref_count=int(BRE13_ref)+int(BRL13_ref)
                            BR13_alt_count=int(BRE13_alt)+int(BRL13_alt)
                            BR13Info=(BRE13_ref+"\t"+BRE13_alt+"\t"+BRL13_ref+"\t"+BRL13_alt+"\t")
                        else:
                            BR13Info=("0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t")
                    
                    if cols[9] != "./.":
                        bases = cols[9].split(":")
                        #print bases[1]
                        BRE14_ref,BRE14_alt = bases[1].split(",")
                        BRE14_count= int(BRE14_ref)+int(BRE14_alt)
                        qtot = (float(BRE14_ref))/float(BRE14_count)
                        if BRE14_count >= Min_read_depth and BRE14_count <= Max_read_depth and qtot>=RareQ and qtot<=(1.0-RareQ):
                            pop_count+=1
                            BR14_ref_count=int(BRE14_ref)
                            BR14_alt_count=int(BRE14_alt)
                            BR14Info=(BRE14_ref+"\t"+BRE14_alt+"\t")
                        else:
                            BR14Info=("0"+"\t"+"0"+"\t")
                    
                    if cols[12] != "./." and cols[14] != "./.":
                        bases = cols[12].split(":")
                        IME14ref,IME14alt = bases[1].split(",")
                        bases = cols[14].split(":")
                        IML14ref,IML14alt = bases[1].split(",")
                        IME14_count= int(IME14ref)+int(IME14alt)
                        IML14_count= int(IML14ref)+int(IML14alt)
                        qtot = (float(IME14ref)+float(IML14ref))/float(float(IME14_count)+float(IML14_count))
                        if IME14_count >= Min_read_depth and IML14_count >= Min_read_depth and IME14_count <= Max_read_depth and IML14_count <= Max_read_depth and qtot>RareQ and qtot<(1.0-RareQ):
                            IM14_Info=(IME14ref+"\t"+IME14alt+"\t"+IML14ref+"\t"+IML14alt+"\t")
                            IM14_ref_count=int(IME14ref)+int(IML14ref)
                            IM14_alt_count=int(IME14alt)+int(IML14alt)
                            pop_count+=1
                        else:
                            IM14_ref_count=0
                            IM14_alt_count=0                         
                            IM14_Info=("0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t")
                        
                    if cols[13] != "./." and cols[15] != "./.":
                        bases = cols[13].split(":")
                        IME13ref,IME13alt = bases[1].split(",")
                        bases = cols[15].split(":")
                        IML13ref,IML13alt = bases[1].split(",")
                        IME13_count= int(IME13ref)+int(IME13alt)
                        IML13_count= int(IML13ref)+int(IML13alt)
                        qtot = (float(IME13ref)+float(IML13ref))/float(float(IME13_count)+float(IML13_count))
                        if IME13_count >= Min_read_depth and IML13_count >= Min_read_depth and IME13_count <= Max_read_depth and IML13_count <= Max_read_depth and qtot>=RareQ and qtot<=(1.0-RareQ):
                            IM13_Info=(IME13ref+"\t"+IME13alt+"\t"+IML13ref+"\t"+IML13alt+"\t")
                            IM13_ref_count=int(IME13ref)+int(IML13ref)
                            IM13_alt_count=int(IME13alt)+int(IML13alt)
                            pop_count+=1
                        else:
                            IM13_ref_count=0
                            IM13_alt_count=0                         
                            IM13_Info=("0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t") 
                            
                    if cols[16] != "./." and cols[17] != "./.":
                         bases = cols[16].split(":")
                         QE14ref,QE14alt = bases[1].split(",")
                         bases = cols[17].split(":")
                         QL14ref,QL14alt = bases[1].split(",")
                         QE14_count= int(QE14ref)+int(QE14alt)
                         QL14_count= int(QL14ref)+int(QL14alt)
                         qtot = (float(QE14ref)+float(QL14ref))/float(float(QE14_count)+float(QL14_count))
                         if QE14_count >= Min_read_depth and QL14_count >= Min_read_depth and QE14_count <= Max_read_depth and QL14_count <= Max_read_depth and qtot>=RareQ and qtot<=(1.0-RareQ):
                             pop_count+=1
                             Q14_ref_count=int(QE14ref)+int(QL14ref)
                             Q14_alt_count=int(QE14alt)+int(QL14alt)
                             Q14_Info=(QE14ref+"\t"+QE14alt+"\t"+QL14ref+"\t"+QL14alt+"\t")
                         else:
                             Q14_ref_count=0
                             Q14_alt_count=0                            
                             Q14_Info=("0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t") 
                             
                    if cols[18] != "./." and cols[19] != "./.":
                         bases = cols[18].split(":")
                         QE13ref,QE13alt = bases[1].split(",")
                         bases = cols[19].split(":")
                         QL13ref,QL13alt = bases[1].split(",")
                         QE13_count= int(QE13ref)+int(QE13alt)
                         QL13_count= int(QL13ref)+int(QL13alt)
                         qtot = (float(QE13ref)+float(QL13ref))/float(float(QE13_count)+float(QL13_count))
                         if QE13_count >= Min_read_depth and QL13_count >= Min_read_depth and QE13_count <= Max_read_depth and QL13_count <= Max_read_depth and qtot>=RareQ and qtot<=(1.0-RareQ):
                             pop_count+=1
                             Q13_ref_count=int(QE13ref)+int(QL13ref)
                             Q13_alt_count=int(QE13alt)+int(QL13alt)
                             Q13_Info=(QE13ref+"\t"+QE13alt+"\t"+QL13ref+"\t"+QL13alt+"\t")
                         else:
                             Q13_ref_count=0
                             Q13_alt_count=0                            
                             Q13_Info=("0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t")
                             
                    if cols[10] != "./." and cols[11] != "./." or cols[9] != "./.":
                        BRE_ref_count=int(BRE13_ref)+int(BRE14_ref)
                        BRE_alt_count=int(BRE13_alt)+int(BRE14_alt)
                        BRL_ref_count=BRL13_ref
                        BRL_alt_count=BRL13_alt
                        BRInfo=(str(BRE_ref_count)+"\t"+str(BRE_alt_count)+"\t"+str(BRL_ref_count)+"\t"+str(BRL_alt_count)+"\t")
                    else:
                        BRInfo=("0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t")
                    if (cols[16] != "./." and cols[17] != "./.") or (cols[18] != "./." and cols[19] != "./."):
                        QE_ref_count=int(QE13ref)+int(QE14ref)
                        QE_alt_count=int(QE13alt)+int(QE14alt)
                        QL_ref_count=int(QL13ref)+int(QL14ref)
                        QL_alt_count=int(QL13alt)+int(QL14alt)
                        QInfo=(str(QE_ref_count)+"\t"+str(QE_alt_count)+"\t"+str(QL_ref_count)+"\t"+str(QL_alt_count))
                    else:
                        QInfo=("0"+"\t"+"0"+"\t"+"0"+"\t"+"0")
                        
                    
                    if (cols[12] != "./." and cols[14] != "./.") or (cols[13] != "./." and cols[15] != "./."):
                        IME_ref_count=int(IME13ref)+int(IME14ref)
                        IME_alt_count=int(IME13alt)+int(IME14alt)
                        IML_ref_count=int(IML13ref)+int(IML14ref)
                        IML_alt_count=int(IML13alt)+int(IML14alt)
                        IMInfo=(str(IME_ref_count)+"\t"+str(IME_alt_count)+"\t"+str(IML_ref_count)+"\t"+str(IML_alt_count)+"\t")
                    else:
                        IMInfo=("0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t")
                    
                    
                    
                    if pop_count>0:
                        if (Q13_ref_count+Q14_ref_count+IM13_ref_count+IM14_ref_count+BR13_ref_count+BR14_ref_count)!=0 and (Q13_alt_count+Q14_alt_count+IM13_alt_count+IM14_alt_count+BR13_alt_count+BR14_alt_count)!=0:
#                                out0.write("BRE13_R"+"\t"+"BRE13_A"+"\t"+"BRE14_R"+"\t"+"BRE14_A"+"\t"+"BRL13_R"+"\t"+"BRL13_A"+"\t"+"IME14_R"+"\t"+"IME14_A"+"\t"+"IML14_R"+"\t"+"IML14_A"+"\t"+"IME13_R"+"\t"+"IME13_A"+"\t"+"IML13_R"+"\t"+"IML13_A"+"\t"+"QE14_R"+"\t"+"QE14_A"+"\t"+"QL14_R"+"\t"+"QL14_A"+"\t"+"QE13_R"+"\t"+"QE13_A"+"\t"+"QL13_R"+"\t"+"QL13_A"+"\t"+"BRE_R"+"\t"+"BRE_A"+"\t"+"BRL_R"+"\t"+"BRL_A"+"\t"+"IME_R"+"\t"+"IME_A"+"\t"+"IML_R"+"\t"+"IML_A"+"\t"+"QE_R"+"\t"+"QE_A"+"\t"+"QL_R"+"\t"+"QL_A"+"\n")
                            out0.write(cols[0]+"\t"+cols[1]+"\t"+cols[3]+"\t"+cols[4]+"\t"+BR13Info+BR14Info+IM13_Info+IM14_Info+Q13_Info+Q14_Info+BRInfo+IMInfo+QInfo+"\n")        
                            out1.write(str(num_reads)+"\t"+str(BRE14_count)+"\t"+str(BRE13_count)+"\t"+str(BRL13_count)+"\t"+str((BRE14_count+BRE13_count+BRL13_count))+"\t"+str(IME14_count)+"\t"+str(IML14_count)+"\t"+str(IME13_count)+"\t"+str(IML13_count)+"\t"+str((IME13_count+IML13_count+IME14_count+IML14_count))+"\t"+str(QE14_count)+"\t"+str(QL14_count)+"\t"+str(QE13_count)+"\t"+str(QL13_count)+"\t"+str((QE13_count+QL13_count+QE14_count+QL14_count))+"\n")
                if line_idx % 100000 == 0:
                    print line_idx
        
    
#    print "Mean read depth across samples ",raw_read_count/depth_cc      
#    print "Valid snps ",valid_snps                        

if __name__ == "__main__":
    main()

