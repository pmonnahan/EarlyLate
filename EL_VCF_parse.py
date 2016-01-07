#-------------------------------------------------------------------------------
# Name:        VCF parse
# Purpose:     Takes VCF and returns allele counts for each sample in the file outfilepath0.
#Notes: Run the program one time with permissive settings for Min_read_depth and Max_read_depth.  Look at the distribution of coverage per site in the file 'read_depth.txt'.  
#       Use the distribution of coverage, particularly the median, to set the min and max read depth.  
#-------------------------------------------------------------------------------

import sys

#AF = allele frequency analysis; RN = read number analysis

def main(inpath1="black0.txt", infilepath="/Volumes/TOSHIBA EXT/EarlyLate/EL.big.dt.vcf", outfilepath5="/Users/patrick/Documents/Research/EarlyLate/black1.txt", outfilepath4="/Users/patrick/Documents/Research/EarlyLate/out4.txt", outfilepath3="/Users/patrick/Documents/Research/EarlyLate/Ordered_scaff_list.txt", outfilepath2="/Users/patrick/Documents/Research/EarlyLate/genotypes_X.txt", outfilepath0="/Users/patrick/Documents/Research/EarlyLate/counts_dt_SandP_AF250_10-11-15.txt", outfilepath1="/Users/patrick/Documents/Research/EarlyLate/read_depth_dt.txt", outfilepath6="/Users/patrick/Documents/Research/EarlyLate/counts_parsed_8-11-15.txt",outfilepath10="/Users/patrick/Documents/Research/EarlyLate/counts_dt_SandP_RN250_10-11-15.txt"):
    src  =open(infilepath, "rU")
#    in0  =open(inpath1, "rU")
    out0 =open(outfilepath0, "w")
    out1 =open(outfilepath1, "w")
#    out2 =open(outfilepath2, "w")
#    out3 =open(outfilepath3, "w")
#    out4 =open(outfilepath4, "w")
#    out5 =open(outfilepath5, "w")
#    out6 =open(outfilepath6, "w")
#    out9=open("Scaff14_bases.txt","w")
#    out9.write("rownum"+"\t"+"scaff"+"\t"+"pos"+"\t"+"IME13R"+"\t"+"IME13A"+"\t"+"IML13R"+"\t"+"IML13A"+"\t"+"IME14R"+"\t"+"IME14A"+"\t"+"IML14R"+"\t"+"IML14A"+"\t"+"lrt_B_im"+"\t"+"df_B_im"+"\t"+"p_B_im"+"\t"+"ref_base"+"\t"+"alt_base"+"\n")
    out10=open(outfilepath10,"w")
#    with open("/Users/patrick/Documents/Research/EarlyLate/Scaff14.txt","r") as scaff14:
#        info=[]
#        poses=[]
#        for i,posy in enumerate(scaff14):
#            if i >0:
#                posy=posy.strip("\n")
#                posy=posy.strip("\r")
#                posx=posy
#                posy=posy.split("\t")
#                poses.append(int(posy[2]))
#                info.append(posx)
#            
    

#       USER DEFINED
    Number_RILs = 22

    Scaffold = []
    BPpos = []
    LineID= []
    genotypes = [[] for i in range(Number_RILs)]
    BLMarkers = {}
    
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
    cutsites = 0
    noncutsites = 0
    out1.write("Scaff"+"\t"+"pos"+"\t"+"DP"+"\t"+"BRE_14"+"\t"+"BRE_13"+"\t"+"BRL_13"+"\t"+"BR_count"+"\t"+"IME_14"+"\t"+"IML_14"+"\t"+"IME_13"+"\t"+"IML_13"+"\t"+"IM_count"+"\t"+"QE_14"+"\t"+"QL_14"+"\t"+"QE_13"+"\t"+"QL_13"+"\t"+"Q_count"+"\n")
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
            out10.write("Scaff"+"\t"+"pos"+"\t"+"ref"+"\t"+"alt"+"\t")
#            out6.write("Scaff"+"\t"+"pos"+"\t"+"ref"+"\t"+"alt"+"\t")
##            for i in range(len(LineID)):
##                out0.write(LineID[i]+"R"+"\t"+LineID[i]+"A"+"\t")
            out0.write("BRE14_R"+"\t"+"BRE14_A"+"\t"+"BRE13_R"+"\t"+"BRE13_A"+"\t"+"BRL13_R"+"\t"+"BRL13_A"+"\t"+"IME14_R"+"\t"+"IME14_A"+"\t"+"IML14_R"+"\t"+"IML14_A"+"\t"+"IME13_R"+"\t"+"IME13_A"+"\t"+"IML13_R"+"\t"+"IML13_A"+"\t"+"QE14_R"+"\t"+"QE14_A"+"\t"+"QL14_R"+"\t"+"QL14_A"+"\t"+"QE13_R"+"\t"+"QE13_A"+"\t"+"QL13_R"+"\t"+"QL13_A"+"\n")
            out10.write("BRE14_R"+"\t"+"BRE14_A"+"\t"+"BRE13_R"+"\t"+"BRE13_A"+"\t"+"BRL13_R"+"\t"+"BRL13_A"+"\t"+"IME14_R"+"\t"+"IME14_A"+"\t"+"IML14_R"+"\t"+"IML14_A"+"\t"+"IME13_R"+"\t"+"IME13_A"+"\t"+"IML13_R"+"\t"+"IML13_A"+"\t"+"QE14_R"+"\t"+"QE14_A"+"\t"+"QL14_R"+"\t"+"QL14_A"+"\t"+"QE13_R"+"\t"+"QE13_A"+"\t"+"QL13_R"+"\t"+"QL13_A"+"\n")
#            out6.write("BRE14_RP"+"\t"+"BRE14_AP"+"\t"+"BRE14_RS"+"\t"+"BRE14_AS"+"\t"+"BRE13_RP"+"\t"+"BRE13_AP"+"\t"+"BRE13_RS"+"\t"+"BRE13_AS"+"\t"+"BRL13_RP"+"\t"+"BRL13_AP"+"\t"+"BRL13_RS"+"\t"+"BRL13_AS"+"\t"+"IME14_RP"+"\t"+"IME14_AP"+"\t"+"IME14_RS"+"\t"+"IME14_AS"+"\t"+"IML14_RP"+"\t"+"IML14_AP"+"\t"+"IML14_RS"+"\t"+"IML14_AS"+"\t"+"IME13_RP"+"\t"+"IME13_AP"+"\t"+"IME13_RS"+"\t"+"IME13_AS"+"\t"+"IML13_RP"+"\t"+"IML13_AP"+"\t"+"IML13_RS"+"\t"+"IML13_AS"+"\t"+"QE14_RP"+"\t"+"QE14_AP"+"\t"+"QE14_RS"+"\t"+"QE14_AS"+"\t"+"QL14_RP"+"\t"+"QL14_AP"+"\t"+"QL14_RS"+"\t"+"QL14_AS"+"\t"+"QE13_RP"+"\t"+"QE13_AP"+"\t"+"QE13_RS"+"\t"+"QE13_AS"+"\t"+"QL13_RP"+"\t"+"QL13_AP"+"\t"+"QL13_RS"+"\t"+"QL13_AS"+"\n")
        else:
#            if line_idx<2000:
            key1 = cols[0]+"_"+cols[1]
            try:
                y=BLMarkers[key1]
                #print "black listed"
            except (KeyError):
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
                    
                    if cols[9] != "./.":#BR 14
                        basesP = cols[9].split(":") # Early
                        BRE14_refp,BRE14_altp = basesP[1].split(",")
                    else:
                        BRE14_refp=BRE14_altp=0
                    if cols[10] != "./.":                         
                        basesS = cols[10].split(":")
                        BRE14_refs,BRE14_alts = basesS[1].split(",")
                    else:
                        BRE14_refs=BRE14_alts=0
                    if cols[11] != "./." :
                        basesP = cols[11].split(":")  #Early
                        BRE13_refp,BRE13_altp = basesP[1].split(",")
                    else:
                        BRE13_refp=BRE13_altp=0
                    if cols[12] != "./.":                                            
                        basesS = cols[12].split(":")
                        BRE13_refs,BRE13_alts = basesS[1].split(",")
                    else:
                        BRE13_refs=BRE13_alts=0
                    if cols[13] != "./.":
                        basesP = cols[13].split(":") #Late
                        BRL13_refp,BRL13_altp = basesP[1].split(",")
                    else:
                        BRL13_refp=BRL13_altp=0
                    if cols[14] != "./.":
                        basesS = cols[14].split(":")
                        BRL13_refs,BRL13_alts = basesS[1].split(",")
                    else:
                        BRL13_refs=BRL13_alts=0
                    if cols[15] != "./.":#IM 13
                        basesP = cols[15].split(":")  #Early
                        IME13_refp,IME13_altp = basesP[1].split(",")  
                    else:
                        IME13_refp=IME13_altp=0
                    if cols[16] != "./.":
                        basesS = cols[16].split(":")
                        IME13_refs,IME13_alts = basesS[1].split(",")
                    else:
                        IME13_refs=IME13_alts=0
                    if cols[17] !="./.":
                        basesP = cols[17].split(":") #Late
                        IML13_refp,IML13_altp = basesP[1].split(",")
                    else:
                        IML13_refp=IML13_altp=0
                    if cols[18] != "./.":                    
                        basesS = cols[18].split(":")
                        IML13_refs,IML13_alts = basesS[1].split(",")
                    else:
                        IML13_refs=IML13_alts=0    
                    if cols[19] != "./.": #IM 14 
                        basesP = cols[19].split(":")  #Early
                        IME14_refp,IME14_altp = basesP[1].split(",")
                    else:
                        IME14_refp=IME14_altp=0
                    if cols[20] != "./.":
                        basesS = cols[20].split(":")
                        IME14_refs,IME14_alts = basesS[1].split(",")
                    else:
                        IME14_refs=IME14_alts=0
                    if cols[21] != "./.":
                        basesP = cols[21].split(":") #Late
                        IML14_refp,IML14_altp = basesP[1].split(",")
                    else:
                        IML14_refp=IML14_altp=0
                    if cols[22] != "./.":                 
                        basesS = cols[22].split(":")
                        IML14_refs,IML14_alts = basesS[1].split(",")
                    else:
                        IML14_refs=IML14_alts=0
                    if cols[23] != "./.":#Q 14 
                        basesP = cols[23].split(":")  #Early
                        QE14_refp,QE14_altp = basesP[1].split(",")
                    else:
                        QE14_refp=QE14_altp=0                        
                    if cols[24] != "./.":
                        basesS = cols[24].split(":")
                        QE14_refs,QE14_alts = basesS[1].split(",")
                    else:
                        QE14_refs=QE14_alts=0
                    if cols[25] != "./.":
                        basesP = cols[25].split(":") #Late
                        QL14_refp,QL14_altp = basesP[1].split(",")
                    else:
                        QL14_refp=QL14_altp=0
                    if cols[26] != "./.":                   
                        basesS = cols[26].split(":")
                        QL14_refs,QL14_alts = basesS[1].split(",")
                    else:
                        QL14_refs=QL14_alts=0                        
                    if cols[28] != "./.": #Q 13
                        basesP = cols[28].split(":")  #Early
                        QE13_refp,QE13_altp = basesP[1].split(",")
                    else:
                        QE13_refp=QE13_altp=0
                    if cols[27] != "./.":
                        basesS = cols[27].split(":")
                        QE13_refs,QE13_alts = basesS[1].split(",")
                    else:
                        QE13_refs=QE13_alts=0
                    if cols[29] != "./.":
                        basesP = cols[29].split(":") #Late
                        QL13_refp,QL13_altp = basesP[1].split(",")
                    else:
                        QL13_refp=QL13_altp=0
                    if cols[30] != "./.":                    
                        basesS = cols[30].split(":")
                        QL13_refs,QL13_alts = basesS[1].split(",")
                    else:
                        QL13_refs=QL13_alts=0
                        
                    BRE13_ref=int(BRE13_refp)+int(BRE13_refs)
                    BRE13_alt=int(BRE13_altp)+int(BRE13_alts)
                    BRL13_ref=int(BRL13_refp)+int(BRL13_refs)
                    BRL13_alt=int(BRL13_altp)+int(BRL13_alts)
                    BRE14_ref=int(BRE14_refp)+int(BRE14_refs)
                    BRE14_alt=int(BRE14_altp)+int(BRE14_alts)
                    IME13_ref=int(IME13_refp)+int(IME13_refs)
                    IME13_alt=int(IME13_altp)+int(IME13_alts)
                    IML13_ref=int(IML13_refp)+int(IML13_refs)
                    IML13_alt=int(IML13_altp)+int(IML13_alts)
                    IME14_ref=int(IME14_refp)+int(IME14_refs)
                    IME14_alt=int(IME14_altp)+int(IME14_alts)
                    IML14_ref=int(IML14_refp)+int(IML14_refs)
                    IML14_alt=int(IML14_altp)+int(IML14_alts)
                    QE13_ref=int(QE13_refp)+int(QE13_refs)
                    QE13_alt=int(QE13_altp)+int(QE13_alts)
                    QL13_ref=int(QL13_refp)+int(QL13_refs)
                    QL13_alt=int(QL13_altp)+int(QL13_alts)
                    QE14_ref=int(QE14_refp)+int(QE14_refs)
                    QE14_alt=int(QE14_altp)+int(QE14_alts)
                    QL14_ref=int(QL14_refp)+int(QL14_refs)
                    QL14_alt=int(QL14_altp)+int(QL14_alts)
                    
                    BRE14_count=BRE14_ref+BRE14_alt
                    BRE13_count=BRE13_ref+BRE13_alt
                    BRL13_count=BRL13_ref+BRL13_alt
                    IME14_count=IME14_ref+IME14_alt
                    IME13_count=IME13_ref+IME13_alt
                    IML13_count=IML13_ref+IML13_alt
                    IML14_count=IML14_ref+IML14_alt
                    QE14_count=QE14_ref+QE14_alt
                    QE13_count=QE13_ref+QE13_alt
                    QL13_count=QL13_ref+QL13_alt
                    QL14_count=QL14_ref+QL14_alt
                    
                    if BRE14_count > 250 or BRE13_count > 250 or BRL13_count > 250 or IME13_count > 250 or IML13_count > 250 or IME14_count > 250 or IML14_count > 250 or QE13_count > 250 or QL13_count > 250 or QE14_count > 250 or QL14_count > 250:
                        out10.write(cols[0]+"\t"+cols[1]+"\t"+cols[3]+"\t"+cols[4]+"\t"+str(BRE14_ref)+"\t"+str(BRE14_alt)+"\t"+str(BRE13_ref)+"\t"+str(BRE13_alt)+"\t"+str(BRL13_ref)+"\t"+str(BRL13_alt)+"\t"+str(IME14_ref)+"\t"+str(IME14_alt)+"\t"+str(IML14_ref)+"\t"+str(IML14_alt)+"\t"+str(IME13_ref)+"\t"+str(IME13_alt)+"\t"+str(IML13_ref)+"\t"+str(IML13_alt)+"\t"+str(QE14_ref)+"\t"+str(QE14_alt)+"\t"+str(QL14_ref)+"\t"+str(QL14_alt)+"\t"+str(QE13_ref)+"\t"+str(QE13_alt)+"\t"+str(QL13_ref)+"\t"+str(QL13_alt)+"\n")        
                        cutsites+=1
                    else:
                        out0.write(cols[0]+"\t"+cols[1]+"\t"+cols[3]+"\t"+cols[4]+"\t"+str(BRE14_ref)+"\t"+str(BRE14_alt)+"\t"+str(BRE13_ref)+"\t"+str(BRE13_alt)+"\t"+str(BRL13_ref)+"\t"+str(BRL13_alt)+"\t"+str(IME14_ref)+"\t"+str(IME14_alt)+"\t"+str(IML14_ref)+"\t"+str(IML14_alt)+"\t"+str(IME13_ref)+"\t"+str(IME13_alt)+"\t"+str(IML13_ref)+"\t"+str(IML13_alt)+"\t"+str(QE14_ref)+"\t"+str(QE14_alt)+"\t"+str(QL14_ref)+"\t"+str(QL14_alt)+"\t"+str(QE13_ref)+"\t"+str(QE13_alt)+"\t"+str(QL13_ref)+"\t"+str(QL13_alt)+"\n")        
                        noncutsites+=1
#                    for jj, scaff14 in enumerate(poses):
#                        if position == scaff14 and scaff == "sNNffold_14":
#                            print "here1"
#                            out9.write(info[jj]+"\t"+ref_base+"\t"+alt_base+"\n")
                    
                    
#                    out6.write(cols[0]+"\t"+cols[1]+"\t"+cols[3]+"\t"+cols[4]+"\t"+str(BRE14_refp)+"\t"+str(BRE14_altp)+"\t"+str(BRE14_refs)+"\t"+str(BRE14_alts)+"\t"+str(BRE13_refp)+"\t"+str(BRE13_altp)+"\t"+str(BRE13_refs)+"\t"+str(BRE13_alts)+"\t"+str(BRL13_refp)+"\t"+str(BRL13_altp)+"\t"+str(BRL13_refs)+"\t"+str(BRL13_alts)+"\t"+str(IME14_refp)+"\t"+str(IME14_altp)+"\t"+str(IME14_refs)+"\t"+str(IME14_alts)+"\t"+str(IML14_refp)+"\t"+str(IML14_altp)+"\t"+str(IML14_refs)+"\t"+str(IML14_alts)+"\t"+str(IME13_refp)+"\t"+str(IME13_altp)+"\t"+str(IME13_refs)+"\t"+str(IME13_alts)+"\t"+str(IML13_refp)+"\t"+str(IML13_altp)+"\t"+str(IML13_refs)+"\t"+str(IML13_alts)+"\t"+str(QE14_refp)+"\t"+str(QE14_altp)+"\t"+str(QE14_refs)+"\t"+str(QE14_alts)+"\t"+str(QL14_refp)+"\t"+str(QL14_altp)+"\t"+str(QL14_refs)+"\t"+str(QL14_alts)+"\t"+str(QE13_refp)+"\t"+str(QE13_altp)+"\t"+str(QE13_refs)+"\t"+str(QE13_alts)+"\t"+str(QL13_refp)+"\t"+str(QL13_altp)+"\t"+str(QL13_refs)+"\t"+str(QL13_alts)+"\n")
                    out1.write(str(scaff)+"\t"+str(position)+"\t"+str(num_reads)+"\t"+str(BRE14_count)+"\t"+str(BRE13_count)+"\t"+str(BRL13_count)+"\t"+str((BRE14_count+BRE13_count+BRL13_count))+"\t"+str(IME14_count)+"\t"+str(IML14_count)+"\t"+str(IME13_count)+"\t"+str(IML13_count)+"\t"+str((IME13_count+IML13_count+IME14_count+IML14_count))+"\t"+str(QE14_count)+"\t"+str(QL14_count)+"\t"+str(QE13_count)+"\t"+str(QL13_count)+"\t"+str((QE13_count+QL13_count+QE14_count+QL14_count))+"\n")
                if line_idx % 100000 == 0:
                    print line_idx
        
    
#    print "Mean read depth across samples ",raw_read_count/depth_cc      
#    print "Valid snps ",valid_snps
    print "number of sites for RN analysis = ", cutsites
    print "number of sites for AF analysis = ", noncutsites                        

if __name__ == "__main__":
    main()

