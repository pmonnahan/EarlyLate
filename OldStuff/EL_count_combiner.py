# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 20:37:44 2014

@author: patrick
"""
import csv
from matplotlib_venn import *
from matplotlib import pyplot as plt
import numpy as np

BR=open("/Users/patrick/Google Drive/counts_EL_Br_min10_max100.txt")
IM=open("/Users/patrick/Google Drive/counts_EL_IM_min10_max100.txt")
Q=open("/Users/patrick/Google Drive/counts_EL_Quarry_min10_max100.txt")

#Determine number of lines in text files.
BR_lines=sum(1 for line in open("/Users/patrick/Google Drive/counts_EL_Br_min10_max100.txt"))
IM_lines=sum(1 for line in open("/Users/patrick/Google Drive/counts_EL_IM_min10_max100.txt"))
Q_lines=sum(1 for line in open("/Users/patrick/Google Drive/counts_EL_Quarry_min10_max100.txt"))

BR_info_list=[]
IM_info_list=[]
Q_info_list=[]
All_info_list=[]

for i in range(BR_lines):
    #BR_info_list.append(BR.readline().split("\t"))
    x=BR.readline().split("\t")
    x.append("BR")
    if len(x)>6:
        All_info_list.append(x)
for i in range(IM_lines):
    #IM_info_list.append(IM.readline().split("\t"))
    x=IM.readline().split("\t")
    x.append("IM")
    if len(x)>6:
        All_info_list.append(x)
for i in range(Q_lines):
    #Q_info_list.append(Q.readline().split("\t"))
    x=Q.readline().split("\t")
    x.append("Q")
    if len(x)>6:
        All_info_list.append(x)
  
print "IM sites = ", IM_lines
print "BR sites = ", BR_lines
print "Q sites = ", Q_lines
print "All sites = ", (len(All_info_list))
  
All_info_list.sort(key=lambda g: (int(g[1]),g[0],g[8])) 


All_three=0
IM_BR=0
IM_Q=0
BR_Q=0
BR_only=0
IM_only=0
Q_only=0
total=0
missed=0

with open("/Users/patrick/Documents/EL_combined_counts.csv","wb") as combined:
    comb_writer=csv.writer(combined,delimiter=",",dialect='excel')
    x=["Chrom","Position","Ref Base","IM Alt Base","IM-E Ref Count", "IM-E Alt Count","IM-L Ref Count","IM-L Alt Count","BR Alt Base", "BR-E Ref Count","BR-E Alt Count", "BR-L Ref Count", "BR-L Alt Count","Q Alt Base","Q-E Ref Count","Q-E Alt Count", "Q-L Ref Count","Q-L Alt Count"]     
    comb_writer.writerow(x)    
    for i,site in enumerate(All_info_list[0:len(All_info_list)-2]):
        if i>0 and site[0]==All_info_list[i-1][0] and int(site[1])==int(All_info_list[i-1][1]):
            pass
        else:
            site_info=[""]*18
            if site[0]==All_info_list[i+1][0] and site[0]==All_info_list[i+2][0] and int(site[1])==int(All_info_list[i+1][1]) and int(site[1])==int(All_info_list[i+2][1]):
                All_three+=1
                #list order: Chr, Pos, Ref Base,IM Alt Base, IM ref Early Count, IM Early alt count, IM Late ref count, IM Late al count, BR alt base, BR early ref count, etc.
                site_info[0]=site[0]
                site_info[1]=site[1]
                site_info[2]=site[2]
                if site[8]=="IM":
                    site_info[3]=site[3]
                    site_info[4]=site[4]
                    site_info[5]=site[5]
                    site_info[6]=site[6]
                    site_info[7]=site[7]
                if All_info_list[i+1][8]=="IM":
                    site_info[3]=All_info_list[i+1][3]
                    site_info[4]=All_info_list[i+1][4]
                    site_info[5]=All_info_list[i+1][5]
                    site_info[6]=All_info_list[i+1][6]
                    site_info[7]=All_info_list[i+1][7]
                if All_info_list[i+2][8]=="IM":
                    site_info[3]=All_info_list[i+2][3]
                    site_info[4]=All_info_list[i+2][4]
                    site_info[5]=All_info_list[i+2][5] 
                    site_info[6]=All_info_list[i+2][6]
                    site_info[7]=All_info_list[i+2][7]
                if site[8]=='BR':
                    site_info[8]=site[3]
                    site_info[9]=site[4]
                    site_info[10]=site[5]
                    site_info[11]=site[6]
                    site_info[12]=site[7]                    
                if All_info_list[i+1][8]=="BR":
                    site_info[8]=All_info_list[i+1][3]
                    site_info[9]=All_info_list[i+1][4]
                    site_info[10]=All_info_list[i+1][5]
                    site_info[11]=All_info_list[i+1][6]
                    site_info[12]=All_info_list[i+1][7]
                if All_info_list[i+2][8]=="BR":
                    site_info[8]=All_info_list[i+2][3]
                    site_info[9]=All_info_list[i+2][4]
                    site_info[10]=All_info_list[i+2][5]
                    site_info[11]=All_info_list[i+2][6]
                    site_info[12]=All_info_list[i+2][7]
                if site[8]=="Q":
                    site_info[13]=site[3]
                    site_info[14]=site[4]
                    site_info[15]=site[5]
                    site_info[16]=site[6]
                    site_info[17]=site[7]
                if All_info_list[i+1][8]=="Q":
                    site_info[13]=All_info_list[i+1][3]
                    site_info[14]=All_info_list[i+1][4]
                    site_info[15]=All_info_list[i+1][5]
                    site_info[16]=All_info_list[i+1][6]
                    site_info[17]=All_info_list[i+1][7]
                if All_info_list[i+2][8]=="Q":
                    site_info[13]=All_info_list[i+2][3]
                    site_info[14]=All_info_list[i+2][4]
                    site_info[15]=All_info_list[i+2][5]
                    site_info[16]=All_info_list[i+2][6]
                    site_info[17]=All_info_list[i+2][7]
                else:
                    print("what")
                comb_writer.writerow(site_info)
            elif site[0]==All_info_list[i+1][0] and int(site[1])==int(All_info_list[i+1][1]):
                missed+=1
                site_info[0]=site[0]
                site_info[1]=site[1]
                site_info[2]=site[2]
                if site[8]=="IM":
                    site_info[3]=site[3]
                    site_info[4]=site[4]
                    site_info[5]=site[5]
                    site_info[6]=site[6]
                    site_info[7]=site[7]
                    if All_info_list[i+1][8]=="BR":
                        IM_BR+=1
                    else:
                        IM_Q+=1
                if All_info_list[i+1][8]=="IM":
                    site_info[3]=All_info_list[i+1][3]
                    site_info[4]=All_info_list[i+1][4]
                    site_info[5]=All_info_list[i+1][5]
                    site_info[6]=All_info_list[i+1][6]
                    site_info[7]=All_info_list[i+1][7]            
                if site[8]=="BR":
                    site_info[8]=site[3]
                    site_info[9]=site[4]
                    site_info[10]=site[5]
                    site_info[11]=site[6]
                    site_info[12]=site[7]
                    if All_info_list[i+1][8]=="IM":
                        IM_BR+=1
                    else:
                        BR_Q+=1                        
                if All_info_list[i+1][8]=="BR":
                    site_info[8]=All_info_list[i+1][3]
                    site_info[9]=All_info_list[i+1][4]
                    site_info[10]=All_info_list[i+1][5]
                    site_info[11]=All_info_list[i+1][6]
                    site_info[12]=All_info_list[i+1][7]  
                if site[8]=="Q":
                    site_info[13]=site[3]
                    site_info[14]=site[4]
                    site_info[15]=site[5]
                    site_info[16]=site[6]
                    site_info[17]=site[7]
                    if All_info_list[i+1][8]=="BR":
                        BR_Q+=1
                    else:
                        IM_Q+=1                        
                if All_info_list[i+1][8]=="Q":
                    site_info[13]=All_info_list[i+1][3]
                    site_info[14]=All_info_list[i+1][4]
                    site_info[15]=All_info_list[i+1][5]
                    site_info[16]=All_info_list[i+1][6]
                    site_info[17]=All_info_list[i+1][7]
                comb_writer.writerow(site_info)
            else:
                #missed+=1
                site_info[0]=site[0]
                site_info[1]=site[1]
                site_info[2]=site[2]
                if site[8]=="IM":
                    site_info[3]=site[3]
                    site_info[4]=site[4]
                    site_info[5]=site[5]
                    site_info[6]=site[6]
                    site_info[7]=site[7]
                    IM_only+=1
                if site[8]=="BR":
                    site_info[8]=site[3]
                    site_info[9]=site[4]
                    site_info[10]=site[5]
                    site_info[11]=site[6]
                    site_info[12]=site[7]
                    BR_only+=1
                if site[8]=="Q":
                    site_info[13]=site[3]
                    site_info[14]=site[4]
                    site_info[15]=site[5]
                    site_info[16]=site[6]
                    site_info[17]=site[7]
                    Q_only+=1
                    
                comb_writer.writerow(site_info)

print "All three = ", All_three, "(Prop IM = ", (float(All_three)/float(IM_lines)), ", Prop BR = ", (float(All_three)/float(BR_lines)), ", Prop Q = ", (float(All_three)/float(Q_lines)), ")"
print "IM_BR = ", IM_BR, "(Prop IM = ", (float(IM_BR)/float(IM_lines)), ", Prop BR = ", (float(IM_BR)/float(BR_lines)),")"
print "IM_Q = ", IM_Q, "(Prop IM = ", (float(IM_Q)/float(IM_lines)), ", Prop Q = ", (float(IM_Q)/float(Q_lines)),")"
print "BR_Q = ", BR_Q, "(Prop BR = ", (float(BR_Q)/float(BR_lines)), ", Prop Q = ", (float(BR_Q)/float(Q_lines)),")"
print "IM_only = ", IM_only, "(Prop = ", float(IM_only)/float(IM_lines), ")"
print "BR_only = ", BR_only, "(Prop = ", float(BR_only)/float(BR_lines), ")"
print "Q_only = ", Q_only, "(Prop = ", float(Q_only)/float(Q_lines), ")"
print "Total sites Accounted for = ", (All_three*3+IM_BR*2+IM_Q*2+BR_Q*2+IM_only+BR_only+Q_only)

plt.figure(figsize=(6,6))
v = venn3(subsets=(IM_only, BR_only, IM_BR, Q_only, IM_Q, BR_Q,All_three), set_labels = ('IM', 'BR', 'Q'))