# matching Patrick for EL data 12/7/15
# 2	IM pops
import math
DIR="/Users/patrick/Documents/Research/EarlyLate/"
src  =open(DIR+"counts_massiveJan16_AF250.txt", "rU")
out1 =open(DIR+"IM.cc.txt","w")
out2 =open(DIR+"Q.cc.txt","w")


# evaluate contents of each line of input file
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t') 
    out1.write(cols[0]+'\t'+cols[1])
    out2.write(cols[0]+'\t'+cols[1])
    if line_idx==0:
        
        for k in range(16,24): #Adjust these numbers to get correct counts for a population  #Quarry is 16 through 23
            out2.write('\t'+cols[k])
            
        for k in range(8,16):
            out1.write('\t'+cols[k]) #IM is 8 through 15
    else:
        for k in range(16,24):
            out2.write('\t'+cols[k])
        for k in range(8,16):
            out1.write('\t'+cols[k])
    out1.write('\n')
    out2.write('\n')
    
    if line_idx % 100000 == 0:
        print line_idx


