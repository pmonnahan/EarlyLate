# 
# Second of three programs whose purpose is to JackKnife the contrasts between
# Takes a file containing counts 
# 2a	take output from 2: IM.snps.txt and perform contrasts


import math
DIR="/Users/patrick/Documents/Research/EarlyLate/"
src  =open(DIR+"Q.cc.txt", "rU")
out1 =open(DIR+"Q.dz.minDP25.txt","w")
min_m= 25 #Minimum coverage
minQ=  0.05 #Min allele frequency
maxQ=  0.95

# evaluate contents of each line of input file
for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t') 


	if line_idx==0:

# Scaff	pos	IME14_R	IME14_A	IML14_R	IML14_A	IME13_R	IME13_A	IML13_R	IML13_A
# sNNffold_1	2812	53	1	55	7	55	3	84	5

		pass

	else:
		R=[]
		A=[]
		m=[]
		for k in range(4):
			R.append(int(cols[2*k+2]))
			A.append(int(cols[2*k+3]))
			m.append(int(cols[2*k+2])+int(cols[2*k+3]))
		outstring=cols[0]+'\t'+cols[1]
		for cx in range(3):
			for cy in range(cx+1,4):
				if m[cx]>=min_m and m[cy]>=min_m:
					p1=float(R[cx])/float(m[cx])
					p2=float(R[cy])/float(m[cy])
					if (p1+p2)/2.0 > minQ and (p1+p2)/2.0 < maxQ:
						z1=2.0*math.asin(p1**0.5)
						z2=2.0*math.asin(p2**0.5)
						outstring+='\t'+str(z1-z2)+'\t'+str(1.0/float(m[cx])+1.0/float(m[cy]))
					
					else:
						outstring+='\t-9\t-9'
				else:
					outstring+='\t-9\t-9'

		out1.write(outstring+'\n')

        if line_idx % 100000 == 0:
            print line_idx


