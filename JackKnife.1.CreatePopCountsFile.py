# matching Patrick for EL data 12/7/15
# 2	IM pops
import math

src  =open("counts_dt_SandP_AF250_10-11-15.txt", "rU")
#out1 =open("IM.snps.txt","w")
out2 =open("Q.cc.txt","w")

# evaluate contents of each line of input file
for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t') 
	out1.write(cols[0]+'\t'+cols[1])
	out2.write(cols[0]+'\t'+cols[1])
	if line_idx==0:
# Scaff	pos	ref	alt	BRE14_R	BRE14_A	BRE13_R	BRE13_A	BRL13_R	BRL13_A	
# IME14_R	IME14_A	IML14_R	IML14_A	IME13_R	IME13_A	IML13_R	IML13_A	QE14_R	QE14_A	QL14_R	QL14_A	QE13_R	QE13_A	QL13_R	QL13_A
		Samples=[]

		for k in range(18,26): #Adjust these numbers to get correct counts for a population 
			Samples.append(cols[k]) #Quarry is 18 through 26
			out1.write('\t'+cols[k]) #IM is 10 through 17
			out2.write('\t'+cols[k])
		
	else:
		rawCC=[]
		for k in range(18,26):
			rawCC.append(int(cols[k]))
			out2.write('\t'+cols[k])
		for k in range(4):
			m=rawCC[2*k]+rawCC[2*k+1]
			if m>0:
				q=float(rawCC[2*k])/float(m)
				z=2.0*math.asin(q**0.5)
			else:
				z=-9
			out1.write('\t'+str(m)+'\t'+str(z))

	out1.write('\n')
	out2.write('\n')

        if line_idx % 100000 == 0:
            print line_idx


