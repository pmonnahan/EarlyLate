#-------------------------------------------------------------------------------
# Name:        VCF parse v2 :: works with AD
# Purpose:     Takes VCF and calls genotypes
# v2:  [TA versus TB] 
#-------------------------------------------------------------------------------

#from __future__ import division
import sys
import math
import random

if len(sys.argv) == 2: # argument should be runID

	S = int(sys.argv[1]) # number of snps per window for B

	inpath= "/abyss/Trichome_project_2015/"

# def main(infilepath="s1_snp.vcf", outfilepath0="Results.x.txt", outfilepath1="SamplesIM.txt"):
        src  =open(inpath+"Ali.vcf", "rU")
	in1  =open(inpath+"chisq.txt", "rU")

        out0 =open(inpath+"Results_TA_v_TB.txt", "w")
        out1 =open(inpath+"Vz"+sys.argv[1]+"_TA_v_TB.txt", "w")
	out2 =open(inpath+"z"+sys.argv[1]+"_TA_v_TB.txt","w")
	out3 =open(inpath+"B"+sys.argv[1]+"_TA_v_TB.txt","w")
	out4 =open(inpath+"zSpecial"+sys.argv[1]+"_TA_v_TB.txt","w")

#       USER DEFINED

	
	ExtraSpecial_B = (6.0*S)
	Min_reads = 20

	LineID= []
	Locations=[]


#	User defined
	Var_neutral=0.0860743795578-0.0659822278933

	df=[]
	bs=[]
	percentiles=[[] for i in range(100)]
	for line_idx, line in enumerate(in1):
		cols = line.replace('\n', '').split('\t')		
		df.append(float(cols[0]))
		# bs.append(float(cols[1]))
		for j in range(9):
			percentiles[line_idx].append(float(cols[j+1]))
		rx=(percentiles[line_idx][2]+percentiles[line_idx][0]-2*percentiles[line_idx][1])/(percentiles[line_idx][2]-percentiles[line_idx][0])
		bs.append(rx)
				
#       Counters
        depth_cc=0
	raw_read_count=0
	num_snps=0
	Var_snp_specific=0.0
	zraw=[]
	z_std=[]
	q_C=[]
	q_T=[]

        # evaluate contents of each line of input file
        for line_idx, line in enumerate(src):
                cols = line.replace('\n', '').split('\t') 

                if len(cols) < 2:               ## This should be header
                        continue
                elif cols[0] == "#CHROM":
			for i in range(len(cols)):
				print i,cols[i]
				if i>8: 
					LineID.append(cols[i])
							

                else:
			scaff=cols[0].split('_') 
			position = int(cols[1])
		        ref_base = cols[3]
		        alt_base = cols[4]
		        if len(alt_base) > 1:	# multiple bases at site
				pass #
		        else:
				C_count=0
				C_dat=[]
				T_count=0
				T_dat=[]

#CA.rmdup.bam	CB.rmdup.bam	S.rmdup.bam	TA.rmdup.bam	TB.rmdup.bam

#changes here specific to sample ordering
				for j in range(9,9+5):
					if len(cols)<(9+5):
						print "whoa ",line_idx,len(cols)
					else:
						info=cols[j].split(':') 
					if len(info)==5 and (j==12):	# data present for TA
						AD=info[1].split(",")
						if int(AD[0])+int(AD[1]) > 0:
							C_count+=1
							C_dat.append(int(AD[0]))
							C_dat.append(int(AD[1]))

						 
					if len(info)==5 and (j==13):	# data present						
						AD=info[1].split(",")
						if int(AD[0])+int(AD[1]) > 0:
							T_count+=1
							T_dat.append(int(AD[0]))
							T_dat.append(int(AD[1]))
						

				# print dip_count,tet_count
				if C_count >= 0 and T_count >= 0:
					#print "here"
					qC=[0.0,0.0]
					qT=[0.0,0.0]
										
					for j in range(len(C_dat)/2):
						#print dip_dat[2*j],dip_dat[2*j+1]
						m = C_dat[2*j]+C_dat[2*j+1]
						qC[0]+= float(C_dat[2*j])	# ref allele freq						
						qC[1]+= float(m)
					for j in range(len(T_dat)/2):
						m = T_dat[2*j]+T_dat[2*j+1]
						qT[0]+= float(T_dat[2*j])	# ref allele freq						
						qT[1]+= float(m)

					#print qT[1],qC[1]						
					if (qT[1] >= Min_reads and qC[1] >= Min_reads):
						Locations.append(cols[0]+"_"+cols[1])
						num_snps+=1
						qC_hat=qC[0]/qC[1]
						qT_hat=qT[0]/qT[1]

						var_C = 1.0/qC[1]
						var_T = 1.0/qT[1]

						Var_snp_specific+=(var_C + var_T)
						vdiv = Var_neutral + var_C + var_T
						diverge=2.0*(math.asin(qT_hat**0.5)-math.asin(qC_hat**0.5))
						out0.write(cols[0]+'\t'+cols[1]+'\t'+str(qC[1])+'\t'+str(qC_hat)+'\t'+str(qT[1])+'\t'+str(qT_hat)+'\t'+str(diverge)+'\t'+str(vdiv)+'\t'+str(diverge/(vdiv**0.5))+'\n')

						if random.randrange(1,3) == 1:
							zraw.append(diverge)
						else:
							zraw.append(-diverge)
						z_std.append(diverge/(vdiv**0.5))
						q_T.append(qT_hat)
						q_C.append(qC_hat)
						out2.write(str(diverge/(vdiv**0.5))+'\n')


                if line_idx % 100000 == 0:
                    print cols[0],line_idx
	print "Included:Number of snps ",num_snps
	print "Sampling/genotyping variance ", Var_snp_specific/num_snps        
	ranked_z=sorted(zraw)
	n25=ranked_z[int(num_snps/4)]
	n50=ranked_z[int(num_snps/2)]
	n75=ranked_z[int(3*num_snps/4)]
	print "Z percentiles ",n25,n50,n75
	print "total variance in Z (based on IQR) ",((n75-n25)/1.349)**2
	out1.write(cols[0]+'\t'+str(Var_snp_specific/num_snps)+'\t'+str(((n75-n25)/1.349)**2)+'\n')					

	ranked_z=sorted(z_std)
	n25=ranked_z[int(num_snps/4)]
	n50=ranked_z[int(num_snps/2)]
	n75=ranked_z[int(3*num_snps/4)]
	print "Zs percentiles ",n25,n50,n75
							
	Braw=[]
	Bloc=[]
	for k in range(S,num_snps):
		if k % S == 0 or k % S == S/2:
			b=0.0
			for j in range(S):
				b+=( z_std[k-j]**2 )
			Bloc.append(Locations[k])
			Braw.append(b)
			if b > ExtraSpecial_B:
				for j in range(S):
					out4.write(Locations[k]+'\t'+str(b)+'\t'+str(z_std[k-S+j+1])+'\t'+str(q_T[k-S+j+1])+'\t'+str(q_C[k-S+j+1])+'\n')

	ranked_B=sorted(Braw)
	n25=ranked_B[int(len(Braw)/4)]
	n50=ranked_B[int(len(Braw)/2)]
	n75=ranked_B[int(3*len(Braw)/4)]
	print "B percentiles ",n25,n50,n75
	b_skew=(n75+n25-2*n50)/(n75-n25)
	print "B Bowley skew ",b_skew

	m = -1
	if b_skew > bs[0]:
		print "Too much skew"
	else:
		for j in range(1,len(bs)):
			if b_skew > bs[j]:
				m = df[j]
				jstar=j
				break

	print "Degrees of freedom ",m
	cIQR=percentiles[jstar][2]-percentiles[jstar][0]
	sigB=(n75-n25)*(2*m)**0.5/cIQR

	for j in range(len(ranked_B)):
		bx=Braw[j]
		bs=m+(bx-S)*((2*m)**0.5)/sigB
		
		if bs < percentiles[jstar][3]: # p greater than 0.05
			p=0.5
		elif bs > percentiles[jstar][8]:
			p = 5.0* 10**(1.0-8.5) # p less than table min
		else:
			for k in range(4,9):
				if bs > percentiles[jstar][k-1] and bs <= percentiles[jstar][k]:
					dx = (bs-percentiles[jstar][k-1])/(percentiles[jstar][k]-percentiles[jstar][k-1])
					x = k-1+dx
					p = 5.0* 10**(1.0-x)
					break

		if j > 0:
			midpoint=str(Bloc[j-1])
		else:
			midpoint=str(Bloc[j])
		out3.write(midpoint+'\t'+str(bx)+'\t'+str(bs)+'\t'+str(p)+'\n')




