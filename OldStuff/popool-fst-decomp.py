# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:40:10 2015

@author: patrick
"""
import csv

INPUT_FILE="/Volumes/TOSHIBA EXT/EarlyLate/AllPopsSync1.fet"
OUT_FILE=INPUT_FILE+".parsed.csv"

out1=open(OUT_FILE, "wb")
out1x=csv.writer(out1,dialect="excel",delimiter=",")
out1x.writerow(["scaff","pos","numsnp","fracCov","avgmin","BRE14_BRE13","BRE14_BRL13","IME13_IML13","IME13_IME14","IME13_IML14","IML13_IME14","IML13_IML14","IME14_IML14","QE14_QL14","QE14_QE13","QE14_QL13","QL14_QE13","QL14_QL13","QE13_QL13"])
with open(INPUT_FILE,"rb") as sites_file:  
    for i,site in enumerate(sites_file):
            site=site.strip("\n")
            site=site.split("\t")
            scaff,pos,numsnp,fracCov,avgmin,BRE14_BRE13,BRE14_BRL13,BRE14_IME13,BRE14_IML13,BRE14_IME14,BRE14_IML14,BRE14_QE14,BRE14_QL14,BRE14_QE13,BRE14_QL13,BRE13_BRL13,BRE13_IME13,BRE13_IML13,BRE13_IME14,BRE13_IML14,BRE13_QE14,BRE13_QL14,BRE13_QE13,BRE13_QL13,BRL13_IME13,BRL13_IML13,BRL13_IME14,BRL13_IML14,BRL13_QE14,BRL13_QL14,BRL13_QE13,BRL13_QL13,IME13_IML13,IME13_IME14,IME13_IML14,IME13_QE14,IME13_QL14,IME13_QE13,IME13_QL13,IML13_IME14,IML13_IML14,IML13_QE14,IML13_QL14,IML13_QE13,IML13_QL13,IME14_IML14,IME14_QE14,IME14_QL14,IME14_QE13,IME14_QL13,IML14_QE14,IML14_QL14,IML14_QE13,IML14_QL13,QE14_QL14,QE14_QE13,QE14_QL13,QL14_QE13,QL14_QL13,QE13_QL13=site
            line=[scaff,pos,numsnp,fracCov,avgmin]
            line1=[]
            for fst in [BRE14_BRE13,BRE14_BRL13,IME13_IML13,IME13_IME14,IME13_IML14,IML13_IME14,IML13_IML14,IME14_IML14,QE14_QL14,QE14_QE13,QE14_QL13,QL14_QE13,QL14_QL13,QE13_QL13]:
                line1.append(float(fst.split("=")[1]))
            if (all(float(k) == 0.0 for k in line1)):
                pass
            else:
                line+=line1
                out1x.writerow(line)
            
            