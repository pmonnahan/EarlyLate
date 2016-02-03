# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 10:29:31 2016

@author: patrick
"""

Key="/Users/patrick/Documents/Research/Mimulus/Collichio_IDLinkKeyfor_MimNewAnnotated.csv"
GeneListDir="/Users/patrick/Documents/Research/Epistasis DEseq/"
GeneListFile="EpiGeneNames.csv"

genes=open(GeneListDir+GeneListFile,"rb")

newFile=open(GeneListDir+"EpiGeneNames_AnnotForm.csv","wb")
newWriter=csv.writer(newFile,dialect="excel",delimiter=",")
with open(Key,"rb") as key: 
    for i,entry in enumerate(key):
        if i==2:
            print entry
        