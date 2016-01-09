# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 10:29:31 2016

@author: patrick
"""
import csv

Key="/Users/patrick/Documents/Research/Mimulus/Collichio_IDLinkKeyfor_MimNewAnnotated.csv"
GeneListDir="/Users/patrick/Documents/Research/Epistasis DEseq/"
GeneListFile="QTL2GeneNames.csv"

newFile=open(GeneListDir+"QTL2GeneNames_AnnotForm.txt","wb")

hits=0
nonhits=0
numgenes=0


with open(GeneListDir+GeneListFile,"rU") as genes:
    for j,gene in enumerate(genes):                
        if j != 0:
            hit=False
            numgenes+=1                    
            gene=gene.strip("\n")
            gene=gene.split(".")
            gene=gene[1].strip("\"")
            with open(Key,"rU") as key: 
                for i,entry in enumerate(key):
                    entry=entry.strip("\n")
                    entry=entry.split(",")
                    keyGutGene=entry[4]
                    keyAnnotGene=entry[0]
                    if gene == keyGutGene:
                        xx=str(keyAnnotGene)+"\n"
                        newFile.write(xx)
                        hit=True
            if hit==True:
                hits+=1
            else:   
                print gene,keyGutGene
                nonhits+=1
newFile.close()
print "num hits = ", hits
print "num Nonhits = ",nonhits
print "num queried = ",numgenes


                        