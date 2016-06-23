#!/usr/bin/env python

import sys
import os

## CHANGE PATH
path="/Volumes/avery/Research/EarlyLate/GeneSearch/Data/"  

dirList=os.listdir(path)

print dirList
print 'Length:', len(dirList)

command_file = open(path+sys.argv[1], 'w')

command_file.write('#!/bin/sh'+'\n'+'\n')
command_file.write('# enable echoing of commands'+'\n')
command_file.write('set -x'+'\n'+'\n'+'\n')

try:
    dirList.remove('Results')
    dirList.remove('.DS_Store')
except ValueError:
    pass
    
for fq_file in dirList:
    file_name = fq_file.split('.')	
    input_file = file_name[0]+"."+file_name[1]
    command_file.write('python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/FT_gene_query.py %s.csv || exit 1 \n' % (path+input_file))
    command_file.write('python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/gff_query.py %s.csv || exit 1 \n' % (path+input_file))
    command_file.write('python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/MimNewAnnotated_query.py %s.csv || exit 1 \n' % (path+input_file))

command_file.close()
