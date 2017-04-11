# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 09:56:14 2015

@author: patrick
"""

#!/usr/bin/env python

import sys
import os

commandsFile = open("/Users/patrick/Documents/Research/EarlyLate/RunPoPoolation1.sh", 'w')
mainDir = '/Volumes/TOSHIBA EXT/EarlyLate/EL_remap_Oct_2015/RG/BAMS' 


in0=os.listdir(mainDir)

for line_idx, line in enumerate(in0[1:]):
    cols=line.split(".")
    
    sampleName= cols[0]+'.'+cols[1]

#    fq1=cols[0]+'.'+cols[1]+'.R1.fq'
#    fq2=cols[0]+'.'+cols[1]+'.R2.fq'

    print sampleName

#    shellScriptName = 'remapEL.%s.sh' % (sampleName)
#    shellScript = open(shellScriptName, 'w' )
#
#    commandsFile.write('qsub %s \n' % (shellScriptName))

#    shellScript.write("#PBS -N map_%s\n" % (sampleName))
#    shellScript.write("#PBS -l nodes=1:ppn=1,mem=5000m,walltime=48:00:00\n")
#    shellScript.write("#PBS -q default\n")
#    shellScript.write("#PBS -M jkk@ku.edu\n")
#    shellScript.write("#PBS -m a\n")
#    shellScript.write("#PBS -d %s/\n" % (mainDir))
#    shellScript.write("#PBS -e %s/mx_%s.err\n" % (mainDir,sampleName))
#    shellScript.write("#PBS -o %s/mx_%s.out\n\n" % (mainDir,sampleName))


    commandsFile.write("samtools pileup %s > %s.pileup || exit 1\n" % (line,sampleName))
    commandsFile.write("perl ~/src/popoolation_1.2.2/Variance-sliding.pl --input %s.pileup --output %s.pi --measure pi --window-size 5000 --step-size 5000 --min-count 2 --min-coverage 4 --max-coverage 200 --min-qual 20 --pool-size 100 || exit 1\n" % (sampleName,sampleName))
    commandsFile.write("perl ~/src/popoolation_1.2.2/Variance-sliding.pl --input %s.pileup --output %s.theta --measure theta --window-size 5000 --step-size 5000 --min-count 2 --min-coverage 4 --max-coverage 200 --min-qual 20 --pool-size 100 || exit 1\n" % (sampleName,sampleName))
    commandsFile.write("perl ~/src/popoolation_1.2.2/Variance-sliding.pl --input %s.pileup --output %s.D --measure D --window-size 5000 --step-size 5000 --min-count 2 --min-coverage 4 --max-coverage 200 --min-qual 20 --pool-size 100 || exit 1\n" % (sampleName,sampleName))
    commandsFile.write("perl ~/src/popoolation_1.2.2/Variance-at-position.pl --pool-size 100 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 200 --pileup %s.pileup --gtf Mguttatus_256_v2.0.just_exons.gtf --output %s.genes.pi --measure pi || exit 1\n" % (sampleName,sampleName))
    commandsFile.write("perl ~/src/popoolation_1.2.2/Variance-at-position.pl --pool-size 100 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 200 --pileup %s.pileup --gtf Mguttatus_256_v2.0.just_exons.gtf --output %s.genes.theta --measure theta || exit 1\n" % (sampleName,sampleName))
    commandsFile.write("perl ~/src/popoolation_1.2.2/Variance-at-position.pl --pool-size 100 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 200 --pileup %s.pileup --gtf Mguttatus_256_v2.0.just_exons.gtf --output %s.genes.D --measure D || exit 1\n" % (sampleName,sampleName))
#    
#    
#    commandsFile.write("scythe -a illumina_adapters1.fa -o sc.%s %s || exit 1\n" % (fq2,fq2))
#    commandsFile.write("rm -f %s || exit 1\n" % (fq2))
#
#    shellScript.write("sickle pe -f sc.%s -r sc.%s -l 50 -t sanger -o ss.%s -p ss.%s -s Singles.%s.fq || exit 1\n" % (fq1,fq2,fq1,fq2,sampleName))
#    shellScript.write("rm -f sc.%s || exit 1\n" % (fq1))
#    shellScript.write("rm -f sc.%s || exit 1\n" % (fq2))
#
#    shellScript.write("bwa mem /research/jkkelly/Mguttatus_v2.0_256_hardmasked.fa ss.%s ss.%s | samtools view -bS - | samtools sort - %s.sorted\n" % (fq1,fq2,sampleName))
#    shellScript.write("java -Xmx2048m -jar /scratch/xiaofei/softwares/picard/picard-tools-1.102/MarkDuplicates.jar INPUT=%s.sorted.bam OUTPUT=%s.rmdup.bam METRICS_FILE=%s.metrics.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT || exit 1\n" % (sampleName,sampleName,sampleName))
#    shellScript.write("java -Xmx2048m -jar /scratch/xiaofei/softwares/picard/picard-tools-1.102/AddOrReplaceReadGroups.jar I=%s.rmdup.bam O=RG/%s.rmdup.bam SO=coordinate RGID=SeqRUN# RGLB=%s.rmdup.bam RGPL=illumina RGPU=%s.rmdup.bam RGSM=%s.rmdup.bam VALIDATION_STRINGENCY=LENIENT || exit 1 \n" % (sampleName,sampleName,sampleName,sampleName,sampleName ))
#    shellScript.write("samtools index  RG/%s.rmdup.bam RG/%s.rmdup.bam.bai  || exit 1\n" % (sampleName,sampleName))
#    shellScript.write("gzip ss.%s || exit 1\n" % (fq1))
#    shellScript.write("gzip ss.%s || exit 1\n" % (fq2))
#    shellScript.write("rm -f %s.sorted.bam || exit 1\n" % (sampleName))
#
#    shellScript.write("bwa mem /research/jkkelly/Mguttatus_v2.0_256_hardmasked.fa Singles.%s.fq | samtools view -bS - | samtools sort - %s.singles\n" % (sampleName,sampleName))
#    shellScript.write("java -Xmx2048m -jar /scratch/xiaofei/softwares/picard/picard-tools-1.102/MarkDuplicates.jar INPUT=%s.singles.bam OUTPUT=%s.srmdup.bam METRICS_FILE=%s.smetrics.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT || exit 1\n" % (sampleName,sampleName,sampleName))
#    shellScript.write("java -Xmx2048m -jar /scratch/xiaofei/softwares/picard/picard-tools-1.102/AddOrReplaceReadGroups.jar I=%s.srmdup.bam O=RG/%s.srmdup.bam SO=coordinate RGID=SeqRUN# RGLB=%s.srmdup.bam RGPL=illumina RGPU=%s.srmdup.bam RGSM=%s.srmdup.bam VALIDATION_STRINGENCY=LENIENT || exit 1 \n" % (sampleName,sampleName,sampleName,sampleName,sampleName ))
#    shellScript.write("samtools index  RG/%s.srmdup.bam RG/%s.srmdup.bam.bai  || exit 1\n" % (sampleName,sampleName))
#    shellScript.write("gzip Singles.%s.fq || exit 1\n" % (sampleName))
#    shellScript.write("rm -f %s.singles.bam || exit 1\n" % (sampleName))




