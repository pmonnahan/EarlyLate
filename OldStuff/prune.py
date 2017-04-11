# -*- coding: utf-8 -*-
"""
Created on Tue May 26 20:06:32 2015

@author: patrick
"""
new=open("/Users/patrick/Google Drive/Research/EarlyLate/counts_5-26-15_pruned.txt","wb")
with open("/Users/patrick/Google Drive/Research/EarlyLate/counts_bothYears_5-26-15.txt","rb") as sites_file:
    dist=0.0    
    for i,site in enumerate(sites_file):
        if i<=1:
            new.write(site)
        line=site
        site=site.strip("\n")
        site=site.split("\t")
        scaff=site[0]
        pos=site[1]
        if i>1:
            if scaff==old_scaff:
                if (int(pos)-int(old_pos))>50:
                    new.write(line)
                    dist+=int(pos)-int(old_pos)
        old_scaff=scaff
        old_pos=pos
    dist=dist/i
    print dist
new.close()