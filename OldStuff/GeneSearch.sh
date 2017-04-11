#!/bin/sh

# enable echoing of commands
set -x


python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/FT_gene_query_PerPop.py .DS_Store || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/gff_query_PerPop.py .DS_Store || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/MimNewAnnotated_query_perpop.py .DS_Store || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/FT_gene_query_PerPop.py BR13.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/gff_query_PerPop.py BR13.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/MimNewAnnotated_query_perpop.py BR13.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/FT_gene_query_PerPop.py IM13.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/gff_query_PerPop.py IM13.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/MimNewAnnotated_query_perpop.py IM13.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/FT_gene_query_PerPop.py IM14.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/gff_query_PerPop.py IM14.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/MimNewAnnotated_query_perpop.py IM14.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/FT_gene_query_PerPop.py Q13.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/gff_query_PerPop.py Q13.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/MimNewAnnotated_query_perpop.py Q13.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/FT_gene_query_PerPop.py Q14.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/gff_query_PerPop.py Q14.SigIndTests || exit 1 
python /Users/patrick/Documents/Research/EarlyLate/EL_Working_Code/MimNewAnnotated_query_perpop.py Q14.SigIndTests || exit 1 
