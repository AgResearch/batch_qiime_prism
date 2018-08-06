#!/usr/bin/env python
import sys
import os
import re


# fasta text is streamed in , and each sequence name is edited using the 
# command line argument, which is a sample filename like
# /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/buildR1/processed_S398Buccal-SD-75MG_S96_L001_R2_001.fastq.trimmed.combined
# or 
# /dataset/public_invermay_scratch/scratch/SL_Vac13/joined/VAC13-003_S3_L001_R1_001.fastq.combined.fasta
# example : 
# cat /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/buildR1/R1_combined.fa | ./add_sample_name.py  /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/buildR1/processed_S398Buccal-SD-75MG_S96_L001_R2_001.fastq.trimmed.combined

def main():
    sample_file_name = sys.argv[1]
    #print os.path.basename(sample_file_name) #processed_S398Buccal-SD-75MG_S96_L001_R2_001.fastq.trimmed.combined
    #match = re.match("^processed_(\S+)_L0", os.path.basename(sample_file_name))
    #match = re.match("^[^_]+_(\S+)_L0", os.path.basename(sample_file_name))
    match = re.match("^[^_]+_(\S+)_", os.path.basename(sample_file_name))
    if match is None:
        raise Exception("unable to parse sample id from %s"%os.path.basename(sample_file_name))
    sample = match.groups()[0]
    sample=re.sub("_","-",sample)
    #print match.groups()
 
    seq_number = 1
    for record in sys.stdin:
        match = re.match("^>(\S+)", record )
        if match is not None:
            illumina_moniker = match.groups()[0]
            sys.stdout.write(">%s_%07d %s\n"%(sample, seq_number, illumina_moniker))
            seq_number += 1
        else:
            sys.stdout.write(record)

if __name__ == "__main__":
   main()

