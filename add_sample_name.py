#!/usr/bin/env python
import sys
import os
import re
import argparse

def run(args):
    sample_file_name = args['sample_file_name'][0]
    match = re.match(args["regexp"], os.path.basename(sample_file_name))
    if match is None:
        raise Exception("unable to parse sample id from %s"%os.path.basename(sample_file_name))
    sample = match.groups()[0]
    sample=re.sub("_","-",sample)

    if args["dry_run"]:
        print "(parsed sample name %s from filename %s)"%(sample, sample_file_name)
        return

 
    seq_number = 1
    for record in sys.stdin:
        match = re.match("^>(\S+)", record )
        if match is not None:
            illumina_moniker = match.groups()[0]
            sys.stdout.write(">%s_%07d %s\n"%(sample, seq_number, illumina_moniker))
            seq_number += 1
        else:
            sys.stdout.write(record)

def get_options():
    description = """
    """

    long_description = """
# fasta text is streamed in , and each sequence name is edited using the
# command line argument, which is a sample filename like
# /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/buildR1/processed_S398Buccal-SD-75MG_S96_L001_R2_001.fastq.trimmed.combined
# or
# /dataset/public_invermay_scratch/scratch/SL_Vac13/joined/VAC13-003_S3_L001_R1_001.fastq.combined.fasta
# example :
# cat /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/buildR1/R1_combined.fa | ./add_sample_name.py  /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/buildR1/processed_S398Buccal-SD-75MG_S96_L001_R2_001.fastq.trimmed.combined
#
# example specifying an alt regexp
# tardis -d /dataset/gseq_processing/scratch/batch_qiime/afm_test/qiime_analysis -c 999999999  cat _condition_fastq2fasta_input_/dataset/Rumen_Livestock/scratch/PRJ0255318_FILES/39_samples_16S_Data/scratch/J293_2.fq.gz \| /dataset/gseq_processing/active/bin/batch_qiime_prism/add_sample_name_new.py -r "\"^([^_]+)_\"" /dataset/Rumen_Livestock/scratch/PRJ0255318_FILES/39_samples_16S_Data/scratch/J293_2.fq.gz \> _condition_uncompressedtext_output_/dataset/gseq_processing/scratch/batch_qiime/afm_test/qiime_analysis/J293_2.fq.gz.combined.fasta
#
# dry run to test parsing
# ./add_sample_name.py -n -r "^([^_]+)_" /dataset/Rumen_Livestock/scratch/PRJ0255318_FILES/39_samples_16S_Data/scratch/J293_2.fq.gz
    """

    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sample_file_name', type=str, nargs=1,help='sample_file_name')  # not opened  - we only use the name
    parser.add_argument('-r', '--regexp', dest='regexp', type=str, metavar='regexp to parse sampleidf', default = "^[^_]+_(\S+)_", help="regexp to parse sampleid")
    parser.add_argument('-n','--dry_run', dest='dry_run', action='store_const', default = False, const=True, help='dry run only - just to test parsing ')
    
    args = vars(parser.parse_args())

    return args


if __name__ == "__main__":
   args=get_options()
   run(args)

