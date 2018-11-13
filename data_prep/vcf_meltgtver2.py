#!/usr/bin/env python
""" Melt a VCF file into a tab delimited set of calls, one per line
VCF files have all the calls from different samples on one line.  This
script reads vcf on stdin and writes all calls to stdout in tab delimited
format with one call in one sample per line.  This makes it easy to find
a given sample's genotype with, say, grep.
"""

import sys
import csv
import vcf

#out = csv.writer(sys.stdout, delimiter='\t')

if len(sys.argv) > 1:
    inp = str(sys.argv[1])
    heto = str(sys.argv[2])
    GTdat = str(sys.argv[3])
else:
    inp = "No file"
#reader = vcf.VCFReader(inp,'rb')

vcf_reader = vcf.Reader(filename=inp)
#'/storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr2.allfilters.vcf.gz')

reader = vcf_reader

formats = reader.formats.keys()
infos = reader.infos.keys()

OutFile1 = open(GTdat,'wt')
OutFile2 = open(heto,'wt')

Out1 = csv.writer(OutFile1,delimiter=',')
Out2 = csv.writer(OutFile2,delimiter=',')


def flatten(x):
    if type(x) == type([]):
        x = ','.join(map(str, x))
    return x

for record in reader:

    for sample in record.samples:
        if str(sample.gt_alleles[0]) != "None":
           if ":" not in str(record.ID):
               if int(sample.gt_alleles[0]) > 0 or int(sample.gt_alleles[1]) > 0:
                   row1 = [sample.sample]
                   row1 += [record.ID]
                   row2 = row1
                   row1 += [sample.gt_alleles[0]]
                   row2 += [sample.gt_alleles[1]]
               
                   Out1.writerow(row1)
               if int(sample.gt_alleles[0]) == 0 and int(sample.gt_alleles[1]) == 0:
                   row = [sample.sample]
                   row += [record.ID]
                   Out2.writerow(row)

OutFile1.close()
OutFile2.close()

