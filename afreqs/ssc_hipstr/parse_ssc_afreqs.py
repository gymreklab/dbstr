#!/usr/bin/env python3
"""
Usage: ./parse_ssc_afreqs.py <VCFFILE>
"""

import vcf
import sys

try:
    VCFFILE = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

reader = vcf.Reader(open(VCFFILE, "rb"))

def GetAlleleCounts(record):
    acounts = {}
    for sample in record:
        gb = sample["GB"]
        if gb is not None:
            alleles = gb.split("|")
            gtalleles = sample["GT"].split("|")
            for item in gtalleles:
                astring = ([record.REF]+record.ALT)[int(item)]
                agb = len(astring)-len(record.REF)
                key = "%s_%s"%(astring, agb)
                acounts[key] = acounts.get(key, 0) + 1
    return acounts

sys.stdout.write(",".join(["STRID","allele_string", "allele_gb", "count"])+"\n")
for record in reader:
    strid = record.ID
    acounts = GetAlleleCounts(record)
    for a in acounts:
        allele = str(a)
        astring = allele.split("_")[0]
        agb = allele.split("_")[1]
        count = acounts[a]
        sys.stdout.write(",".join([strid, astring, agb, str(count)])+"\n")

    
