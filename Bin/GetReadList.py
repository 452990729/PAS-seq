#!/ysr/bin/env python2

import sys
import re
from glob import glob

fqs = glob(sys.argv[1]+'/*.fastq.gz')

i = 0
for fq in sorted(fqs):
    i += 1
    if i %2 == 0:
        a = a+','+fq
        group, sample = re.findall('(.*)_(.*)_', re.split('/', fq)[-1])[0]
        print '\t'.join([sample, group, a])
    else:
        a = fq


