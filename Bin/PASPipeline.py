#!/usr/bin/env python2

import sys
import re
import os
import argparse
import ConfigParser


BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../Bin/config.ini')

#### SOFT
FASTP = config.get('SOFTWARE', 'fastp')
BOWTIE2 = config.get('SOFTWARE', 'bowtie2')
SAMTOOLS = config.get('SOFTWARE', 'samtools')
BEDTOOLS = config.get('SOFTWARE', 'bedtools')
PYTHON = config.get('SOFTWARE', 'python')

### SCRIPT
FilterPositions = config.get('SCRIPT', 'FilterPosition')
MergeDepths = config.get('SCRIPT', 'MergeDepth')

#### DATABASE
HG19 = config.get('DATABASE', 'hg19')

class ReadList(object):
    def __init__(self, line_in):
        list_split = re.split('\s+', line_in)
        self.Sample = list_split[0]
        self.Group = list_split[1]
        self.Name = '_'.join(list_split[:2])
        list_fq = re.split(',', list_split[2])
        self.fq1 = list_fq[0]
        self.fq2 = list_fq[1]

class Snake(object):
    def __init__(self, process):
        self.process = process
        self.input = ''
        self.output = ''
        self.params = ''
        self.log = ''
        self.threads = ''
        self.shell = ''

    def UpdateInput(self, line_in):
        self.input = line_in

    def UpdateOutput(self, line_in):
        self.output = line_in

    def UpdateParams(self, line_in):
        self.params = line_in

    def UpdateLog(self, line_in):
        self.log = line_in

    def UpdateThreads(self, line_in):
        self.threads = line_in

    def UpdateShell(self, line_in):
        self.shell = line_in

    def WriteStr(self, fn):
        fn.write('rule '+self.process+':\n')
        fn.write('\tinput:\n\t\t'+self.input+'\n')
        if self.output:
            fn.write('\toutput:\n\t\t'+self.output+'\n')
        if self.params:
            fn.write('\tparams:\n\t\t'+self.params+'\n')
        if self.log:
            fn.write('\tlog:\n\t\t'+self.log+'\n')
        if self.threads:
            fn.write('\tthreads: '+self.threads+'\n')
        if self.shell:
            fn.write('\tshell:\n\t\t'+self.shell+'\n')
        fn.write('\n')


def main():
    parser = argparse.ArgumentParser(description="PAS pipeline")
    parser.add_argument('-c', help='the input fasta list', required=True)
    parser.add_argument('-o', help='the output path', required=True)
    parser.add_argument('-a1', help='the read1 adapter', default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA')
    parser.add_argument('-a2', help='the read2 adapter', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    argv=vars(parser.parse_args())
    outpath = argv['o']
    snakefile = open(os.path.join(outpath, 'snakefile.txt'), 'w')
    RawData = os.path.join(outpath, 'RawData')
    list_ob = []
    if not os.path.exists(RawData):
        os.mkdir(RawData)
    else:
        os.system('rm -rf '+RawData)
        os.mkdir(RawData)
    with open(argv['c'], 'r') as f:
        os.chdir(RawData)
        for line in f:
            if not line.startswith('#'):
                ob = ReadList(line.strip())
                list_ob.append(ob)
                if ob.fq1.endswith('.gz'):
                    lb = '.fq.gz'
                else:
                    lb = '.fastq'
                os.system('ln -s {} {}'.format(ob.fq1, ob.Name+'_1'+lb))
                os.system('ln -s {} {}'.format(ob.fq2, ob.Name+'_2'+lb))
    os.chdir(outpath)
    ### config file
    snakefile.write('Samples = "{}".split()\n'.format(' '.join([i.Name for i in\
                                                                list_ob])))
    snakefile.write('adapter1 = "{}"\nadapter2 = "{}"\n'.format(argv['a1'], argv['a2']))
    snakefile.write('HG19 = "{}"\n'.format(HG19+'/hg19.fa'))
    snakefile.write('HG19Bed = "{}"\n'.format(HG19+'/hg19.genome'))

    ###all
    All = Snake('All')
    All.UpdateInput('expand(dir("4.Coverage/{sample}"), sample=Samples)')
    All.WriteStr(snakefile)

    ###QC
    QC = Snake('QC')
    QC.UpdateInput('"RawData/{sample}_1.fq.gz"')
    QC.UpdateOutput('"1.QC/{sample}_1.clean.fq.gz"')
    QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
    QC.UpdateShell(r'"'+FASTP+r' -i RawData/{wildcards.sample}_1.fq.gz -o 1.QC/{wildcards.sample}_1.clean.fq.gz -I RawData/{wildcards.sample}_2.fq.gz -O 1.QC/{wildcards.sample}_2.clean.fq.gz --adapter_sequence {adapter1} --adapter_sequence_r2 {adapter2} -j 1.QC/{wildcards.sample}_QC_report.json -h 1.QC/{wildcards.sample}_QC_report.html"')
    QC.WriteStr(snakefile)

    ### Align
    Align = Snake('Align')
    Align.UpdateInput('"1.QC/{sample}_1.clean.fq.gz"')
    Align.UpdateOutput('"2.Align/{sample}.bam"')
    Align.UpdateThreads('5')
    Align.UpdateLog('e = "logs/{sample}.align.e", o = "logs/{sample}.align.o"')
    Align.UpdateShell(r'"'+BOWTIE2+r' --sensitive-local -p 5 -x {HG19} -1 1.QC/{wildcards.sample}_1.clean.fq.gz -2 1.QC/{wildcards.sample}_2.clean.fq.gz |'+SAMTOOLS+r' view -Sb -@ 4 -T {HG19} -o {output}"')
    Align.WriteStr(snakefile)

    ### bamtobed
    Bamtobed = Snake('Bamtobed')
    Bamtobed.UpdateInput('"2.Align/{sample}.bam"')
    Bamtobed.UpdateOutput('"3.Bed/{sample}.bed"')
    Bamtobed.UpdateLog('e = "logs/{sample}.bamtobed.e", o = "logs/{sample}.bamtobed.o"')
    Bamtobed.UpdateShell(r'"'+BEDTOOLS+r' bamtobed -i 2.Align/{wildcards.sample}.bam -cigar > {output}"')
    Bamtobed.WriteStr(snakefile)

    ### FilterPosition
    FilterPosition = Snake('FilterPosition')
    FilterPosition.UpdateInput('"3.Bed/{sample}.bed"')
    FilterPosition.UpdateOutput('"3.Bed/{sample}.bed.filter"')
    FilterPosition.UpdateLog('e = "logs/{sample}.bedfilter.e", o = "logs/{sample}.bedfilter.o"')
    FilterPosition.UpdateShell(r'"'+PYTHON+r' '+FilterPositions+r' 3.Bed/{wildcards.sample}.bed 2 > {output}"')
    FilterPosition.WriteStr(snakefile)

    ### Sort
    Sort = Snake('Sort')
    Sort.UpdateInput('"3.Bed/{sample}.bed.filter"')
    Sort.UpdateOutput('"3.Bed/{sample}.bed.filter.sort"')
    Sort.UpdateLog('e = "logs/{sample}.sort.e", o = "logs/{sample}.sort.o"')
    Sort.UpdateShell(r'"'+BEDTOOLS+r' sort -i 3.Bed/{wildcards.sample}.bed.filter > {output}"')
    Sort.WriteStr(snakefile)

    ### Genomecov
    Genomecov = Snake('Genomecov')
    Genomecov.UpdateInput('"3.Bed/{sample}.bed.filter.sort"')
    Genomecov.UpdateOutput('"4.Coverage/{sample}.depth"')
    Genomecov.UpdateLog('e = "logs/{sample}.depth.e", o = "logs/{sample}.depth.o"')
    Genomecov.UpdateShell(r'"'+BEDTOOLS+r' genomecov -bg -i 3.Bed/{wildcards.sample}.bed.filter.sort -g {HG19Bed} > {output}"')
    Genomecov.WriteStr(snakefile)

    ### MergeDepth
    MergeDepth = Snake('MergeDepth')
    MergeDepth.UpdateInput('expand("4.Coverage/{sample}.depth", sample=Samples)')
    MergeDepth.UpdateOutput('expand(dir("4.Coverage/{sample}/"), sample=Samples)')
#    MergeDepth.UpdateLog('e = "logs/{sample}.merge.e", o = "logs/{sample}.merge.o"')
    MergeDepth.UpdateShell(r'"'+PYTHON+' '+MergeDepths+r' 4.Coverage/{wildcards.sample}.depth {output} 24"')
    MergeDepth.WriteStr(snakefile)


if __name__ == '__main__':
    main()






