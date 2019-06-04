#!/usr/bin/env python2


import sys
import re
import os
import ConfigParser
from copy import deepcopy

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../Bin/config.ini')
Bedtools = config.get('SOFTWARE', 'bedtools')

class Bed(object):
    def __init__(self, line_in):
        self.line = line_in
        list_split = re.split('\t', line_in)
        self.chr = list_split[0]
        self.start = int(list_split[1])
        self.end = int(list_split[2])
        self.count = int(list_split[3])
        self.range = self.chr+'-'+str(self.start/100000)

    def __eq__(self, other):
        return self.line == other.line

def CompareFileLen(file1, file2):
    a = 0
    b = 0
    with open(file1, 'r') as f:
        for line in f:
            a += 1
    with open(file2, 'r') as f:
        for line in f:
            b += 1
    if a == b:
        return 1
    else:
        return 0

def GetNewDepth(old, merge, outpath):
    out = open(os.path.join(outpath, 'tmp.txt'), 'w')
    dict_old = {}
    with open(old, 'r') as f:
        for line in f:
            ob = Bed(line.strip())
            if ob.range not in dict_old:
                dict_old[ob.range] = [ob,]
            else:
                dict_old[ob.range] += [ob,]
    for key in dict_old:
        dict_old[key] = set(dict_old[key])
    with open(merge, 'r') as f:
        for line in f:
            list_split = re.split('\t', line.strip())
            list_tmp = []
            start = int(list_split[1])/100000
            end = int(list_split[2])/100000
            for i in range(end-start+1):
                for ob in dict_old[list_split[0]+'-'+str(start+i)]:
                    if ob.start <= int(list_split[1]) and ob.end <= int(list_split[2]):
                        list_tmp.append(ob)
            ob_max = sorted(list_tmp, key=lambda x:x.count, reverse=True)[0]
            out.write('\t'.join([ob_max.chr, str(ob_max.start), str(ob_max.end), list_split[3]])+'\n')
    out.close()
    os.system('mv {} {}'.format(os.path.join(outpath, 'tmp.txt'), old))

def All(depth, outpath, dist):
    os.system('{} merge -i {} -d {} -c 4 -o sum > {}'.
             format(Bedtools, depth, dist, os.path.join(outpath, 'merge.txt')))
    GetNewDepth(depth, os.path.join(outpath, 'merge.txt'), outpath)
    os.system('{} merge -i {} -d {} -c 4 -o sum > {}'.
              format(Bedtools, depth, dist, os.path.join(outpath, 'merge2.txt')))
    if CompareFileLen(os.path.join(outpath, 'merge.txt'),\
                      os.path.join(outpath, 'merge2.txt')):
        GetNewDepth(depth, os.path.join(outpath, 'merge.txt'), outpath)
        os.system('mv {} {}/final.merge.point.depth'.format(depth, outpath))
        os.system('mv {}/merge.txt {}/final.merge.depth'.format(outpath, outpath))
        os.system('rm {}/merge2.txt'.format(outpath))
    else:
        All(depth, outpath, dist)

def main():
    if os.path.exists(sys.argv[2]):
        os.system('rm -rf {}'.format(sys.argv[2]))
    os.mkdir(sys.argv[2])
    base = os.path.basename(sys.argv[1])
    os.system('cp {} {}'.format(sys.argv[1], os.path.join(sys.argv[2], base)))
    All(os.path.join(sys.argv[2], base), sys.argv[2], sys.argv[3])


if __name__ == '__main__':
    main()

