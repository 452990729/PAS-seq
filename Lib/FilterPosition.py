#!/usr/bin/env python2


import sys
import re
import os
import ConfigParser


BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../Bin/config.ini')
hg19 = config.get('DATABASE', 'hg19')+'/hg19.fa'


def MakeDict(file_in):
    dict_tmp = {}
    seq = ''
    with open(file_in, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if len(seq) == 0:
                    pass
                else:
                    dict_tmp[label] = seq
                label = re.split('\s+', line.strip())[0].lstrip('>')
                seq = ''
            else:
                seq += line.strip()
        dict_tmp[label] = seq
    return dict_tmp

def GetGigar(list_split, clip,  dict_hg19):
    string = list_split[-1]
    list_pos = [int(i) for i in re.findall('(\d+)', string)]
    if sum(list_pos) > 99:
        pass
    else:
        if string[-1] == 'M':
            if dict_hg19[list_split[0]][int(list_split[2]):int(list_split[2])+2] != 'AA':
                print '\t'.join([list_split[0], list_split[2],]+list_split[2:])
        elif string[-1] == 'S':
            if list_pos[-1] <= clip:
                pos = int(list_split[2])-list_pos[-1]
                if dict_hg19[list_split[0]][pos:pos+2] != 'AA':
                    print '\t'.join([list_split[0], str(pos), str(pos)]+list_split[3:])

def IdentifyPos(line, dict_hg19, clip):
    list_split = re.split('\t', line.strip())
    if list_split[-2] == '+' and int(list_split[-3]) >= 30:
        GetGigar(list_split, clip, dict_hg19)

def main():
    dict_hg19 = MakeDict(hg19)
    with open(sys.argv[1], 'r') as f:
        for line in f:
            IdentifyPos(line, dict_hg19, sys.argv[2])


if __name__ == '__main__':
    main()
