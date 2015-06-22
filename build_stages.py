#!/usr/bin/env python

import sys
import os
import argparse
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Build ChIP-seq analysis stages")
    parser.add_argument('-stage', dest='stage')
    parser.add_argument('-species', dest='species')
    parser.add_argument('-gene', dest='gene')
    parser.add_argument('-condition', dest='condition')
    parser.add_argument('-bowtie_library', dest='bowtie_library')
    parser.add_argument('-genomic_feature', dest='genomic_feature')
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)
    path = parsed.species + '/' + parsed.gene + '/' + parsed.condition + '/'

    # build stage 1
    if parsed.stage == '1':
        filename_sample_index_key = path + 'samples_index_key.txt'
        if not os.path.isfile(filename_sample_index_key):
            sys.exit("No sample index key found.")
        else:
            index_key = numpy.loadtxt(filename_sample_index_key, dtype=str, skiprows=1)
            if not os.path.isdir(path + 'analysis/'):
                os.mkdir(path + 'analysis/')
            writer = open(path + 'analysis/rBowtie', 'w')
            for i in range(len(index_key)):
                writer.write('bowtie %s ' % (parsed.species + '/' + parsed.bowtie_library))
                writer.write('-q %s -m 1 ' % (path + 'samples/' + index_key[i,0]))
                writer.write('> %s\n' % (path + 'analysis/' + index_key[i,1] + '.1.hits'))
            writer.close()

    # build stage 2
    elif parsed.stage == '2':
        filename_sample_index_key = path + 'samples_index_key.txt'
        if not os.path.isfile(filename_sample_index_key):
            sys.exit("No sample index key found.")
        else:
            index_key = numpy.loadtxt(filename_sample_index_key, dtype=str, skiprows=1)
            if len(index_key) == 6:
                gene_tag = [x for x in index_key[:,1] if (x.endswith('IP') and x!='WT.IP')][0].split('.')[0]
                pairs = []
                pairs.append([gene_tag + '.IP', gene_tag + '.IN'])
                pairs.append([gene_tag + '.IP', gene_tag + '.MOCK'])
                pairs.append([gene_tag + '.IP', 'WT.IP'])
                pairs.append([gene_tag + '.IP', 'WT.IN'])
                pairs.append([gene_tag + '.IP', 'WT.MOCK'])
                pairs.append([gene_tag + '.IN', gene_tag + '.MOCK'])
                pairs.append(['WT.IP', 'WT.IN'])
                pairs.append(['WT.IP', 'WT.MOCK'])
                for i in range(len(pairs)):
                    peaks_dir = path + 'analysis/peaks_' + pairs[i][0] + '_' + pairs[i][1] + '/'
                    if not os.path.isdir(peaks_dir):
                        os.mkdir(peaks_dir)
                    writer = open(peaks_dir + 'rMACS', 'w')
                    writer.write('macs14 --gsize=19000000 ')
                    writer.write('-t %s ' % (path + 'analysis/' + pairs[i][0] + '.1.hits'))
                    writer.write('-c %s ' % (path + 'analysis/' + pairs[i][1] + '.1.hits'))
                    writer.write('-n %s\n' % (peaks_dir + 'MACS'))
                    writer.close()
            else:
                sys.exit("Sample count should be 6, i.e. [gene-tag].IN, [gene-tag].IP, [gene-tag].MOCK, WT.IN, WT.IP, WT.MOCK")

    # build stage 3
    elif parsed.stage == '3':
        peaks_dirs = [path + 'analysis/' + x for x in os.listdir(path + 'analysis/') if (os.path.isdir(path + 'analysis/' + x))]
        for peaks_dir in peaks_dirs:
            writer = open(peaks_dir + '/rAnnotatePeaks', 'w')
            writer.write('./annotatePeaks.rb ')
            writer.write('-g %s ' % (parsed.species + '/' + parsed.genomic_feature))
            writer.write('-b %s ' % (peaks_dir + '/MACS_summits.bed'))
            writer.write('> %s\n' % (peaks_dir + '/MACS_summits.bed.neighbors'))
            writer.write('grep CNAG %s > %s\n' % (peaks_dir + '/MACS_summits.bed.neighbors', peaks_dir + '/MACS_summits.bed.neighbors_genes'))
            writer.close()

if __name__ == "__main__":
    main(sys.argv)
