#! /usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="A program to graph the quality score distribution by nucleotide")
    parser.add_argument('-f', type = str, help = 'filename', required = True)
    parser.add_argument('-l', type = int, help = 'read length', required = True)
    parser.add_argument('-o', type = str, help = 'output plot name', required = True)
    return parser.parse_args()
args = get_args()

fq_file = args.f
read_length = args.l
fout = args.o

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return(ord(letter)-33)

line_count = 0
sum_lst: list = []
sum_lst = ([0.0] * read_length)
with gzip.open(fq_file, "rt") as fh:   
        for line in fh:
            line_count += 1
            line = line.strip('\n')
            if line_count%4 == 0:
                for pos, qual_score in enumerate(line):
                    sum_lst[pos] += convert_phred(qual_score)
        for pos, sum_scores in enumerate(sum_lst):
            mean = float(sum_scores/(line_count/4))
            sum_lst[pos] = mean

fig, ax = plt.subplots()
for i, mean in enumerate(sum_lst):
    ax.bar(i, mean, facecolor='#580F41')
    ax.set_xlabel('Position')
    ax.set_ylabel('Mean Quality Score')
    ax.set_title('Mean Quality Score by Position in Read')
plt.savefig(f"{fout}_dist_plot.png", ) 
