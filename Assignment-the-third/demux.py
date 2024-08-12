#!/usr/bin/env python

import argparse
import gzip
import itertools
import bioinfo


def get_args():
    parser = argparse.ArgumentParser(description="Demultiplexing of dual matched indexes")
    parser.add_argument('-i', type = str, help = 'indexes filename', default = '/projects/bgmp/shared/2017_sequencing/indexes.txt')
    parser.add_argument('-f1', type = str, help = 'file 1 input', required = True)
    parser.add_argument('-f2', type = str, help = 'file 2 input', required = True)
    parser.add_argument('-f3', type = str, help = 'file 3 input', required = True)
    parser.add_argument('-f4', type = str, help = 'file 4 input', required = True)
    parser.add_argument('-o', type = str, help = 'summary info filename', default = '/projects/bgmp/carlyham/bioinfo/Bi621/Demultiplex/Assignment-the-third/output_files/01summary_stats.txt')
    return parser.parse_args()
args = get_args()

f1 = args.f1
f2 = args.f2
f3 = args.f3
f4 = args.f4


def reverse_complement(seq: str) -> str:
    """Given a DNA sequence, will find the reverse complement of that
    sequence."""
    rev_seq = ''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    for base in reversed(seq.upper()):
        rev_seq = rev_seq + complement[base]
    return(rev_seq)


#parse indexes file and save only the index sequence to the indexes list
indexes = []
matched_counts = {}
#create a dictionary with indexes as keys, 0s as values to count matched indexes later
lines = 0
with open(args.i, "r") as fh:
    for line in fh:
        lines += 1
        line = line.strip("\n")
        line = line.split()
        if line[4] not in indexes and lines != 1: 
            indexes.append(line[4])
            matched_counts[line[4]] = 0
            #append index sequence line to indexes, create key-value pair in dictionary

matched_files = {}
#assign open files for each index to dictionary keys for later in demultiplexing when files can't be opened within loop
for index in indexes:
    matched_files[index] = (open(f"{index}_R1.fq", "w"), open(f"{index}_R2.fq", "w"))


index_hopped_count = 0
unknown_count = 0
total_reads = 0
#initialize hopped, unknown, and total read counts
seq_indexes = []
qual_score = []
#create lists for indexes from reads and quality scores from reads
with open('unknown_R1.fq', 'w') as unk1, open('unknown_R2.fq', 'w') as unk2, \
    open('index_hopped_R1.fq', 'w') as hop1, open('index_hopped_R2.fq', 'w') as hop2:
    #open unknown and index hopped output files for writing
    with gzip.open(f1, 'rt') as R1, gzip.open(f2, 'rt') as R2, gzip.open(f3, 'rt') as R3, gzip.open(f4, 'rt') as R4:
        #open sequence files for reading
        while True:
        #while reading this file...    
            R1_lines = []
            for line in itertools.islice(R1, 4):
                R1_lines.append(line)

            R2_lines = []
            for line in itertools.islice(R2, 4):
                R2_lines.append(line)

            R3_lines = []
            for line in itertools.islice(R3, 4):
                R3_lines.append(line)

            R4_lines = []
            for line in itertools.islice(R4, 4):
                R4_lines.append(line)
            #use islice itertools function to read 4 lines from each file without saving whole file to memory. 
            #create a list, each item is a line in the 4 lines
            
            if not R1_lines or not R2_lines or not R3_lines or not R4_lines:
                break
            #if not reading 4 file lines, break the loop.
            
            total_reads += 1
            #iterate total reads count by one

            seq_indexes = []
            qual_score = []
            #set/reset index and quality score lists to empty

            seq_indexes.append(R2_lines[1].strip())
            #append index1 to seq_index list
            qual_score.append(R2_lines[3].strip())
            #append index1 quality score line to qual_score list

            seq_indexes.append(reverse_complement(R3_lines[1].strip()))
            #reverse complement and append index2 to list
            qual_score.append(R3_lines[3].strip())
            #append quality score line to list

            R1_lines[0] = (f"{R1_lines[0].strip()} {seq_indexes[0]}-{seq_indexes[1]}\n")
            R4_lines[0] = (f"{R4_lines[0].strip()} {seq_indexes[0]}-{seq_indexes[1]}\n")
            #add both indexes to header line

            if seq_indexes[0] not in indexes or seq_indexes[1] not in indexes:
                for line in R1_lines:
                    unk1.write(f"{line}")
                for line2 in R4_lines:
                    unk2.write(f"{line2}")
                unknown_count += 1
                #if either seq_index is not in the known index list, write the full read to the unknown file
                
            elif seq_indexes[0] in indexes and seq_indexes[1] in indexes and seq_indexes[0] != seq_indexes[1]:
                for line in R1_lines:
                    hop1.write(f"{line}")
                for line2 in R4_lines:
                    hop2.write(f"{line2}")
                index_hopped_count += 1
                #if the indexes are known, but don't match add them to hopped indexes for either read
            
            elif seq_indexes[0] in indexes and seq_indexes[1] in indexes and seq_indexes[0] == seq_indexes[1]:
                #check that both indexes match and are in the indexes list
                lowest_score = 100
                #initialize a variable lowest quality score, resets every loop
                for score in qual_score:
                    for pos in score:
                        if bioinfo.convert_phred(pos) < lowest_score:
                            lowest_score = bioinfo.convert_phred(pos)
                        #conver phred by position, if score is lower than current lowest_score, replace lowest_score
                    
                if lowest_score < 30:
                    for line in R1_lines:
                        unk1.write(f"{line}")
                    for line2 in R4_lines:
                        unk2.write(f"{line2}")
                    unknown_count += 1
                    #if quality score is below 30, write reads to unknown files and iterate unknown count
                            
                else:
                    #if lowest quality score is >= 30...
                    files = matched_files.get(str(seq_indexes[0]))
                    #find pre-opened files for the appropriate index
                    idx_count = matched_counts.get(seq_indexes[0], 0)  # Default to 0 if not found
                    idx_count += 1
                    #iterate assoctated count by 1
                    matched_counts[seq_indexes[0]] = idx_count
                    #add new count back to dictionary
                    (match1, match2) = files
                    for line in R1_lines:
                        match1.write(f"{line}")
                    for line2 in R4_lines:
                        match2.write(f"{line2}")
                    #write each read to appropriate output file
                
percent_matched = []
total_percent = 0
for index, index_count in matched_counts.items():
    percent = (index_count/total_reads)*100
    #calculate percentage matched reads out of total reads
    total_percent += percent
    #add matched index percentage to running sum
    percent_matched.append((index, index_count, percent))
    #add index, count of reads, and percent of reads to percent matched list
print(percent_matched)

idx_hopped_percent = (index_hopped_count/total_reads)*100
unknown_percent = (unknown_count/total_reads)*100
#calculate hopped and unknown read percentage
total_percent += idx_hopped_percent
total_percent +=unknown_percent
#add unknown and hopped read percentage to running total

with open(args.o, "w") as stats:
    stats.write(f"Index\tCount\tPercent of Reads\n")
    for position, index in enumerate(percent_matched):
        stats.write(f"{index[0]}\t{index[1]}\t{index[2]}\n")
    stats.write(f"index hopped\t{index_hopped_count}\t{idx_hopped_percent}\n")
    stats.write(f"unknown\t{unknown_count}\t{unknown_percent}\n")
    stats.write(f"total\t{total_reads}\t{total_percent}\n")
