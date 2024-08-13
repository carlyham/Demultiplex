# Demultiplexing Lab Notebook

## Thursday 07/25/2024

### Initial Data Exploration:

path to files on Talapas: /projects/bgmp/shared/2017_sequencing/

What do we want to know? 

- File size and number of lines per file
- Which files contain barcodes and which contain sequences
- Is R3 the reverse complement of the barcode or the barcode?

ls -lah

![Untitled](Demultiplexing%20Lab%20Notebook%205b25f830949e4e6fb4ed5f33a796495c/Untitled.png)

zcat 1294_S1_L008_R1_001.fastq.gz | wc -l  → 1452986940*

1452986940/4 = 363246735 records per file

*also ran on R3 to confirm num of lines

Head sequence files (R1 and R4):

![Screenshot 2024-07-25 at 11.35.44 AM.png](Demultiplexing%20Lab%20Notebook%205b25f830949e4e6fb4ed5f33a796495c/Screenshot_2024-07-25_at_11.35.44_AM.png)

Head index files (R2 and R3):

![Screenshot 2024-07-25 at 11.39.02 AM.png](Demultiplexing%20Lab%20Notebook%205b25f830949e4e6fb4ed5f33a796495c/Screenshot_2024-07-25_at_11.39.02_AM.png)

Data actually from 08/6 but I wanted it to be with all other exploratory data:

![Screenshot 2024-08-06 at 11.22.06 AM.png](Demultiplexing%20Lab%20Notebook%205b25f830949e4e6fb4ed5f33a796495c/Screenshot_2024-08-06_at_11.22.06_AM.png)

Determine the number of ‘N’ containing barcodes per file:

```bash
(base) [carlyham@login2 Bi622]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | grep 'N' | awk 'NR % 4 == 2' | wc -l 
91805837
(base) [carlyham@login2 Bi622]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | grep 'N' | awk 'NR % 4 == 2' | wc -l 
91643697
```

Apparently supposed to be R2: 3976613 R3: 3328051 total: 7304664 . Determine how to fix this.

## Tuesday 07/30/24

### Creating Quality Score Distributions per Nucleotide

Largely this is the same code as PS4, but modified to run in Talapas as a batch script:

```bash
import argparse
import matplotlib.pyplot as plt
#import bioinfo
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="A program to graph the quality score distribution by nucleotide")
    parser.add_argument('-f', type = str, help = 'filename', required = True)
    parser.add_argument('-l', type = int, help = 'read length', required = True)
    return parser.parse_args()
args = get_args()

fq_file = args.f
read_length = args.l

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
plt.savefig(f"{args.f}_plot.png")
```

My original plan was to import bioinfo and use convert_phred from there as we’ve been doing, but I was unable to get it sbatch to work with my “homemade” module. I used argparse to get the filepath and length of the read to create the “matrix” for the mean. 

## Thursday 08/1/2024

Adding to existing Python 3.12 environment:

```bash
conda activate bgmp_py312
conda install matplotlib
```

Submit sbatch jobs for each read:

```bash
(bgmp_py312) [carlyham@n0349 Demultiplex]$ sbatch --account=bgmp --partition=bgmp --cpus-per-task=4 --mem=16GB ./dist_plot.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" -l 101 -o R1
Submitted batch job 7846252
(bgmp_py312) [carlyham@n0349 Demultiplex]$ sbatch --account=bgmp --partition=bgmp --cpus-per-task=4 --mem=16GB ./dist_plot.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" -l 8 -o R2
Submitted batch job 7846259
(bgmp_py312) [carlyham@n0349 Demultiplex]$ sbatch --account=bgmp --partition=bgmp --cpus-per-task=4 --mem=16GB ./dist_plot.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" -l 8 -o R3
Submitted batch job 7846311
(bgmp_py312) [carlyham@n0349 Demultiplex]$ sbatch --account=bgmp --partition=bgmp --cpus-per-task=4 --mem=16GB ./dist_plot.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" -l 101 -o R4
Submitted batch job 7846313
```

Ahhhhh I forgot to run my scripts with usr bin time -v so I don’t have a good idea of what happened while my plots were generated 😟 but the index plots took less than 30 mins each with these parameters. The reads took much longer than that, over 1:30 each.  Luckily I still feel confident that the plots generated are accurate, I ran my script on the test file from PS4 and got the same results as I did for that assignment. Next time I de this, I would change my graph names to be more informative with each run, I failed to think about needing to add them to [answers.bd](http://answers.bd) in markdown and had to add another header.

## Tuesday 08/06/2024

create function reverse complement to output string as opposed to list or other data type

```bash
def reverse_complement(seq: str) -> str:
    """Given a DNA sequence, will find the reverse complement of that
    sequence."""
    rev_seq = ''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    for base in reversed(seq.upper()):
        rev_seq = rev_seq + complement[base]
    return(rev_seq)
```

Check to make sure function is working:

```bash
#run script with just this function present:
(reverse_complement('ATTNGccGTA')) -> TACGGCNAAT 
(base) (base) [carlyham@n0353 Assignment-the-third]$ ./demux.py  
TACGGCNAAT
type((reverse_complement('ATTNGccGTA'))) -> TACGGCNAAT
(base) (base) [carlyham@n0353 Assignment-the-third]$ ./demux.py 
<class 'str'>
```

Reverse complement seems to be working just fine. Yay! My plan doesn’t require any more functions outside of those in [bioinfo.py](http://bioinfo.py) so I will continue to work on the rest of the script.

## Thursday 08/08/2024

Continued following pseudocode plan for demultiplexing. I am using itertools islice to read the file. My hope with this is that I won't need to store much of anything in memory which will save time and computational power.

## Saturday 08/10/2024

Finished demultiplexing script! When tested on text files everything worked really well which is good. I just need to remember to switch from open to [gzip.open](http://gzip.open) before running on real files. Path to script:

```bash
/projects/bgmp/carlyham/bioinfo/Bi621/Demultiplex/Assignment-the-third/demux.py
```

Path to shell script to run demux.py:

```bash
/projects/bgmp/carlyham/bioinfo/Bi621/Demultiplex/Assignment-the-third/demux_sbatch.sh
```

Usr/bin/time output:

```bash
[('GTAGCGTA', 5774439, 1.5896740269392924), ('CGATCGAT', 4237854, 1.1666599012927124), ('GATCAAGG', 4628196, 1.2741190915315452), ('AACAGCGA', 6368144, 1.753118028714009), ('TAGCCATG', 7148153, 1.9678505850850938), ('CGGTAATC', 2393021, 0.6587866508972201), ('CTCTGGAT', 24515042, 6.748867818454032), ('TACCGGAT', 49686878, 13.678547723216287), ('CTAGCTCA', 13034311, 3.5882802910809373), ('CACTTCAC', 2577666, 0.709618491133857), ('GCTACTCT', 4301318, 1.1841312214409856), ('ACGATCAG', 5933528, 1.6334704288532695), ('TATGGCAC', 7651472, 2.106411775456151), ('TGTTCCGT', 11450554, 3.152279950981528), ('GTCCTAAG', 6200133, 1.7068654450534841), ('TCGACAAG', 2644260, 0.7279514845467228), ('TCTTCGAC', 30089661, 8.283532403945765), ('ATCATGCG', 6927867, 1.907206956725984), ('ATCGTGGT', 4730009, 1.3021476985883988), ('TCGAGAGT', 7448072, 2.050416778006277), ('TCGGATTC', 2874320, 0.7912858459691317), ('GATCTTGC', 2636332, 0.7257689460030522), ('AGAGTCCA', 7602663, 2.09297490313299), ('AGGATAGC', 5861709, 1.6136990192079772)]
	Command being timed: "./demux.py -f1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -f2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -f3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -f4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
	User time (seconds): 3920.17
	System time (seconds): 47.41
	Percent of CPU this job got: 90%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:13:17
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 248996
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 38770
	Voluntary context switches: 48628
	Involuntary context switches: 2260
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

Looks like I didn’t need to allocate as many resources as I did to run this script, it only used a quarter of the allotted CPUs. I also accidentally printed my list of summary stats to std out, whoops!

Summary stats from 01summary_stats.txt:

```bash
Index	Count	Percent of Reads
GTAGCGTA	5774439	1.5896740269392924
CGATCGAT	4237854	1.1666599012927124
GATCAAGG	4628196	1.2741190915315452
AACAGCGA	6368144	1.753118028714009
TAGCCATG	7148153	1.9678505850850938
CGGTAATC	2393021	0.6587866508972201
CTCTGGAT	24515042	6.748867818454032
TACCGGAT	49686878	13.678547723216287
CTAGCTCA	13034311	3.5882802910809373
CACTTCAC	2577666	0.709618491133857
GCTACTCT	4301318	1.1841312214409856
ACGATCAG	5933528	1.6334704288532695
TATGGCAC	7651472	2.106411775456151
TGTTCCGT	11450554	3.152279950981528
GTCCTAAG	6200133	1.7068654450534841
TCGACAAG	2644260	0.7279514845467228
TCTTCGAC	30089661	8.283532403945765
ATCATGCG	6927867	1.907206956725984
ATCGTGGT	4730009	1.3021476985883988
TCGAGAGT	7448072	2.050416778006277
TCGGATTC	2874320	0.7912858459691317
GATCTTGC	2636332	0.7257689460030522
AGAGTCCA	7602663	2.09297490313299
AGGATAGC	5861709	1.6136990192079772
index hopped	707740	0.19483726398807136
unknown	135823393	37.39149726975523
total	363246735	100.0
```

Percentages add up to 100 which is promising, and number of reads the same as calculated during file exploration. File is tab separated, can be parsed on command line if needed.