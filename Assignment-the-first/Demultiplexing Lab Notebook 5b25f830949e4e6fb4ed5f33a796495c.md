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

Head sequence files:

![Screenshot 2024-07-25 at 11.35.44 AM.png](Demultiplexing%20Lab%20Notebook%205b25f830949e4e6fb4ed5f33a796495c/Screenshot_2024-07-25_at_11.35.44_AM.png)

Head barcode files:

![Screenshot 2024-07-25 at 11.39.02 AM.png](Demultiplexing%20Lab%20Notebook%205b25f830949e4e6fb4ed5f33a796495c/Screenshot_2024-07-25_at_11.39.02_AM.png)

Determine the number of ‘N’ containing barcodes per file:

```bash
(base) [carlyham@login2 Bi622]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | grep 'N' | awk 'NR % 4 == 2' | wc -l 
91805837
(base) [carlyham@login2 Bi622]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | grep 'N' | awk 'NR % 4 == 2' | wc -l 
91643697
```

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