# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here: ```https://github.com/carlyham/Demultiplex/blob/master/dist_plot.py```

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 | +33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.

       Read 1 Plot:
       
       ![](https://github.com/carlyham/Demultiplex/blob/master/Assignment-the-first/plots/R2_dist_plot.png)

       Index 1 Plot:
       
       ![](https://github.com/carlyham/Demultiplex/blob/master/Assignment-the-first/plots/R3_dist_plot.png)

       Index 2 Plot:
       
       ![](https://github.com/carlyham/Demultiplex/blob/master/Assignment-the-first/plots/R1_dist_plot.png)

       Read 2 Plot:
       
       ![](https://github.com/carlyham/Demultiplex/blob/master/Assignment-the-first/plots/R4_dist_plot.png)

    3. Justify qual score cutoff:
       Generally, a quality score threshold of 30 is considered acceptable and represents a 99.9% base call accuracy. Comparing the lowest quality score to this cutoff would generate better data for the index reads. For these reads, we want a strict threshold because we want good filtering out of indexes that could match eachother by mistake (although this is very unlikely).
    For biological reads (R1 and R4) a quality score cutoff is most likely not necessary. The aligner will take care of these reads in a future step. Sequences with mistakes (of low quality) will likely not align to the genome anywhere.

    5. 

        zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | grep 'N' | awk 'NR % 4 == 2' | wc -l 
            91805837
        zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | grep 'N' | awk 'NR % 4 == 2' | wc -l 
            91643697

## Part 2
1. Define the problem
    
   Reads from multiple places were run on the same flow cell and indexes were added to ensure we could determine which reads were which. Now they need to be sorted according to wether or not the inexes match, as well as if they are in the library of indexes we are working with. If both of these criteria are met, a quality score check should be done to ensure that data is up to par.
   
3. Describe output

   Fastq files will be created for each output type for each read (index hopped and unknown indexes). In the case of matching, a file will be created for each set of matching barcodes. Reads that meet the criteria of index matched, hopped, and unknown will be in the appropriate files.
   
5. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
6. Pseudocode
7. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
