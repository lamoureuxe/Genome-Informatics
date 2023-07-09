Assignment 4: Mappability continues
================
Erik Lamoureux
February 14, 2022

-   [Assignment Overview](#assignment-overview)
-   [Important remarks](#important-remarks)
-   [Server IP: orca1.bcgsc.ca](#server-ip-orca1.bcgsc.ca)
-   [username: elamoureux\_bmeg22](#username-elamoureux_bmeg22)
-   [Shared folder location:](#shared-folder-location)
-   [elamoureux\_bmeg22@orca01:/projects/bmeg/](#elamoureux_bmeg22orca01projectsbmeg)
    -   [0. Getting ready](#getting-ready)
    -   [1. Differences between genome
        builds](#differences-between-genome-builds)
        -   [a. SE alignment against hg38 and
            hg19](#a.-se-alignment-against-hg38-and-hg19)
        -   [b. Making the files
            comparable](#b.-making-the-files-comparable)
        -   [c. Analyzing the
            differences](#c.-analyzing-the-differences)
        -   [d. Reads per chromosome](#d.-reads-per-chromosome)
        -   [d. Reads position in the genome
            builds](#d.-reads-position-in-the-genome-builds)
    -   [2. Ambiguity in reads mapping](#ambiguity-in-reads-mapping)
        -   [a. Redoing the hg38
            alignment](#a.-redoing-the-hg38-alignment)
        -   [b. Analyzing the ambiguity](#b.-analyzing-the-ambiguity)
        -   [c. Non-deterministic seeds](#c.-non-deterministic-seeds)
        -   [d. Analyzing the changes](#d.-analyzing-the-changes)
-   [Authors and contributions](#authors-and-contributions)

# Assignment Overview

The goal of this assignment is to get you acquainted with how the
different ways to analyze a file can change the results of the reads’
alignment against the reference genome. We will be using only one file:
**SRR12506919\_subset.fastq.gz**, that can be found under the following
path: **/projects/bmeg/A4/**. It will be divided into two parts:

Part 1: Analyzing different genome builds

Part 2: Ambiguity in reads mapping

# Important remarks

-   Remember to be as clear as possible with your answers.

-   Please make sure to explain your thought process behind your code
    and answers.

-   If you have used methods suggested on forums, websites, make sure to
    cite them.

-   If you have not been able to find the answer to a random bug with
    reasonable effort, please ask on Piazza! We are here to help, the
    assignments should be a safe environment for you to try new things
    and learn.

## 0. Getting ready

As always, before we get started we will install the tools required for
the assignment. This time, we only need to add:

-   LiftOver (<https://genome.sph.umich.edu/wiki/LiftOver>). These is a
    package designed to change files from a specific coordinates system
    (i.e., genome build) to another.

-   bedtools (<https://bedtools.readthedocs.io/en/latest/>). It’s a
    powerful tool to compare genomic positions.

``` bash
#?# Add liftOver to your conda environment created on A1, type the command you used below - 1 pt
conda activate bmeg591e
conda install -c bioconda ucsc-liftover 
 
#?# Add bedtools to your conda environment created on A1, type the command you used below - 1 pt
conda install -c bioconda bedtools 
```

## 1. Differences between genome builds

Your Professor informs you that the Information and Technology (IT)
department was able to recover part of your data from the server.
Unfortunately, they were not able to recover your pipelines or scripts.
Since you were using a pre-made index of the reference genome, you don’t
remember which genome build was used to map your sequences (hg19 or
hg38). You cannot decide if it would be a big deal to use different
genome builds for different alignments, at the end you could just make
sure they are in the same genome build when you compare them, right?
Thus, you decide to investigate if it would make a lot of difference to
use a different version to assess what varies when you align the same
reads to different genome-builds.

### a. SE alignment against hg38 and hg19

``` bash
## Pre-made indexes can be found here: 
## hg19 index: /projects/bmeg/indexes/hg19_bowtie2_index
## hg38 index: /projects/bmeg/indexes/hg38_bowtie2_index
## Recall that this is the fastq to be used throughout: /projects/bmeg/A4/SRR12506919_subset.fastq.gz

#?# Perform a single-end alignment using bowtie2 against the pre-made index of the hg38 genome build - 2 pt
cd assignments
mkdir assignment_4
cd
cp /projects/bmeg/A4/SRR12506919_subset.fastq.gz /home/elamoureux_bmeg22/assignments/assignment_4/
cp /projects/bmeg/A4/hg19ToHg38.over.chain.gz /home/elamoureux_bmeg22/assignments/assignment_4/
cd assignments/assignment_4/
gunzip SRR12506919_subset.fastq.gz
bowtie2 -x /projects/bmeg/indexes/hg38/hg38_bowtie2_index SRR12506919_subset.fastq -S bowtie2_out_hg38.sam

#?# Perform a single-end alignment using bowtie2 against the pre-made index of the hg19 genome build - 2 pt
bowtie2 -x /projects/bmeg/indexes/hg19/hg19_bowtie2_index SRR12506919_subset.fastq -S bowtie2_out_hg19.sam
```

### b. Making the files comparable

Before you can start comparing the files, you realize you need to
translate them to the same genomic coordinate system. If you don’t do
this and try to find mismatches between the files you will find a ton,
but that wouldn’t mean that the reads are mapping to different parts of
the genome, just that the coordinates are different (e.g. if there is a
1 bp insertion in hg38 relative to hg19, every base after that insertion
will have different coordinates). Hence, you go ahead and use your
favorite genome build: hg38. To translate hg19 to hg38, we need to do a
couple of steps:

-   Sam to Bam: Convert the alignemnt file (sam) to binary format (bam),
    this will facilitate the manipulaiton of the files and will decrease
    the disk space used substantially.

-   Bam to bed: Convert the bam alignment file to bed format, enabling
    the comparison of the genomic posisions where the reads mapped.

-   Change genomic coordinates: Use liftOver to change the alignment
    file made using the hg19 index to the hg38 coordinates of the
    genome.

``` bash
## Sam to Bam -------------
#?# Convert the SE alignment performed against hg19  (hg19 alignment) to bam, type the command you used below -1 pt
samtools view -S -b  bowtie2_out_hg19.sam > bowtie2_out_hg19.bam

#?# Convert the SE alignment performed against hg38 (hg38 alignment) to bam, type the command you used below -1 pt
samtools view -S -b  bowtie2_out_hg38.sam > bowtie2_out_hg38.bam

## Bam to bed -------------
## Tip: Look into the bedtools bamtobed command
#?# Use bedtools to convert the hg19 alignment bam file to bed format, type the command you used below - 1 pt 
bedtools bamtobed -i bowtie2_out_hg19.bam > bowtie2_out_hg19.bed

#?# Use bedtools to convert the hg38 alignment bam file to bed format, type the command you used below - 1 pt 
bedtools bamtobed -i bowtie2_out_hg38.bam > bowtie2_out_hg38.bed

## LiftOver --------------
#?# Use liftOver to change the hg19 alignment bed file to the hg38 coordinate system, type the command/s you used below - 2 pt
## To do this, you will need the "chain file": /projects/bmeg/A4/hg19ToHg38.over.chain.gz
## Tip: Look at the liftOver documentation! 
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
liftOver bowtie2_out_hg19.bed hg19ToHg38.over.chain.gz liftover_out_hg19.bed unlifted_hg19.bed
```

### c. Analyzing the differences

Now that both alignments are on the same coordinate system, they are
comparable and ready to be analyzed. What you really want to see how
individual reads mapped against the two genome builds. Did they map to
the same place or different places? To answer this, you need to sort
your bed files by read name so that you can identify which bed entries
in each file correspond to the same original read.

``` bash
#?# Using bash commands to sort the transformed hg19 alignment file bed alignment file by read name (column 4), type 
#?# the command you used below - 2 pt
## Tip: Look at the sort command!
sort -k 4 -o liftover_out_hg19_col.bed liftover_out_hg19.bed
 
#?# Using bash commands, sort the hg38 bed alignment file by read name (column 4), type the command you used below - 2 pt
## Tip: Look at the sort command!
sort -k 4 -o bowtie2_out_hg38_col.bed bowtie2_out_hg38.bed 
head -10 bowtie2_out_hg38_col.bed
```

You were really happy to see a visual representation of your data the
last time you talked to your advisor about mapping parameters. You
decide to give it a try this time with your merged bed file to answer
your two main questions:

-   How many reads are there per chromosome and does this differ between
    genome builds?

-   Do the reads mapped to the same genome region?

### d. Reads per chromosome

Before you get started, you discover that a labmate of yours was
comparing the number of reads per chromosome under different conditions
and they created a function to make this process more robust (function
is below). You are really happy that this seems like the perfect
function to plot the diferent number of reads per chromosome in the
different genome builds, but there is one problem. The bed files need to
be merged into one, before you can use the function. Plus, you realize
that the function is very poorly documented and your labmate is AWOL due
to midterms, so there is no way he can explain you how the function
works. Your Professor asks you to go through the function and document
as much as possible so future people can use it too (also because she
wants to make sure you know what you are doing).

``` bash
## Merging the files: ---------------
#?# Using the join command on bash, merge the two bed files, so they follow the following format: 
## read_id  chr_hg38  start_hg38  end_hg38  strand_hg38 chr_hg19  start_hg19  end_hg19  strand_hg19 
join -j 2 -o 1.4,1.1,1.2,1.3,1.6,2.1,2.2,2.3,2.6 <(sort -k 4 bowtie2_out_hg38.bed) <(sort -k 4 liftover_out_hg19.bed) > merged_bed.bed

#?# Type the command you used to merge the files below - 2pt 

#?# Use the head command to view the first 3 rows of your merged file, copy the output below: - 2pt 
# SRR12506919.1 chr16 1379150 1379301 + chr16 1379150 1379301 +
# SRR12506919.10 chr1 109467028 109467179 + chr1 109467028 109467179 +
# SRR12506919.100 chr19 3556952 3557103 + chr19 3556952 3557103 +

## Copy the merged bed file to your local computer for analysis
# On local command prompt: pscp elamoureux_bmeg22@orca1.bcgsc.ca:assignments/assignment_4/merged_bed.bed c:\users\eriks\bmeg_591e\merged_bed.bed
```

Now that you have the files in the right order, you move your files to
your local computer to work on your personal RStudio!

``` r
#?# Go through the function line by line using your merged bed file and your chosen parameters, as if it weren't a 
#?# function (e.g. set "merged_bed" to the data.frame containing your data, and run each line of the function (you will 
#?# also need to set the parameters)). Explain in a concise way how each line is changing the data. Use functions like head 
#?# and tail to visualize the data as it's changing. - 4 pt
merged_bed <- read.table("C:/Users/eriks/bmeg_591e/merged_bed.bed", quote="\"", comment.char="")

install.packages("reshape")
library(reshape)

## reads.per.chr:
# This function takes a merged bed file of two conditions A and B and gives a data.frame of 3 columns: Chr, variable 
# (condition), value (how many reads per chromosome are when using that condition)

## Parameters: 
# merged_bed: refers to the bed file you created on the previous section
# cols2compare=c(2,6): default is column 2 versus 6, which if you followed the format specified when you merged the files, 
# they should correspond to the chromosome column of each read for the two conditions (e.g., hg38 and hg19)
# type.a=c("hg38", "redo"): you should specify a string, that states what is condition A. Defaults are "hg38" and "redo"
# type.b=c("hg19", "noDet"): you should specify a string, that states what is condition B. Defaults are "hg19" and "noDet"

reads.per.chr <- function(merged_bed, cols2compare=c(2,6), type.a=c("hg38", "redo"), type.b=c("hg19", "noDet")){
  
  ## Create canonical chromosomes array to filter out contigs and scaffolds for simplicity
  canonical_chromosomes <- paste0("chr", 1:22)
  
  ## For column 1
  # The line below creates a dataframe of the 2nd column (chr_hg38)
  chr_subset <- merged_bed[,c(cols2compare[1])]
  # table() translates this dataframe into a table
  table_chrs1 <- table(chr_subset)

  ## For column 2
  # The line below creates a dataframe of the 6th column (chr_hg19)
  chr_subset <- merged_bed[,c(cols2compare[2])]
  # table() translates this dataframe into a table
  table_chrs2 <- table(chr_subset)
  
  # The below line creates a data frame with columns: [line_number chr_hg38 hg38_freq chr_hg19 hg19_freq]
  compare.df <- data.frame(column1=table_chrs1[names(table_chrs1) %in% canonical_chromosomes],
                           column2=table_chrs2[names(table_chrs2) %in% canonical_chromosomes])
  
  # The below line cuts off the chr_hg19 (3rd) column: [line_number chr_hg38 hg38_freq hg19_freq]
  compare.df <- compare.df[,c(1,2,4)]
  # The below line renames the columns: "chr", "hg38_reads", "hg19_reads"
  colnames(compare.df) <- c("chr",paste0(type.a, "_reads"), paste0(type.b, "_reads"))
  
  # melt() stacks the set of columns into single columns of data: line_number, chr, variable(hg38_reads at top, 
  # hg19_reads at bottom), value
  compare.df <- melt(compare.df)
  
  return(compare.df)
  
}
```

``` r
#?# Copy the files from the server to your local computer - 1pt
# On local command prompt: pscp elamoureux_bmeg22@orca1.bcgsc.ca:assignments/assignment_4/merged_bed.bed c:\users\eriks\bmeg_591e\merged_bed.bed

#?# Load your merged bed file into R using the *read.csv* function and save it into a data.frame
#?# Type the command you used below  - 1pt
merged_bed <- read.csv("C:/Users/eriks/bmeg_591e/merged_bed.bed", quote="\"", comment.char="")

## Change the column names of your merged bed data.frame to: 
# read_id  chr_hg38  start_hg38  end_hg38  strand_hg38 chr_hg19  start_hg19  end_hg19  strand_hg19 
#?# Type the command you used below:
colnames(merged_bed) <- c("read_id", "chr_hg38", "start_hg38", "end_hg38", "strand_hg38", "chr_hg19", "start_hg19", "end_hg19", "strand_hg19")

## Load the reshape2 library, install it if you don't already have it! 
## Tip: Use the "packages" tab on the left bottom screen 
library(reshape2)

#?# Run the reads.per.chr on your genome builds merged bed (previously loaded), specify all the parameters following the 
# instructions of the function, type the command used below: - 1.5 pt 
compare.df <- reads.per.chr(merged_bed, cols2compare=c(2,6), type.a=c("hg38", "redo"), type.b=c("hg19", "noDet"))

#?# How many reads were mapped to two different chromosomes? What percent of reads is this? Type the code and the answers 
# for each below. 2 pt
nrow(compare.df)
nrow(merged_bed)
# The number of rows in compare.df is 44, meaning there are 44 reads mapped to two different chromosomes. With a merged_bed
# length of 310909, 44/310909 = 0.014% 

## Using the output data.frame you got from running the reads.per.chr function on your merged bed, create a barplot that: 
## Uses the Chr column for the x-axis
## Useds the value (number of reads) column for the y-axis
## Uses the variable (conditions, also known as different genome builds in this case) column to "fill in" the color 
## Each build should have their own bar (next to each other), they shouldn't be stacked!!
#?# Type the command you used below: - 1.5 pt
b_plot <- ggplot(compare.df, aes(x=chr, y=value)) + geom_bar(stat='identity')
b_plot
```

Which chromosome has the biggest difference between reads? Which genome
build had more reads for this chromosome? Answer below - 1 pt Chromosome
1 had the greatest difference between reads (29040). The hg38 genome had
more reads for this chromosome.

### d. Reads position in the genome builds

``` r
## Using the start position of the reads on both genome builds, create a scatterplot using ggplot2 that: 
## Has the start in the hg38 genome build in the x-axis
## Has the start in the hg19 genome build in the y-axis
## Plots each chromosome in its own subplot (panel) (e.g. see facet_wrap())
## Plots only cases where both reads mapped to the same chromosome
#?# Type the command you used below: - 3 pt
scatter <- ggplot(compare.df, aes (x = variable, y = value)) + geom_point()
scatter + facet_wrap(~chr)
```

## 2. Ambiguity in reads mapping

You are glad that you have answered most of your burning questions about
read mapping and identified some of the things that can go wrong. So,
you decide to share your knowledge with your friend. They tell you that
they ran the SE alignment following your instructions and were about to
share their results, only to find that when repeating the alignment for
the same file their results changed! They come to you to help them with
your wisdom. Your vast experience leads you to believe that something
must have happened when the alignment was performed.

### a. Redoing the hg38 alignment

``` bash
#?# Re-run the SE alignment that you performed on 1a against the hg38 genome build, use exactly the same parameters, 
#?# just change the output name  - 0.5 pt
bowtie2 -x /projects/bmeg/indexes/hg38/hg38_bowtie2_index SRR12506919_subset.fastq -S bowtie2_out_hg38_rerun.sam

## Change both sam output to bam. Remember to remove the sam files right after it's done!
#?# Type the commands you used to convert the file below  - 0.5 pt
samtools view -S -b bowtie2_out_hg38_rerun.sam > bowtie2_out_hg38_rerun.bam

#?# Change the bam file to bed, using the betdools bedtobam function, type the command you used for the file below - 0.5 pt
bedtools bamtobed -i bowtie2_out_hg38_rerun.bam > bowtie2_out_hg38_rerun.bed

#?# Sort the file by read name (same as you did on part 1, using column 4), type the command you used below - 1 pt
sort -k 4 -o hg38_rerun_col.bed bowtie2_out_hg38_rerun.bed

## Because what you really want to see is if and what changed between these bowtie2 runs compared to your first run on 
## Part 1b, you decide to merge each new run file 
## with the original:

#?# Merge the "redo" bed file and the "original" hg38 alignment bed (from part 1c) using the join command, as in part 1c, this time follow this format: 1 pt
## read_id chr_ori  start_ori  end_ori  strand_ori chr_redo  start_redo  end_redo  strand_redo
## NOTE: Remember to save the output!
join -j 2 -o 1.4,1.1,1.2,1.3,1.6,2.1,2.2,2.3,2.6 <(sort -k 4 bowtie2_out_hg38.bed) <(sort -k 4 bowtie2_out_hg38_rerun.bed) > merged_bed_redo.bed


## Copy the merged bed file to your local computer for analysis
# On local command prompt: pscp elamoureux_bmeg22@orca1.bcgsc.ca:assignments/assignment_4/merged_bed_redo.bed c:\users\eriks\bmeg_591e\merged_bed_redo.bed
```

### b. Analyzing the ambiguity

Your last analysis on the differences between genome build turn out so
well, that you want to do the same. You have prepared the files so they
are in the same format as needed to run your labmate’s
*reads.per.chromosome* function, and are ready to see the graph.

``` r
#?# Load your merged bed file into R using the *read.csv* function and save it into a data.frame
#?# Type the command you used below  - 1pt
merged_bed_nds <- read.csv("C:/Users/eriks/bmeg_591e/merged_bed_redo.bed", quote="\"", comment.char="")

## Change the column names of your merged bed data.frame to: 
## read_id chr_ori  start_ori  end_ori  strand_ori chr_redo  start_redo  end_redo  strand_redo
#?# Type the command you used below:
colnames(merged_bed_redo) <- c("read_id", "chr_ori", "start_ori", "end_ori", "strand_ori", "chr_redo", "start_redo", "end_redo", "strand_redo")

#?# Run the reads.per.chr on your genome builds merged bed (previously loaded), specify all the parameters following the 
# instructions of the function, type the command used below: - 1.5 pt 
compare.df <- reads.per.chr(merged_bed_redo, cols2compare=c(2,6), type.a=c("hg38", "redo"), type.b=c("hg19", "noDet"))


#?# How many reads were mapped to two different chromosomes? What percent of reads is this? Type the code and the 
# answers for each below. 2 pt
nrow(compare.df)
nrow(merged_bed_redo)

# 44/617011 = 7.13e-03%

## Using the output data.frame you got from running the reads.per.chr function on your merged bed, do a barplot that: 
## Uses the Chr column for the x-axis
## Useds the value (number of reads) column for the y-axis
## Uses the variable (conditions, also known as different runs in this case) column to "fill in" the color 
## Each condition must have their own bar, they shouldn't be stacked!!

#?# Type the command you used below: - 1.5 pt
b_plot <- ggplot(compare.df, aes(x=chr, y=value)) + geom_bar(stat='identity')
b_plot


#?# Do you see differences among the number of reads per chromosome between the two runs? Answer yes or no - 0.5 pt
# There are slight variations in the number of reads per chromosome between the two runs. For example, chr2, chr19, ch20, ch21
# are different between the runs. However, overall these differences are minor. 
```

You are intrigued by the results of your graph and decide to go deeper
into the alignment to get a better idea of where the reads mapped within
the genome.

``` r
## Subtract the start position of the original bed from the start position of the redo for all the reads
#?# Type the command used below: - 0.5 pt
data$Diff<-(data$start_ori - data$start_redo)

## Use the *table* command to tabulate the results from the previous question. Ex. table(a-b)
#?# Type the command you used below: - 0.5 pt
class(data)
install.packages("data.table")          
library("data.table")
class(data)
apply(data, 2, function(c)sum(c!=0))

#?# What do you see? How many have a non zero difference in position start? - 0.5 pt
# I find that there are no positions with non zero position start differences. 

#?# Describe how would you expect a scatterplot comparing the start ends in both runs would look like - 0.5 pt
## x-axis: original run
## y-axis: re-run 
# Since there are no positions with non zero differences, I expect a scatter plot would looke like y = x. 
```

### c. Non-deterministic seeds

You are confused by your friend’s results, you don’t seem to have the
same problem. You ask her for the command she used to run her alignment
and you notice a key difference. She included the following flags:
**–non-deterministic –seed 3** . You decide to explore what is this
command doing and if it would change your data.

``` bash
#?# Re-run the SE alignment that you performed on 1a against the hg38 genome build, change the output name and add this 
#?# parameter:* --non-deterministic --seed 3 * - 1 pt
bowtie2 -x /projects/bmeg/indexes/hg38/hg38_bowtie2_index SRR12506919_subset.fastq -S bowtie2_out_hg38_nds.sam --non-deterministic --seed 3

## Change both sam outputs to bam. Remember to remove the sam files right after it's done!
#?# Type the commands you used to convert the file below  - 0.5 pt
samtools view -S -b bowtie2_out_hg38_nds.sam > bowtie2_out_hg38_nds.bam

#?# Change the bam file to bed, using the betdools bedtobam function, type the command you used for the file below  - 0.5 pt
bedtools bamtobed -i bowtie2_out_hg38_nds.bam > bowtie2_out_hg38_nds.bed

#?# Sort the files by read name (same as you did on part 1, using column 4), type the command you used below - 1 pt
sort -k 4 -o hg38_nds_col.bed bowtie2_out_hg38_nds.bed

#?# Merge the "non deterministic" bed file and the "original" hg38 alignment bed (part 1c) using the join command, as in 
#?# part 1c, this time follow this format: - 1 pt
## read_id  chr_ori  start_ori  end_ori  strand_ori chr_nonDet  start_nonDet  end_nonDet  strand_nonDet 
## NOTE: Remember to save the output!
join -j 2 -o 1.4,1.1,1.2,1.3,1.6,2.1,2.2,2.3,2.6 <(sort -k 4 bowtie2_out_hg38.bed) <(sort -k 4 bowtie2_out_hg38_nds.bed) > merged_bed_nds.bed

## Copy the merged bed file to your local computer for analysis
# On local command prompt: pscp elamoureux_bmeg22@orca1.bcgsc.ca:assignments/assignment_4/merged_bed_nds.bed c:\users\eriks\bmeg_591e\merged_bed_nds.bed
```

### d. Analyzing the changes

``` r
#?# Load your merged bed file onto R using the *read.csv* function and save it into a data.frame
#?# Type the command you used below  - 1 pt
merged_bed_nds <- read.csv("C:/Users/eriks/bmeg_591e/merged_bed_nds.bed", quote="\"", comment.char="")

## Change the column names of your merged bed data.frame to: 
## read_id  chr_ori  start_ori  end_ori  strand_ori chr_nonDet  start_nonDet  end_nonDet  strand_nonDet 
#?# Type the command you used below:
colnames(merged_bed_nds) <- c("read_id", "chr_ori", "start_ori", "end_ori", "chr_nonDet", "chr_nonDet", "start_nonDet", "end_nonDet", "strand_nonDet")


#?# How many reads were mapped to two different chromosomes? What percent of reads is this? Type the code and the answers 
#?# for each below. 2 pt
compare.df <- reads.per.chr(merged_bed_nds, cols2compare=c(2,6), type.a=c("hg38", "redo"), type.b=c("hg19", "noDet"))
nrow(compare.df)
nrow(merged_bed_nds)

# 44/617011 = 7.13e-03%

## Using the start position of the reads on both alignment runs do a scatterplot in ggplot that: 
## Has the start in the hg38_original in the x-axis
## Has the start in the hg38_nonDet build in the y-axis
## Plots each chromosome in its own subplot (panel) (e.g. see facet_wrap())
## Plots only cases where both reads mapped to the same chromosome
#?# Type the command you used below: - 2 pt
scatter <- ggplot(compare.df, aes (x = variable, y = value)) + geom_point()
scatter + facet_wrap(~chr)


#?# Explain why this changes when you add the --non-deterministic --seed 3 flags. What is are these flags doing? Why did 
#?# you get the result you saw in 2b?- 2 pt
## Tip: Look at the bowtie2 documentation!

# The --non-deterministic flag re-initializes bowtie2's pseudo-random generator for each read using the current time. 
# This means that bowtie2 will not necessarily report the same alighment for two identical reads. The --seed flag 
# specifies the seed for the pseudo-random number generator. 
# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-non-deterministic


#?# How do the number of off-diagonal reads and reads mapping to different chromosomes compare between where we mapped to 
#?# two different genome versions (and then lifted over), versus the use of non-deterministic alignment? What fraction 
#?# of reads that you found aligned to different chromsomes when using hg19 vs hg38 result from the differences 
#?# between these two versions? - 3 pts

# Using hg38, there was 92.11% overall alignment rate, and using hg19 there was 91.24% overall alignment rate. Using
# --non-deterministic and --seed 3 there was 92.07% overall alignment rate. Most of the alignment differnces between hg19,  
# hg38, and --non-deterministic and --seed 3 result from the randomness/stochasticity for ambiguously mapped reads. 
```

Please knit your *Rmd* file to github\_document (*md document*) and
include both in your submission.

Successful knitting to github\_document - 2 pts

# Authors and contributions

Following completion of your assignment, please fill out this section
with the authors and their contributions to the assignment. If you
worked alone, only the author (e.g. your name and student ID) should be
included.

Authors: Erik Lamoureux

Contributions: Assignment completed individually.
