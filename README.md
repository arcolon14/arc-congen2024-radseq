# RADseq data analysis with *Stacks* 2 hands-on lecture

Documentation and examples for the RADseq analysis with *Stacks* 2 hands-on 
exercise at ConGen 2024, on Aug 27, 2024.

(C) Angel G. Rivera-Colón (<ariverac@uoregon.edu>)

## Associated readings

The RADseq data used for this exercise comes from the following publication:

> Long, KM, Rivera-Colón, AG, Bennett, KFP, *et al*. (2024) **Ongoing introgression 
> of a secondary sexual plumage trait in a stable avian hybrid zone**. *Evolution*, 
> 2024, qpae076, [DOI: 10.1093/evolut/qpae076](https://doi.org/10.1093/evolut/qpae076)

The data is used here with the permission of the authors. Thanks to Kira M. Long!

This analysis follows the general guidelines described in the *Stacks* 2 
protocol manuscript:

> Rivera-Colón, AG, Catchen, JM (2022). **Population Genomics Analysis with RAD, 
> Reprised: *Stacks* 2**. In: Verde, C., Giordano, D. (eds) *Marine Genomics. 
> Methods in Molecular Biology*, vol 2498. Humana, New York, NY. 
> [DOI: 10.1007/978-1-0716-2313-8_7](https://doi.org/10.1007/978-1-0716-2313-8_7)

For an algorithmic description of the *Stacks* 2 software, please check the 2019 
*Mol Ecol* manuscript:

> Rochette, NC, Rivera-Colón, AG, Catchen, JM. (2019) ***Stacks* 2: Analytical methods 
> for paired-end sequencing improve RADseq-based population genomics**. *Molecular 
> Ecology*, 28, 4737–4754. [DOI: 10.1111/mec.15253](https://doi.org/10.1111/mec.15253)

For more information regarding PCR duplicates, please check the 2023 *Mol Ecol Resour* 
publication:

> Rochette, NC, Rivera-Colón, AG, Walsh, J, *et al*. (2023) **On the causes, 
> consequences, and avoidance of PCR duplicates: Towards a theory of library 
> complexity.** *Molecular Ecology Resources*, 23, 1299–1318. 
> [DOI: 10.1111/1755-0998.13800](https://doi.org/10.1111/1755-0998.13800)

For information on the download, installation, and documentation of *Stacks*, 
please visit the official [website](https://catchenlab.life.illinois.edu/stacks/).

## Repository for the exercise

A copy of this document and the associated data can be found in:
[https://github.com/arcolon14/arc-congen2024-radseq](https://github.com/arcolon14/arc-congen2024-radseq).

## Preparing the environment

**NOTE:** the directory hierarchy in the commands below refer to the 
ConGen2024 server.

Make a directory for the Stacks hands-on assignment.

```sh
$ mkdir stacks-radseq
$ cd stacks-radseq/
```

Copy the raw data from the instructors directory

```sh
$ cp /data/instructor_materials/Angel_Rivera-Colon/2024/arc-radseq-data.congen24.tar.gz .
```

Uncompress this directory

```sh
$ tar xvf arc-radseq-data.congen24.tar.gz
```

Check the contents of the directory

```sh
$ ls arc-radseq-data.congen24/*
  alignments:
  CG_10_067.bam   PR_09_096.bam   RO_05_206.bam   SO_03_468.bam
  CG_10_069.bam   PR_09_102.bam   RO_05_221.bam   SO_03_469.bam
  CG_10_070.bam   PR_09_111.bam   RO_05_222.bam   SO_03_471.bam
  CG_10_093.bam   PR_09_112.bam   RO_05_223.bam   SO_03_477.bam
  CG_10_094.bam   PR_09_115.bam   RO_05_254.bam   SO_03_482.bam
  FC_04_510.bam   QP_06_190.bam   RU_08_203.bam   SS_02_071.bam
  FC_04_511.bam   QP_06_191.bam   RU_08_204.bam   SS_02_081.bam
  FC_04_514.bam   QP_06_195.bam   RU_08_205.bam   SS_02_082.bam
  FC_04_550.bam   QP_06_226.bam   RU_08_245.bam   SS_02_085.bam
  FC_04_553.bam   QP_06_228.bam   RU_08_249.bam   SS_02_090.bam

  catalog:
  catalog.calls      catalog.fa.gz   gstacks.log.distribs
  catalog.chrs.tsv   gstacks.log

  info:
  barcodes.tsv   popmap.tsv   popmap_catalog.tsv

  processed-samples:
  CG_10_067.1.fq.gz   PR_09_112.2.fq.gz   RU_08_205.1.fq.gz
  CG_10_067.2.fq.gz   PR_09_115.1.fq.gz   RU_08_205.2.fq.gz
  CG_10_069.1.fq.gz   PR_09_115.2.fq.gz   RU_08_245.1.fq.gz
  CG_10_069.2.fq.gz   QP_06_190.1.fq.gz   RU_08_245.2.fq.gz
  CG_10_070.1.fq.gz   QP_06_190.2.fq.gz   RU_08_249.1.fq.gz
  CG_10_070.2.fq.gz   QP_06_191.1.fq.gz   RU_08_249.2.fq.gz
  CG_10_093.1.fq.gz   QP_06_191.2.fq.gz   SO_03_468.1.fq.gz
  CG_10_093.2.fq.gz   QP_06_195.1.fq.gz   SO_03_468.2.fq.gz
  CG_10_094.1.fq.gz   QP_06_195.2.fq.gz   SO_03_469.1.fq.gz
  CG_10_094.2.fq.gz   QP_06_226.1.fq.gz   SO_03_469.2.fq.gz
  FC_04_510.1.fq.gz   QP_06_226.2.fq.gz   SO_03_471.1.fq.gz
  FC_04_510.2.fq.gz   QP_06_228.1.fq.gz   SO_03_471.2.fq.gz
  FC_04_511.1.fq.gz   QP_06_228.2.fq.gz   SO_03_477.1.fq.gz
  FC_04_511.2.fq.gz   RO_05_206.1.fq.gz   SO_03_477.2.fq.gz
  FC_04_514.1.fq.gz   RO_05_206.2.fq.gz   SO_03_482.1.fq.gz
  FC_04_514.2.fq.gz   RO_05_221.1.fq.gz   SO_03_482.2.fq.gz
  FC_04_550.1.fq.gz   RO_05_221.2.fq.gz   SS_02_071.1.fq.gz
  FC_04_550.2.fq.gz   RO_05_222.1.fq.gz   SS_02_071.2.fq.gz
  FC_04_553.1.fq.gz   RO_05_222.2.fq.gz   SS_02_081.1.fq.gz
  FC_04_553.2.fq.gz   RO_05_223.1.fq.gz   SS_02_081.2.fq.gz
  PR_09_096.1.fq.gz   RO_05_223.2.fq.gz   SS_02_082.1.fq.gz
  PR_09_096.2.fq.gz   RO_05_254.1.fq.gz   SS_02_082.2.fq.gz
  PR_09_102.1.fq.gz   RO_05_254.2.fq.gz   SS_02_085.1.fq.gz
  PR_09_102.2.fq.gz   RU_08_203.1.fq.gz   SS_02_085.2.fq.gz
  PR_09_111.1.fq.gz   RU_08_203.2.fq.gz   SS_02_090.1.fq.gz
  PR_09_111.2.fq.gz   RU_08_204.1.fq.gz   SS_02_090.2.fq.gz
  PR_09_112.1.fq.gz   RU_08_204.2.fq.gz

  raw-reads:
  MAVIRAD2_NoIndex_L002_R1_001.fastq.gz
  MAVIRAD2_NoIndex_L002_R2_001.fastq.gz

  stacks-logs:
  denovo        populations_base     process_radtags
  gstacks_ref   populations_strict
```

## Breakdown of the available data

This data directory contains several subdirectories describing this RADseq dataset at 
different stages of the analysis.

### `raw-reads`

The `raw-reads` directory contains the raw sequencing files (in `FASTQ` format) from 
this library. As described in [Long et al. 2024](https://doi.org/10.1093/evolut/qpae076), 
a total of 130 manakin (genus *Manacus*) samples were processed using a single-digest 
RADseq library ([Baird et al. 2008](https://doi.org/10.1371/journal.pone.0003376); 
[Etter et al. 2011](https://doi.org/10.1371/journal.pone.0018561)) using the *SbfI* 
restriction enzyme and paired-end sequenced on an Illumina NovaSeq6000, generating 
2x150bp reads.

```sh
$ ls raw-reads/
  MAVIRAD2_NoIndex_L002_R1_001.fastq.gz   MAVIRAD2_NoIndex_L002_R2_001.fastq.gz
```

**Note:** Due to limited time, the files here only contain a fraction (~1 million 
reads) of the total reads in the real library.

### `processed-samples`

```sh
$ ls processed-samples/
  CG_10_067.1.fq.gz
  CG_10_067.2.fq.gz
  CG_10_069.1.fq.gz
  CG_10_069.2.fq.gz
  ...
  SS_02_082.2.fq.gz
  SS_02_085.1.fq.gz
  SS_02_085.2.fq.gz
  SS_02_090.1.fq.gz
  SS_02_090.2.fq.gz
```

The contents of the `processed-samples` directory represent the sequencing reads 
(in `FASTQ` format) for 40 manakin individuals. They are the demultiplexed, cleaned 
reads generated by running `process_radtags` on the raw data. In the excersices below, 
we are just using 40 out of the total 130 samples.

See `stacks-logs/process_radtags/process_radtags.MAVI2.log` for additional details 
about how these reads were processed.

**NOTE:** Due to limited time, the files here only contain a subset of the total reads 
for each individual, those corresponding to a single chromosome-level scaffold.

### `alignments`

```sh
$ ls alignments/
  CG_10_067.bam
  CG_10_069.bam
  CG_10_070.bam
  ...
  SS_02_082.bam
  SS_02_085.bam
  SS_02_090.bam
```

The `alignments` directory contains the aligned reads (in `bam` format) from 40 manakin
individuals. The data for these 40 samples was run previously through `process_radtags` 
and aligned to the golden-collared manakin (*Manacus vitellinus*) RefSeq assembly (NCBI 
accession [GCF_001715985.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001715985.3/)) 
using *BWA* `mem` ([Li 2013](https://doi.org/10.48550/arXiv.1303.3997)) and processed 
with the `view` and `sort` commands from *SAMtools* 
([Li et al. 2009](https://doi.org/10.1093/bioinformatics/btp352)).

For example:

```sh
$ bwa mem manVit.db PR_09_096.1.fq.gz PR_09_096.2.fq.gz | \
      samtools view -b -h | \
      samtools sort -o PR_09_096.bam
```

Remember, `bam`s store the alignments in binary format (see `bam` 
[documentation](https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf)). In order to 
view the alignments as text, we have to run the *SAMtools* `view` command:

```sh
$ samtools view PR_09_096.bam | less
```

For more information on how to process and view alignments in the `bam` format, see the 
*SAMtools* `view` [documentation](http://www.htslib.org/doc/samtools-view.html).

**NOTE:** Due to limited time, the `bam` files provided here have been filtered to 
include data from only one chromosome-level scaffold.

### `catalog`

```sh
$ ls catalog/
  catalog.calls     gstacks.log
  catalog.chrs.tsv  gstacks.log.distribs
  catalog.fa.gz
```

The files in the `catalog` directory describe a reference-based *Stacks* catalog. After 
aligning the reads to the *M. vitellinus* reference, the alignments were processed by 
the `gstacks` software to assemble RAD loci, genotype individuals, and remove 
PCR-duplicate reads.

See `stacks-logs/gstacks_ref/gstacks.log` for additional details regarding the 
generation of this catalog.

**NOTE:** For the sake of time, this catalog only contains a subset of the real data, 
representing the loci corresponding to a single chromosome-level scaffold (~5000 raw loci).

### `info`

```sh
$ ls info/
  barcodes.tsv  popmap.tsv  popmap_catalog.tsv
```

The `info` directory contains several miscellaneous files used across the *Stacks* pipeline.

#### Barcodes

```sh
$ head info/barcodes.tsv
  TCGCCAG<tab>SS_02_071
  CGTGTGA     SS_02_081
  TGTGCGG     SS_02_082
  CGATACC     SS_02_085
  ACATAGA     SS_02_090
  TGACATG     SO_03_468
  GAGTTGC     SO_03_469
  CGTCCGC     SO_03_471
  AGTGAAG     SO_03_477
  GGCCTCG     SO_03_482
```

The file `barcodes.tsv` is a tab-delimited file containing the individual barcode sequences 
used to demultiplex each sample from the raw reads in `process_radtags`. Since this library 
was sequenced using the standard single-digest RADseq protocol 
([Baird et al. 2008](https://doi.org/10.1371/journal.pone.0003376); 
[Etter et al. 2011](https://doi.org/10.1371/journal.pone.0018561)), it only uses a single 
7-bp barcode per individual. See the *Stacks* 
[manual](https://catchenlab.life.illinois.edu/stacks/manual/#specbc) for additional 
information about the barcodes specification.

**NOTE:** Due to limited time, this file only contains barcodes for 40 out of the total 
130 manakin samples sequenced in the library.

#### The popmap files

The `info` directory also contains two population map (i.e., "popmap") files. Each popmap is 
a tab-delimited file describing the group assignment of each sequenced individual.

The first popmap file (`popmap.tsv`), contains the information for 40 manakins. Each is 
assigned to a specific sampling site along a transect in Panama, as described by 
[Long et al. 2024](https://doi.org/10.1093/evolut/qpae076). In total, we have 40 samples 
across 8 total populations, each representing a different sampling site.

```sh
$ cat info/popmap.tsv
  SS_02_071<tab>02SS
  SS_02_081     02SS
  SO_03_468     03SO
  SO_03_469     03SO
  PR_09_096     09PR
  PR_09_102     09PR
  CG_10_067     10CG
  CG_10_069     10CG
  ...
```

As described in the *Stacks* [manual](https://catchenlab.life.illinois.edu/stacks/manual/#popmap), 
the first column contains the individual's ID, while the second column contains the ID of 
its respective group or population. In this example, individual `SS_02_081` is assigned to 
population `02SS` (shortened ID for transect site #2, San San Drury in 
[Long et al. 2024](https://doi.org/10.1093/evolut/qpae076)), while individual `PR_09_096` 
belongs to population `09PR` (site #9, Palma Real).

The second popmap, `popmap_catalog.tsv`, only contains a subset of 15 individuals, all 
assigned to the sample population (labelled `opt`).

```sh
$ cat info/popmap_catalog.tsv
  SS_02_071<tab>opt
  SS_02_081     opt
  PR_09_096     opt
  PR_09_102     opt
  CG_10_067     opt
  CG_10_069     opt
  ...
```

This smaller subset of the data will be used for some specific applications, e.g., the 
parameter optimization of the *de novo* assembly of loci.

### `stacks-logs`

The `stacks-logs` directory contains several subdirectories, each containing the log and 
distribution files for several, pre-run steps of the *Stacks* pipeline.

```sh
$ ls stacks-logs/*
  stacks-logs/denovo:
  cov_per_sample.tsv   denovo_map.log

  stacks-logs/gstacks_ref:
  gstacks.log   gstacks.log.distribs

  stacks-logs/populations_base:
  populations.log   populations.log.distribs

  stacks-logs/populations_strict:
  populations.log   populations.log.distribs

  stacks-logs/process_radtags:
  process_radtags.MAVI2.log
```

For the sake of time, in this exercise we are running the *Stacks* pipeline on reduced 
datasets. While fast, the logs from these reduced runs might not always be the most 
informative. Thus, this `stacks-logs` directory contains logs from the real, full dataset 
and will be used to explore the analysis in more detail. 

## Running the pipeline

### Inspect the raw data

Before proceeding, it is always important to manually inspect the raw data. We have hundreds 
of millions of reads, so it is impossible to *completely* check our files. However, it is 
always good practice to check that it follows the the design of our RAD library.

Here, we are looking at the first few lines of the forward reads.

```sh
$ zcat arc-radseq-data.congen24/raw-reads/MAVIRAD2_NoIndex_L002_R1_001.fastq.gz | head -n 12
  @A00419:507:HM23CDRXY:2:2101:1108:1000 1:N:0:1
  TCATGAGTGCAGGACTGACACATTTTCCAGGGCAGTGTCCACAGAACATCATCATGCCAGGGCT...
  +
  FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF:FFFFFFFF...
  @A00419:507:HM23CDRXY:2:2101:1090:1000 1:N:0:1
  CTCGAAGTGCAGGAAAACCAGCCCCTGCCACGCTCANCAGGCACGCCTGGCATCTCATTTTCCT...
  +
  :FFFFFFF,F:FFFFF:FFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF...
  @A00419:507:HM23CDRXY:2:2101:1597:1000 1:N:0:1
  GTAGCTTTGCAGGAAGCCGAGGCAATGTGCCCATGTCACCAAGGCTGGGGATGAGGGTCCGCAT...
  +
  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFF...
```

The raw reads in `FASTQ` format ([Cock et al. 2009](https://doi.org/10.1093/nar/gkp1137)). 
The records on a `FASTQ` file are broken down in four lines (or fields):

1. The sequence ID line, starting with an `@`. This can contain information about 
   index barcodes.
2. The sequence line, containing the nucleotide sequence.
3. The third field starts with a `+`. It is often empty, but it can also repeat the 
   information in the header line.
4. The quality (or `PHRED`) scores (as indicated by `F`, `:`, `,`, etc.)

This library was constructed using a single-digest RADseq protocol 
([Baird et al. 2008](https://doi.org/10.1371/journal.pone.0003376); 
[Etter et al. 2011](https://doi.org/10.1371/journal.pone.0018561)),
so we expect the first 7 bases to be the in-line barcode. This means our first
record corresponds to barcode `TCATGAG`, the second record to barcode `CTCGAAG`,
the third to barcode `GTAGCTT`, etc.

After the barcode, we expect to see the cutsite. This library was constructed 
with the *SbfI* restriction enzyme, so we expect to see the corresponding 
overhang (`TGCAGG`).

Let's now look at the paired reads.

```sh
$ zcat arc-radseq-data.congen24/raw-reads/MAVIRAD2_NoIndex_L002_R2_001.fastq.gz | head -n 12
  @A00419:507:HM23CDRXY:2:2101:1108:1000 2:N:0:1
  AGTGACAATAGCTGTAATAGCCTGGTTCTTGGTCAGTGGTGACCTGCTAAGACAGGGATTGAAT...
  +
  F:FFFFF:,:F,FFFFFFFFFFFFFFF,F:FFF:FFFFFF,,FFF::FFFF,FFFFFFFFF:FF...
  @A00419:507:HM23CDRXY:2:2101:1090:1000 2:N:0:1
  CAGGAAGGCAATGGGCTCTCCAAGGGTGCTCCAAAGGGCTTTCAACCTCCATAATGTCTCATTT...
  +
  :FF:FF:F,FFFFFFFFFF:FFF,FF:FFFFFFFF,F,FF,FF:F::,FF:FFFF,:F:FF,FF...
  @A00419:507:HM23CDRXY:2:2101:1597:1000 2:N:0:1
  AGAGCCGGTGGGAAAGTCTGTGGCTGTTGCTATTTTCTATCTGTTATCAAGGCCTCTTATCACT...
  +
  FFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF,FFFFFF:FFFFFFFFFFFF:F:FF:FFFF:...
```

Following the library configuration (in this case a paired-end, 
single-digest RADseq), we do not expect to see a barcode nor a restriction 
enzyme cutsite.

### Processing raw data

After inspecting the raw data, we can start processing it using the *Stacks*
`process_radtags` program. We can demultiplex the combined reads into the 
data for the individual samples, and we can filter the reads to remove 
those of low quality (as defined by the `PHRED` score).

First, we will create a new directory to store the data for the processed 
samples.

```sh
$ mkdir processed_samples
```

The raw paired-end reads are located in the `arc-radseq-data.congen24/raw-reads/` 
directory. The barcodes describing the individual assignment of each sample are 
in `arc-radseq-data.congen24/info/barcodes.tsv`. We want to remove any reads 
with uncalled bases (`N`s), we want to rescue barcodes and cutsites, and remove 
reads with low quality. Also, we want to specify that this library was generated 
using a single restriction enzyme, *SbfI*. Other parameters can be left default, 
for example the library uses `inline-null` barcodes.

For more information on the program, please see the `process_radtags` 
[documentation](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php).

Here is an example of our `process_radtags` command.

```sh
$ process_radtags \
      --in-path ./arc-radseq-data.congen24/raw-reads/ \
      --paired \
      --barcodes ./arc-radseq-data.congen24/info/barcodes.tsv \
      --out-path ./processed_samples/ \
      --clean \
      --rescue \
      --quality \
      --renz-1 sbfI
```

Upon completion, `process_radtags` will report statistics about the filtered reads. 
Remember, we are running a subset of the original reads, so the results might not be 
the most representative. We can look at the results from a full run in the log files 
located in `arc-radseq-data.congen24/stacks-logs/process_radtags/`.

We can inspect the logs from `process_radtags` using the *Stacks* 
`stacks-dist-extract` utility, used to extract specific subsets of information from 
log and distribution files (see official 
[documentation](https://catchenlab.life.illinois.edu/stacks/comp/stacks_dist_extract.php)).

```sh
  $ stacks-dist-extract process_radtags.MAVI2.log total_raw_read_counts
    Total Sequences        822561114
    Barcode Not Found      37391352   4.5%
    Low Quality            866995     0.1%
    Poly-G Runs            2982449    0.4%
    RAD Cutsite Not Found  4241295    0.5%
    Retained Reads         777079023  94.5%
    Properly Paired        384920302  93.6%
```

From this, we can observe that 94.5% of reads in this library were retained, 93.6% 
of which were properly paired (i.e., both forward and paired reads were kept). Most 
of the reads removed were due to not having a proper barcode sequence (4.5%); 
however, losses due to low quality and missing cutsites were minimal (0.1% and 
0.5%, respectively).

We can also look at a per-sample breakdown of the reads. Once again, we can access this 
information using the `stacks-dist-extract` command.

```sh
# Removing some of the columns for readability
$ stacks-dist-extract process_radtags.MAVI2.log per_barcode_raw_read_counts | \
cut -f2,3,7,8,9,10 
  Filename     Total      Pct Retained   Pct Properly Paired   Pct Total Reads
  RU_080_003   2718084    98.9%          97.9%                 0.3%
  RU_080_005   5905912    98.8%          97.6%                 0.7%
  RU_080_007   5016924    99.1%          98.2%                 0.6%
  RU_080_009   6368388    98.7%          97.6%                 0.8%
  RU_080_010   5832120    98.9%          97.9%                 0.7%
  RU_080_013   5640074    98.9%          97.9%                 0.7%
  RU_080_021   4250646    97.9%          95.9%                 0.5%
  RU_080_244   14129520   99.2%          98.5%                 1.7%
  RU_080_245   7057982    99.1%          98.4%                 0.9%
  ...
```

For simplicty, we are just looking at some of the columns here (selected with the 
`cut` command). For each sample, we can observe the proportion of reads retained, 
and properly paired, as well as the total proportion of reads of that sample in 
the total library. Moreover, the output of `stacks-dist-extract` is by default 
tab-delimited, so it can be redirected to a file. This new file can be loaded into 
a R/Excel/Google Sheets/etc. for additional visualization and analysis.

```sh
$ stacks-dist-extract process_radtags.MAVI2.log per_barcode_raw_read_counts > per_barcode_raw_read_counts.tsv
```

### Create a *de novo* catalog

Next, we will see how to assemble a *Stacks* catalog *de novo*--assembling loci 
and genotyping individuals without the need of a reference genome. The *de novo* 
pipeline in *Stacks* consists of several programs used to 1. cluster reads into 
loci in each individual (`ustacks`), 2. combine the loci of several individuals 
into a catalog of loci of the metapopulation (`cstacks`), and 3. match the loci 
of all sequenced individuals to this catalog (`sstacks`). Loci are then further 
further processed and genotyped in `gstacks`. Since this involves running multiple 
programs in *Stacks*, some of them several times across all the individuals, it 
is common to run the *de novo* pipeline using one of the wrapper scripts, in this 
case `denovo_map.pl`, which automates this process.

#### Optimizing a *de novo* catalog generation

Before generating a *de novo* catalog on the whole data, a common question is 
which parameter should be use to generate this catalog. These parameters define 
the thresholds used to cluster alleles into loci in a single individual, and then 
the threshold to match loci across individuals. Improper selection of these 
thresholds means that, e.g., different alleles can be erroneously processed as 
separate loci, or non-homologuos loci can be combined as a single locus across 
sequenced individuals.

Within *Stacks*, there are to main options that control these two threshols 
during *de novo* catalog generation:

1. `M`: used in `ustacks` to define the number of mismatches allowed between 
   the stacks (in this case alleles) within an individual.
2. `n`: used in `cstacks` to define the number of mismatches allowed between 
   the statcks (in this case loci) between individuals.

There are several methods used for optimize these two parameters in *Stacks*, the 
most common one being the "r80" method, described by 
[Paris et al. 2017](https://doi.org/10.1111/2041-210X.12775) and further detailed 
in [Rochette & Catchen 2017](https://doi.org/10.1038/nprot.2017.123). This method 
sets the values of `M` by maximizing the number of loci that are recovered in 80% 
of individuals in the analysis (defined by the `r` parameter in `populations`). 
The value of `n` is set by default equal to that of `M`, with the assumption that 
the genetic distance within and between individuals is roughly the same. The 
general idea is to chose the paramter that lead to the highest number of usable 
loci and the least proportion of missing data.

**NOTE:** The r80 method just provides general guideliness and is just one of 
many valid alternatives. Some dataset might "respond" differently to contransting 
optimization approaches, and specific experimental designs might seek to maximize 
these parameters differently. For example, when working in phylogenetics, it might 
be more important to fine tune the number of mistmatches between individuals (`n`) 
of different species. Other datasets might have several "optimal" valies of `M`. 
There is no "one size fits all" when discussing parameter optimization.

To run the r80 method on our dataset, we will run the `denovo_map.pl` wrapper 
several times, changing the value of `M` within a range of values. Here, we will 
do five simple runs over a smaller range of values (5 to 9). On real dataset, testing over 
a larger range of values is recommended (e.g., 1 to 12).

Additionally, will generate a catalog using only a subset of samples, with the 
assumption that these samples are representative of the whole dataset. This is 
often recommended to provide a general speedup of the process. For this purpose, we 
have provided a secondary popmap (`info/popmap_catalog.tsv`) just containing 15 
manakin samples, all assigned to a single group labelled `opt`.

We will also generate an output directory to store the outputs for all opmtimization
runs. 

```sh
$ mkdir denovo_param_opt
```

Then, we will name this directory for the specific run according to the value 
of `M`, 5 in this example.

```sh
$ mkdir denovo_param_opt/param_opt_M5
```

We then want to run `denovo_map.pl`, specfying the value of `M`, the smaller 
popmap, and our new output directory. Also, we want to specify that our input 
data contains paired reads, that we want to remove PCR duplicate reads, and that 
we only want to keep loci/SNPs present in 80% of samples per-population.

As our input data, we will specify the processed the FASTQ files in 
`arc-radseq-data.congen24/processed-samples/`. Note that these are the processed 
samples provided in the shared data, not the ones we just generated using 
`processed_radtags`. For more information, see the `denovo_map.pl` 
[documentation](https://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php).

Here is an example command:

```sh
$ denovo_map.pl \
  --samples ./arc-radseq-data.congen24/processed-samples/ \
  --popmap ./arc-radseq-data.congen24/info/popmap_catalog.tsv \
  --out-path ./denovo_param_opt/param_opt_M5 \
  -M 5 \
  -n 5 \
  -r 0.8 \
  --paired \
  --rm-pcr-duplicates -r 0.8
```

After `denovo_map.pl` completes, we can check the number of retained loci by 
checking at the log from `populations`. We can easily find this number on the 
command line by searching for a line starting with "Kept" using the `grep` 
program. This will give us information by the total number of loci and variant 
sites retained after applying the 80% missing data filter.

```sh
$ cat param_opt_M5/populations.log | grep '^Kept'
  Kept 2681 loci, composed of 1724932 sites; 
    461635 of those sites were filtered, 
    17346 variant sites remained.
```

Additionally, we can look at the data in the SUMSTATS table `
(`populations.sumstats.tsv`), which will provide us information of the number 
of total remained loci that polymorphic (i.e., that contain variant sites). 
We can contain this by counting the number of loci seen on the first column 
of the SUMSTATS.

```sh
$ cat param_opt_M5/populations.sumstats.tsv | \
    grep -v '^#' | cut -f 1 | sort -n -u | wc -l
  2599
```

After running this for the full range of `M`, we obtain the following results:


| `M` | `n` | Kept Loci | Polymorphic Loci | % Polymorphic | Variant Sites | Change in Loci |
| --- | --- | --------- | ---------------- | ------------- | ------------- | -------------- |
| 1   | 1   | 2,325     | 2,204            | 94.80%        | 10,818        | --             |
| 2   | 2   | 2,534     | 2,439            | 96.25%        | 14,226        | 235            |
| 3   | 3   | 2,641     | 2,551            | 96.59%        | 16,124        | 112            |
| 4   | 4   | 2,668     | 2,583            | 96.81%        | 17,014        | 32             |
| 5   | 5   | 2,681     | 2,599            | 96.94%        | 17,346        | 16             |
| 6   | 6   | 2,690     | 2,606            | 96.88%        | 17,488        | 7              |
| 7   | 7   | 2,689     | 2,605            | 96.88%        | 17,591        | -1             |
| 8   | 8   | 2,692     | 2,609            | 96.92%        | 17,668        | 4              |
| 9   | 9   | 2,693     | 2,609            | 96.88%        | 17,650        | 0              |
| 10  | 10  | 2,691     | 2,609            | 96.95%        | 17,655        | 0              |
| 11  | 11  | 2,689     | 2,606            | 96.91%        | 17,664        | -3             |
| 12  | 12  | 2,686     | 2,603            | 96.91%        | 17,628        | -3             |

Following the guidelines from the r80 method, we see that the increase of new loci 
plateaus at `M=n=6`. After that, there is a minimal change in the number of loci that 
are gained or lossed after increasing the `M=n` threshold. An `M` of 6 is likely the
optimal value for this parameter under these *specific* conditions.

##### Optional: Optimizing `n`

This optimization is likely sufficient. However, if we wanted to optimize this 
dataset further, we could separately optimize `n` indepdendent of `M`. Depending 
on the genetic distance between the individuals, it might be reasonable to expect 
`n` (the mistmatches *between* individuals) to be larger than `M` (the mismatches 
*within* individuals).

To do this, we can test of a range of `n` values against a fixed value of `M` (6 in 
this example). We could start at `n` equal to `M-1` (an unlikely value, but we can 
use it to "anchor" our results) and end, for example, at `n` equals two times `M` 
(again, an unlikely value, but provides a good comparison).

Like before, we will make a new output directory and run `denovo_map.pl` with the 
corresponding `M` and `n` values.

```sh
$ mkdir ./denovo_param_opt/param_opt_M6n8/

$ denovo_map.pl \
  --samples ./arc-radseq-data.congen24/processed-samples/ \
  --popmap ./arc-radseq-data.congen24/info/popmap_catalog.tsv \
  --out-path ./denovo_param_opt/param_opt_M6n8/ \
  -M 6 \
  -n 8 \
  -r 0.8 \
  --paired \
  --rm-pcr-duplicates
```

Once each runs complete, we can extract the values of kept loci and variant 
sites from the `populations.log` and `populations.sumstats.tsv` files, 
as shown previously.

| `M` | `n` | Kept Loci | Polymorphic Loci | % Polymorphic | Variant Sites | Change in Loci |
| --- | --- | --------- | ---------------- | ------------- | ------------- | -------------- |
| 6   | 5   | 2,686     | 2,602            | 96.87%        | 17,424        | --             |
| 6   | 6   | 2,690     | 2,606            | 96.88%        | 17,488        | 4              |
| 6   | 7   | 2,692     | 2,608            | 96.88%        | 17,578        | 2              |
| 6   | 8   | 2,694     | 2,610            | 96.88%        | 17,640        | 2              |
| 6   | 9   | 2,692     | 2,608            | 96.88%        | 17,649        | -2             |
| 6   | 10  | 2,693     | 2,609            | 96.88%        | 17,642        | 1              |
| 6   | 11  | 2,692     | 2,608            | 96.88%        | 17,638        | -1             |
| 6   | 12  | 2,693     | 2,608            | 96.84%        | 17,636        | 0              |

From these results, we can see that the varying the values of `n` does not 
generate major differences in the number of loci and variant sites kept. 
Increasing `n` to 8 does provide the largest number of polymorphic loci kept, 
but the differences between runs are relatively small. In this example, we 
*could* proceed with `M=6` and `n=8`; however, both `M` and `n` set to 6 is 
likely more than sufficient for an optimized analysis.

##### Takeways on paramater optimization

There are additional parameters beyond `M` and `n`. While there might be 
scenarios in which further optimization is required, it is easy to get lost in
parameter space, endlessly permuting values to obtain dimishing returns of 
just a handful of retained loci. We should never forget that paramter 
optimization, while important, is not the end goal of the analysis. The role 
of parameter optimization is to obtain aproximately appropiate parameters to 
then implement in the real analysis, from which to get our biological results.

#### Generating a *de novo* catalog

Now that we have optmized our parameters, we can proceed with generating a 
*de novo* catalog containing loci and genotypes from all of our individuals.

We can apply generally the same commands used during parameter optimization, 
with some key differences. First, we will use our full popmap containing our 
40 manakins. Also, we will **not** provide any missing data filters at this 
stage. The reason for this is that we want to create a static "master" 
catalog containing all the data. After inspecting this catalog, we can then 
run `populations` to apply the desired filters and exports. In fact, we could 
re-run `populations` multiple times over the same catalog, each time providing 
different filtering parameters fitting a particular downstream analysis. 
There is no need to rerun the whole *de novo* pipeline each time we need to 
new filters or generate a new export of the data.

Like before, we generate a new output directory for this *de novo* run. We 
will name this directory according to the paramters used to generate the 
catalog, in this case `M=6`.

```sh
$ mkdir denovo_M6
```

Our `denovo_map.pl` command then looks like:

```sh
$ denovo_map.pl \
      --samples ./arc-radseq-data.congen24/processed-samples/ \
      --popmap ./arc-radseq-data.congen24/info/popmap.tsv \
      --out-path denovo_M6/ \ 
      -M 6 \
      -n 6 \
      --paired \
      --rm-pcr-duplicates
```

#### Assessing the results of the *de novo* pipeline

Like we mentioned previously, the *de novo* pipeline runs several 
programs to generate loci within and between individuals. The generated 
log files allow us to track the data across the several stages of the 
pipeline, and it is good practice to check some of these summaries 
to better understand any issues that may arise.

##### Assessing per-individual assembly

The first stage we can do this is after running `ustacks`. At this stage,
`ustacks` clusters the reads of each individual into sets of loci (using the
`M` parameter we optimize before). The software will report the number of 
loci assembled at this stage, their depth of coverage, and the proportion 
of the total reads that were successfully incorportated into these initial 
set of loci. We can see these for each individual by scrolling through the 
`denovo_map.log` file, but this can be difficult when we have a large 
number of samples. Instead, we can extract a summary using the very useful 
`stacks-dist-extract` utility.

```sh
# Removing comment lines for readability...
$ stacks-dist-extract denovo_M6/denovo_map.log cov_per_sample | grep -v '^#'
  sample     loci assembled  depth of cov  max cov  number reads incorporated  % reads incorporated
  SS_02_071  3074            20.85         87       62947                      92.9
  SS_02_081  3019            17.77         67       52551                      91.8
  SS_02_082  3057            19.38         86       58052                      92.2
  SS_02_085  3092            24.28         109      73630                      92.6
  SS_02_090  3151            24.62         115      76104                      92.4
  SO_03_468  3113            24.40         96       74596                      92.9
  SO_03_469  3119            31.18         126      95644                      93.9
  SO_03_471  3108            25.59         99       77986                      92.7
  SO_03_477  3164            28.43         118      88303                      93.1
  SO_03_482  3165            23.92         92       74233                      92.8
  FC_04_510  3134            27.09         129      83572                      94.1
  FC_04_511  3174            25.62         95       79946                      93.6
  FC_04_514  3133            27.18         118      83664                      93.3
  ...
```

At this stage, we can see that the average depth of coverage of our samples 
is larger than 20x. This is great! Also, we see that >90% of reads are being 
incorporated for each sample, which is also a good result. Observing samples 
with a low proportion of reads incorporated might reflect problems with the 
data that the user should explore further (e.g., high repeat content, 
contamination, etc.)

##### Assessing the final catalog

Once the data of all individuals has been combined, the final catalog is generated 
by the `gstacks` program. It is an important place to assess the data, since 
it incorporates the final coverage after the consensus locus has been assembled 
after the removal of PCR duplicates. At this stage, we can also take a look at 
the final genotypes that we will pass to the final stages of the analysis.

We can start by taking a quick look at the `gstacks.log` file. Here, we are 
taking a quick look in the command line, by extracting a few lines around the 
line beggining with "`Genotyped`" using the program `grep`.

```sh
$ cat denovo_M6/gstacks.log | grep -A 3 -B 3 '^Genotyped' 
  Removed 24891 unpaired (forward) reads (0.9%); kept 2858492 read pairs in 4242 loci.
  Removed 1097150 read pairs whose insert length had already been seen in the same sample as putative PCR duplicates (38.4%); kept 1761342 read pairs.
  
  Genotyped 4242 loci:
    effective per-sample coverage: mean=14.9x, stdev=3.0x, min=9.9x, max=22.2x
    mean number of sites per locus: 610.7
    a consistent phasing was found for 33178 of out 34045 (97.5%) diploid loci needing phasing
```

From these few lines alone, we can can obtain some very useful information 
about our catalog. We can see that we assembled and genotyped 4,242 loci across 
all individuals. Note here that this number might not represent the real data 
in the genome, as it includes both the "real" loci that are present in the 
population, in addition to junk loci coming from repeats and error in each 
sample. This is expected, since we haven't applied any filters to this data.

We can also obtain information about our final depth of coverage. We see a mean 
coverage of 14.9x. This is good, but lower than the >20x than we saw in `ustacks`.
Why is that? This is because `gstacks` is removing PCR duplicate reads. Some of 
those initial reads making the 20x coverage were duplicates (in fact, 38.4% of 
them, as we see in the log) and were removed, leaving us with a remaining 
non-redundant depth of coverage of ~15x. PCR duplicates are a complex subject, 
but generally, we care about the non-redundant coverage than the proportion 
of duplicates themselves. We can obtain good results if we have good (>10x) 
non-redundant coverage even if the proportion of duplicates is relatively high.
For more information on PCR duplicates see 
[Rochette et al. 2023](https://doi.org/10.1111/1755-0998.13800).

Additionally, we can also look at the phasing rate (97.5%), which is indicative 
of how well data at heterozygous sites is resolved into two haplotypes at a 
locus. A higher phasing rate means that the two haplotypes are getting resolved 
properly. Seeing low phasing rates might indicate larger problems with, e.g., 
repeats and paralog sequences, both affecting how variant sites are resolved 
into two discrete haplotypes. Note also, that sites that are unable to be phased 
are filtered at this stage, so they should not affect downstream analyses.

Last, we can also take a look at the number of sites per locus (610.7 bp). When 
using paired-end data, `gstacks` will assembled the paired reads into a longer 
contig. The lenght and coverage distribution of this contig depends on the library 
configuration; however this metric can give us information about the effects of 
the molecular protocol (e.g., the size selection), and give us an idea of the 
distance at which we can observe variant sites at each locus. For example, here 
our average locus is 611 bp length, using 2x150bp reads. This means that any site 
past the 150 read length comes from the assembled paired end contig, and that we 
can potentially generate microphaplotypes around 600 bp length.

The `gstacks.log` file provides a summary of the data across all the indivuals; 
however, we can obtain more detailed info using the handy `stacks-dist-extract` 
utility script, this time on the `gstacks.log.distribs` file. This distributions 
file contains additional breakdowns of the data per each sample.

For example, we can look at the per-individual coverage and PCR duplicate rate 
using the following command:

```sh
# Removing headers and filtering some columns for readibility...
$ stacks-dist-extract denovo_M6/gstacks.log.distribs effective_coverages_per_sample | \
    grep -v '^#' | cut -f 1,2,5,8
  sample     n_loci  mean_cov_ns  pcr_dupl_rate
  CG_10_067  3080    15.181       0.382
  CG_10_069  3029    16.924       0.384
  CG_10_070  3017    13.495       0.378
  CG_10_093  3101    14.946       0.383
  CG_10_094  2994    12.997       0.379
  FC_04_510  3027    16.714       0.391
  FC_04_511  3058    15.989       0.387
  FC_04_514  3025    16.909       0.385
  FC_04_550  3077    15.518       0.388
  ...
```

Here, for each sample, we see the number of loci processed (`n_loci`), the adjusted 
average coverage per-locus (`mean_cov_ns`), and the proportion of PCR duplicate 
reads (`pcr_dupl_rate`). From these results, we can indentify outlier samples, 
e.g., those abnormally low coverage or high PCR duplicate rate that we might want 
to remove from the analysis at this stage.

The `gstacks.log.distribs` does contain additional useful information, e.g., 
about phasing. We will revisit this file 
[later](#evaluating-a-reference-based-catalog) in this protocol.

Similarly, we will look into filtering the catalog, generating exports, an 
exploring the log files of `populations` in a later 
[section](#filtering-acatalog-using-populations).

### Generating a reference-based catalog

When a reference genome of the species of interest (or a relatively close 
one) is available, the best approach is likely generating reference-based 
catalog. Instead of clusting the reads based on similarity to define the 
loci as done *de novo*, the reference-based approach uses the location of 
the aligned reads in the genome to define the loci. The reads of a locus 
are still assembled and the indiviuals are genotyped using the same key 
algorithms as a *de novo* analysis, but instead the final coordinate 
system used is that of the reference genome provided. Note that *Stacks* 
does not use the information from the reference sequence when considering 
the alleles seen in the data. This is used to prevent issues when aligning 
to a genome from a distant population or from a closely related species. Instead of a reference vs alternative alleles on the strict sense, *Stacks* 
is only using the information seen *in the sequenced data* to define what 
the alleles at a site are. Alleles are instead labelled as major and minor,
according to their frequency in the metapopulation.

The reference-based approach starts by aligning the processed readed, e.g., 
those generated by `process_radtags` to the reference genome using standard 
short-read aligners and processesing (see the [alignments](#alignments)
section for an example). The alignment must be sorted, but filtering 
the alignments is not necessary. *Stacks* will automatically filter reads 
that are unaligned and poorly aligned.

Similar to the `denovo_map.pl` wrapper script for the *de novo* pipeline, 
the reference-based pipeline also has a wrapper script used to automate 
the process: `ref_map.pl`. However, unlike the *de novo* pipeline, the 
generation of a catalog using a reference genome requires just a single 
*Stacks* program (`gstacks`). The `ref_map.pl` wrapper only runs `gstacks` 
followed by `populations`. In this case, we will just run `gstacks` to 
generate the catalog, and we can then applied the filters separately 
using `populations`.

#### Running `gstacks`

To run `gstacks` we can provide the software with a path to a directory 
contained all of our processed alignment files. In this case, we will use 
the pre-generated alignments located in `arc-radseq-data.congen24/alignments/`.
We will also specify a populations map file containing all the individuals 
that we wish to add to the catalog (`arc-radseq-data.congen24/info/popmap.tsv`).
Also, we will specify an output directory. In this case, we will create a 
new directory to store this new reference-based catalog.

```sh
$ mkdir ref_catalog
```

Lastly, there are several options for paired-end data that we could specify.
Note that these depend on the specific RAD protocol used, since not all 
options are compatible with all library configurations. In this case, since 
we have paired-end, single-digest library, we can specify that we want to 
remove PCR duplicate reads. This option is recommended whenever possible, 
but be careful if its not compatible with the library and protocol used 
(e.g., single-end sequencing, ddRAD, or GBS).

From these, we can construct our `gstacks` command:

```sh
$ gstacks \
      -I ./arc-radseq-data.congen24/alignments/ \
      -M ./arc-radseq-data.congen24/info/popmap.tsv \
      -O ./ref_catalog/ \
      --rm-pcr-duplicates 
```

Here, we are generating a reference-based catalog using only a subset 
of the whole data. In the next section, we can evaluate the generation 
of a catalog using data for the whole genome.

#### Evaluating a reference-based catalog

The `gstacks` log and distrib files for a catalog generated from the full 
dataset of these 40 manakins is available in 
`stacks-radseq/arc-radseq-data.congen24/stacks-logs/gstacks_ref/`. 
It was generated using the same `gstacks` command as before, except that 
the input alignments have not been subset to just one chromosome-level 
scaffold.

The `gstacks.log` file contains some additional information when running 
the analysis in reference mode. We can obtain some general summaries 
about the quality of our alignments.

```sh
$ cat stacks-radseq/arc-radseq-data.congen24/stacks-logs/gstacks_ref/gstacks.log | \
    grep -A 10 'BAM records'
   Read 240520637 BAM records:
     kept 205461235 primary alignments (86.0%), of which 99898265 reverse reads
     skipped 16636389 primary alignments with insufficient mapping qualities (7.0%)
     skipped 9947117 excessively soft-clipped primary alignments (4.2%)
     skipped 6957823 unmapped reads (2.9%)
     skipped some suboptimal (secondary/supplementary) alignment records

     Per-sample stats (details in 'gstacks.log.distribs'):
       read 6013015.9 records/sample (3625648-9382301)
       kept 84.2%-86.3% of these
```

For example, we see that out of 240 million alignment records, `gstacks` 
kept ~205 million (86%) as high quality, primary alignments. This is good, 
as we would ideally expect to find >80% of alignments kept at this stage; 
however, note that this is highly dependent on the quality and evolutionary 
distance of the reference genome used. A highly fragment and a distant 
reference genome would both impact the general proportions of records kept 
at this stage. Important also to note the major factor resulting in 
alignments being discarded,as it might be helpful during troubleshooting.

As done previously for the *de novo* catalog, we can use the information 
available on on `gstacks.log` to take a quick look at the overall summary 
of the data:

```sh
$ cat stacks-radseq/arc-radseq-data.congen24/stacks-logs/gstacks_ref/gstacks.log | \
    grep -A 3 -B 3 '^Genotyped' 
  Removed 9433175 unpaired (forward) reads (8.9%); 
    kept 96129795 read pairs in 132017 loci.
  Removed 41477723 read pairs whose insert length had already been seen in the 
    same sample as putative PCR duplicates (43.1%); kept 54652072 read pairs.

  Genotyped 132017 loci:
    effective per-sample coverage: mean=15.2x, stdev=3.2x, min=9.6x, max=22.9x
    mean number of sites per locus: 663.4
    a consistent phasing was found for 1265417 of out 
      1304473 (97.0%) diploid loci needing phasing
```

As we saw before, important things to look at are the average coverage (15.2x 
in this case) and the proportion of PCR duplicate reads (43.1%). Here, while 
we see over 40% PCR duplicates, the resulting non-redundant coverage should be 
sufficient to confidently call genotypes on most sites. Additional information 
on phasing and locus length is also available, as we saw 
[previously](#assessing-the-final-catalog).

Similarly, we can obtain additional per-sample information from the 
`gstacks.log.distribs` file, which we can extract using the handy 
`stacks.dist.extract` utility program. Here, we extract a per-sample 
breakdown of the depth of coverage and PCR duplicate rate.

```sh
# Removing headers and filtering some columns for readibility...
$ stacks-dist-extract gstacks.log.distribs effective_coverages_per_sample | \
    grep -v '^#' | cut -f 1,2,5,8
  sample     n_loci  mean_cov_ns  pcr_dupl_rate
  CG_10_067  94340   15.211       0.431
  CG_10_069  93509   17.590       0.435
  CG_10_070  94046   13.699       0.424
  CG_10_093  94236   14.988       0.430
  CG_10_094  93795   13.174       0.424
  FC_04_510  92295   17.485       0.436
  FC_04_511  92103   16.282       0.435
  FC_04_514  92466   17.858       0.435
  FC_04_550  92321   16.103       0.436
  FC_04_553  92160   15.626       0.435
  ...
```

Similar to the log, generating a reference-based catalog also adds 
per-sample information about the alignments to the `gstacks.logs.distribs`
file. We can also extract this using `stacks-dist-extract`:

```sh
# Removing some columns for readibility...
$ stacks-dist-extract gstacks.log.distribs bam_stats_per_sample | \
    cut -f 1,2,3,4,6,8
  sample     records  primary_kept  kept_frac  primary_disc_mapq  unmapped
  CG_10_067  6104763  5200722       0.852      435756             173531
  CG_10_069  7023636  6030505       0.859      479803             195674
  CG_10_070  5346987  4613356       0.863      352197             142061
  CG_10_093  6015970  5125528       0.852      422268             176579
  CG_10_094  5123883  4422487       0.863      339292             132274
  FC_04_510  6939233  5955713       0.858      475659             189808
  FC_04_511  6522193  5535588       0.849      465661             198273
  FC_04_514  7121135  6096359       0.856      479805             214871
  FC_04_550  6461796  5480373       0.848      479604             187477
  FC_04_553  6279269  5319959       0.847      458873             193620
  ...
```

Here, we can see for each sample the number of initial alignment records
(`records`), and the number (`primary_kept`) and proportion (`kept_frac`) 
of kept alignments. We can also obtain a per-sample breakdown of the 
different discard classes, such as alignments discarded from low mapping 
quality (`primary_disc_mapq`) or unmapped reads (`unmapped`).

Other information, for example regarding phasing, is also available. This 
might be useful to identify if any particular samples are displaying lower 
phasing rates than the rest of the library. This could be the product of, 
for example, contamination during sample preparation. At early stages of the library, mixing the DNA and/or barcode between two individuals could result 
with the data from this new sample as having more than two haplotypes, 
resulting phasing failure.

```sh
$ stacks-dist-extract gstacks.log.distribs phasing_rates_per_sample
  sample     n_gts  n_multisnp_hets  n_phased  misphasing_rate
  CG_10_067  91723  38313            37264     0.027
  CG_10_069  91499  38878            37703     0.030
  CG_10_070  91481  37093            36266     0.022
  CG_10_093  91695  37726            36774     0.025
  CG_10_094  91240  36617            35875     0.020
  FC_04_510  90805  31809            30705     0.035
  FC_04_511  90528  30720            29668     0.034
  FC_04_514  91118  32067            30955     0.035
  FC_04_550  90595  31477            30354     0.036
  FC_04_553  90733  30616            29579     0.034
```

#### Takeaways from the reference-based catalog

Overall, a good application of these distribution is to try and obtain a 
sense of the general changes in the data as it goes through the pipeline.
For example, if the samples have low depth of coverage, it might be a 
good idea to double check the alignment records to inspect if the low 
coverage is due a low number of reads or by a low proportion of reads 
mapped to the genome. Each case indiciates a different issues with the 
data and different experimental approaches to solve said issue. Here, a
low read count might indicate that the samples were undersequenced, and 
that additional sequencing might fix any issues with the coverage. In 
contrast, low proportion of alignments reads might indicate a problem 
with the reference and shows that, e.g., reads should be mapped to a 
different reference or that *de novo* approach might be better altogether.

### Filtering a catalog using `populations`

After generating the catalog, it can then be filtered to remove loci 
and/or variant sites with low quality, high missing data, or absent 
from specific populations. Here, we are filtering our reference-based 
catalog, the one we just generated in

```sh
  ~/stacks-radseq/ref_catalog
```

However, the approach is the same for catalogs generated 
*de novo*.

Note that like parameter optimization, filtering the data is also a 
complex parameter space. There is no single catch-all filters for all 
datasets. It all depends on the properties of the data and downstream 
application. Below you will see just two examples, showing two very 
generic filtering schemes for this dataset. Your data is likely to use 
a different combination of parameters. You will hear more about filtering 
and processing genotypes later in the workshop as well. 

#### Setting up our environment

First, we want to make an output directory to store all runs of 
`populations`. This will help to keep things organized as we set up analysis
using different filtering parameters.

```sh
$ mkdir filter-catalog
$ cd filter-catalog/
```

The `populations` program is designed to be modular. The idea is that the 
catalog is generated once for a set of samples. After that, we can then 
filter this catalog several times, applying the specific parameters for the 
corresponding analysis downstream and using parameters that correspond to the 
properties observed in the data. This is why we didn't apply any missing data 
filters when first generating the catalogs in the previous steps.

#### General filters

A good practice after generating the catalog, it is a good idea to apply a set 
of broad or baseline level filters. We might not use these for exporting data 
downstream, but they are helpful for observing the general quality of the data.

Let's apply some baseline levels filters. First, let's create a new output 
directory. I like naming output directories based on the parameters I am using 
in the specific run. The specific parameters will make sense in the next step.

```sh
$ mkdir populations_R50_p1_mac3
```

Let's now run `populations` on our catalog. We want to specify the paths to our 
input catalog, our output directory, and to our population map describing our 
40 manakins. Next, we want to filter data to include only loci present in 50% of
samples overall, in a single population, and with a minimum allele count of 3.

```sh
$ populations \
      --in-path ../ref_catalog/ \
      --out-path ./populations_R50_p1_mac3/ \
      --popmap ../arc-radseq-data.congen24/info/popmap.tsv \
      --min-populations 1 \
      --min-samples-overall 0.5 \
      --min-mac 3
```

Moving to full log runs in:
`~/stacks-radseq/arc-radseq-data.congen24/stacks-logs/populations_base`

```sh
$ cat populations.log | grep -B 1 'Kept'
  Removed 36681 loci that did not pass sample/population constraints 
    from 132017 loci.
  Kept 95336 loci, composed of 74654314 sites; 25156688 of those sites 
    were filtered, 689153 variant sites remained.
```

We find 95,336 loci, and 689,153 variant sites (SNPs).

#### Strict filters

```sh
$ mkdir populations_r80_p8_maf05_dp5
```

```sh
$ populations \
      --in-path ../ref_catalog/ \
      --out-path ./populations_r80_p8_maf05_dp5/ \
      --popmap ../arc-radseq-data.congen24/info/popmap.tsv \
      --min-populations 8 \
      --min-samples-per-pop 0.8 \
      --min-maf 0.05 \
      --min-gt-depth 5 \
      --ordered-export \
      --vcf
```

Moving to full log runs in:
`~/stacks-radseq/arc-radseq-data.congen24/stacks-logs/populations_strict`

```sh
$ cat populations.log | grep -B 1 'Kept'
  Removed 44254 loci that did not pass sample/population constraints 
    from 132017 loci.
  Kept 87763 loci, composed of 69491349 sites; 36665715 of those sites 
    were filtered, 199049 variant sites remained.
```

With more strict filters, 87,783 loci and 199,049 variant sites (SNPs).

This highlights some more detailed per-sample distributions, in this case 
the proportion of variant sites per sample. This distribution shows missing 
data at a SNP level. We can obtain this from the `populations.log.distribs` 
file using the `stacks-dist-extract` tool:

```sh
$ stacks-dist-extract --pretty populations.log.distribs variant_sites_per_sample | \
    grep -v '^#'
  sample     n_sites  present_sites  missing_sites  frequency_missing
  SS_02_071  199049   194259         4790           0.0241
  SS_02_081  199049   187171         11878          0.0597
  SS_02_082  199049   192134         6915           0.0347
  SS_02_085  199049   196782         2267           0.0114
  SS_02_090  199049   196763         2286           0.0115
  SO_03_468  199049   196045         3004           0.0151
  SO_03_469  199049   197677         1372           0.0069
  SO_03_471  199049   196645         2404           0.0121
  SO_03_477  199049   197425         1624           0.0082
  SO_03_482  199049   195030         4019           0.0202
  ...
```

The `frequency_missing` column, highlights the proportion of missing 
data in each sample and it is a good approach to remove low-quality 
samples, also known as the "bad apples" as described by
[Cerca et al. 2021](https://doi.org/10.1111/2041-210X.13562).

Other useful breakdowns as well, such as the distribution in the number 
of SNPs per locus:

```sh
$ stacks-dist-extract --pretty populations.log.distribs snps_per_loc_postfilters | \
    grep -v '^#' | head -n12
  n_snps  n_loci
  0       38009
  1       10607
  2       9517
  3       7492
  4       5808
  5       4429
  6       3238
  7       2422
  8       1766
  9       1278
  10      896
  ...
```

#### Takeaways

* `populations` was designed to be modular
  * Catalog is generated once, and filters are independently applied 
    later
* Each downstream application has different requirements
  * Independently run `populations` with the required filters for 
    each application
* There’s no exact or catch-all filtering parameters to use – it always 
  depends on the data and downstream application


## Acknowledgements

Thanks to Kira M. Long for permissions to use the manakin 
RADseq data, and for reviewing and editing this document.

Thanks to Alex Bangs for reviewing this document.

Thanks to Julian M. Catchen for the development of *Stacks* 
and general discussion about the software.

## Authors

**Angel G. Rivera-Colón**  
Institute of Ecology and Evolution  
University of Oregon  
Eugene, OR, USA  
<ariverac@uoregon.edu>  
<https://github.com/arcolon14>
