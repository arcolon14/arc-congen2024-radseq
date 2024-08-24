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

The analysis described follow the general guidelines described in the *Stacks* 2 
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

**NOTE:** the directory hierarchy in the commands below referes to the 
ConGen2024 server.

Make a directory for the Stacks assignment.

```sh
$ mkdir stacks-radseq
$ cd stacks-radseq/
```

Copy the raw data from the instructors directory

```sh
$ cp ../instructor_materials/Angel_Rivera-Colon/2024/arc-radseq-data.congen24.tar.gz .
```

Uncompress this directory

<!---Check this!--->
```sh
$ tar xvf arc-radseq-data.congen24.tar.gz
```

Check the contents of the directory

```sh
$ ls *
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
total of 130 manakins (genus *Manacus*) samples were processed using a single-digest 
RADseq library ([Baird et al. 2008](https://doi.org/10.1371/journal.pone.0003376); 
[Etter et al. 2011](https://doi.org/10.1371/journal.pone.0018561)) using the *SbfI* 
restriction enzyme and paired-end sequenced on an Illumina NovaSeq6000, generating 
2x150bp reads.

```sh
$ ls raw-reads/
  MAVIRAD2_NoIndex_L002_R1_001.fastq.gz   MAVIRAD2_NoIndex_L002_R2_001.fastq.gz
```

**Note:** Due to limited time, the files here only contain a fraction (~1 million 
reads) of the files in the real library.

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

The `alignments` directory contains the aligned reads (in `bam` format) from 40 manakins
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

The file `barcodes.tsv` is a tab-delimted file containing the individual barcode sequences 
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
distribution files for several , pre-run steps of the *Stacks* pipeline.

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
4. The quality (or `PHRED`) scores.

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

Following the library configuration, we do not expect to see a barcode nor 
a restriction enzyme cutsite.

### Processing raw data

After inspecting the raw data, we can start processing it using the *Stacks*
`process_radtags` program. We can demultiplex the combined reads into the 
data for the individual samples, and we can filter the reads to remove 
those of low quality.

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
using a single restriction enzyme, *SbfI*. Other paramters can be left default, 
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
`stacks-dist-extract` program.

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
of the reads were due to not having a proper barcode sequence (4.5%); however, 
losses due to low quality and missing cutsites were minimal (0.1% and 0.5%, 
respectively).

We can also look at a per-sample breakdown of the reads. Once again, we can acess this 
information using the `stacks-dist-extract` command.

```sh
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
do four simple runs, ranging from `M=3` to `M=6`. On real dataset, testing over 
a larger range of values is recommended (e.g., 1 to 12).

Additionally, will generate a catalog using only a subset of samples, with the 
assumption that these samples are representative of the whole dataset. This is 
often recommended to provide a general speedup of the process. For this purpose, we 
have provided a secondary popmap (`info/popmap_catalog.tsv`) just containing 15 
manakin samples, all assigned to a single group labelled `opt`.

We will also generate an output directory to store the outputs of this run. We 
will name this directory according to the value of `M`, 3 in this example.

```sh
$ mkdir param_opt_M3
```

We then want to run `denovo_map.pl`, specfying the value of `M`, the smaller 
popmap, and our new output directory. As our input data, we will specify the 
processed the FASTQ files in `arc-radseq-data.congen24/processed-samples/`. These 


$ denovo_map.pl       --samples ./arc-radseq-data.congen24/processed-samples/       --popmap ./arc-radseq-data.congen24/info/popmap_catalog.tsv       --out-path ./param_opt_M4/      -M 4       -n 4       --paired       --rm-pcr-duplicates -r 0.8

$ cat param_opt_M4/populations.log | grep '^Kept'Kept 2668 loci, composed of 1715600 sites; 460105 of those sites were filtered, 17014 variant sites remained.

####

```sh
$ mkdir denovo_M4
```

Remember, we are not using the raw data generated in the previous step.

```sh
$ denovo_map.pl \
      --samples ./arc-radseq-data.congen24/processed-samples/ \
      --popmap ./arc-radseq-data.congen24/info/popmap.tsv \
      --out-path denovo_M4/ \
      -M 4 \
      -n 4 \
      --paired \
      --rm-pcr-duplicates
```

### Create a reference-based catalog

```sh
$ mkdir ref_catalog
```

```sh
$ gstacks \
      -I ./arc-radseq-data.congen24/alignments/ \
      -M ./arc-radseq-data.congen24/info/popmap.tsv \
      -O ./ref_catalog/ \
      --rm-pcr-duplicate
```

### Filter the catalog using `populations`

#### General filters

```sh
$ mkdir populations_R50_p1_mac3
```

```sh
$ populations \
      --in-path ./ref_catalog/ \
      --out-path ./populations_R50_p1_mac3/ \
      --popmap ./arc-radseq-data.congen24/info/popmap.tsv \
      --min-populations 1 \
      --min-samples-overall 0.5 \
      --min-mac 3
```

#### Strict filters

```sh
$ mkdir populations_r80_p8_maf05_dp5
```

```sh
$ populations \
      --in-path ./ref_catalog/ \
      --out-path ./populations_r80_p8_maf05_dp5/ \
      --popmap ./arc-radseq-data.congen24/info/popmap.tsv \
      --min-populations 8 \
      --min-samples-per-pop 0.8 \
      --min-maf 0.05 \
      --min-gt-depth 5 \
      --ordered-export \
      --vcf
```

<!---Delete after this--->

## Authors

**Angel G. Rivera-Colón**  
Institute of Ecology and Evolution  
University of Oregon  
Eugene, OR, USA  
<ariverac@uoregon.edu>  
<https://github.com/arcolon14>
