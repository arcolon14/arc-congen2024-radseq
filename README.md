# RADseq data analysis with *Stacks* 2 hands-on lecture

Documentation and examples for the RADseq analysis with *Stacks* 2 hands-on 
exercise at ConGen 2024, on Aug 27, 2024.

(C) Angel G. Rivera-Colón (<ariverac@uoregon.edu>)

## Associated readings

The RADseq data used for this exercise comes from the following publication:

> Long, KM, Rivera-Colón, AG, Bennett, KFP, *et al*. (2024) **Ongoing introgression 
> of a secondary sexual plumage trait in a stable avian hybrid zone**. *Evolution*, 
> 2024, qpae076, [DOI: 10.1093/evolut/qpae076](https://doi.org/10.1093/evolut/qpae076)

The data is used here with the permission of the authors.

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

For information on the download, installation, and documentation of the software, 
please visit the *Stacks* [website](https://catchenlab.life.illinois.edu/stacks/).

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

<!---Check this!--->
```sh
$ cp /data/instructor_materials/Angel_Rivera-Colon/2024/arc-radseq-data.congen24.tar.gz .
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
  CG_10_067.bam  PR_09_096.bam  RO_05_206.bam  SO_03_468.bam
  CG_10_069.bam  PR_09_102.bam  RO_05_221.bam  SO_03_469.bam
  CG_10_070.bam  PR_09_111.bam  RO_05_222.bam  SO_03_471.bam
  CG_10_093.bam  PR_09_112.bam  RO_05_223.bam  SO_03_477.bam
  CG_10_094.bam  PR_09_115.bam  RO_05_254.bam  SO_03_482.bam
  FC_04_510.bam  QP_06_190.bam  RU_08_203.bam  SS_02_071.bam
  FC_04_511.bam  QP_06_191.bam  RU_08_204.bam  SS_02_081.bam
  FC_04_514.bam  QP_06_195.bam  RU_08_205.bam  SS_02_082.bam
  FC_04_550.bam  QP_06_226.bam  RU_08_245.bam  SS_02_085.bam
  FC_04_553.bam  QP_06_228.bam  RU_08_249.bam  SS_02_090.bam

  catalog:
  catalog.calls     catalog.fa.gz  gstacks.log.distribs
  catalog.chrs.tsv  gstacks.log

  info:
  barcodes.tsv  popmap.tsv  popmap_catalog.tsv

  processed-samples:
  CG_10_067.1.fq.gz  PR_09_112.2.fq.gz  RU_08_205.1.fq.gz
  CG_10_067.2.fq.gz  PR_09_115.1.fq.gz  RU_08_205.2.fq.gz
  CG_10_069.1.fq.gz  PR_09_115.2.fq.gz  RU_08_245.1.fq.gz
  CG_10_069.2.fq.gz  QP_06_190.1.fq.gz  RU_08_245.2.fq.gz
  CG_10_070.1.fq.gz  QP_06_190.2.fq.gz  RU_08_249.1.fq.gz
  CG_10_070.2.fq.gz  QP_06_191.1.fq.gz  RU_08_249.2.fq.gz
  CG_10_093.1.fq.gz  QP_06_191.2.fq.gz  SO_03_468.1.fq.gz
  CG_10_093.2.fq.gz  QP_06_195.1.fq.gz  SO_03_468.2.fq.gz
  CG_10_094.1.fq.gz  QP_06_195.2.fq.gz  SO_03_469.1.fq.gz
  CG_10_094.2.fq.gz  QP_06_226.1.fq.gz  SO_03_469.2.fq.gz
  FC_04_510.1.fq.gz  QP_06_226.2.fq.gz  SO_03_471.1.fq.gz
  FC_04_510.2.fq.gz  QP_06_228.1.fq.gz  SO_03_471.2.fq.gz
  FC_04_511.1.fq.gz  QP_06_228.2.fq.gz  SO_03_477.1.fq.gz
  FC_04_511.2.fq.gz  RO_05_206.1.fq.gz  SO_03_477.2.fq.gz
  FC_04_514.1.fq.gz  RO_05_206.2.fq.gz  SO_03_482.1.fq.gz
  FC_04_514.2.fq.gz  RO_05_221.1.fq.gz  SO_03_482.2.fq.gz
  FC_04_550.1.fq.gz  RO_05_221.2.fq.gz  SS_02_071.1.fq.gz
  FC_04_550.2.fq.gz  RO_05_222.1.fq.gz  SS_02_071.2.fq.gz
  FC_04_553.1.fq.gz  RO_05_222.2.fq.gz  SS_02_081.1.fq.gz
  FC_04_553.2.fq.gz  RO_05_223.1.fq.gz  SS_02_081.2.fq.gz
  PR_09_096.1.fq.gz  RO_05_223.2.fq.gz  SS_02_082.1.fq.gz
  PR_09_096.2.fq.gz  RO_05_254.1.fq.gz  SS_02_082.2.fq.gz
  PR_09_102.1.fq.gz  RO_05_254.2.fq.gz  SS_02_085.1.fq.gz
  PR_09_102.2.fq.gz  RU_08_203.1.fq.gz  SS_02_085.2.fq.gz
  PR_09_111.1.fq.gz  RU_08_203.2.fq.gz  SS_02_090.1.fq.gz
  PR_09_111.2.fq.gz  RU_08_204.1.fq.gz  SS_02_090.2.fq.gz
  PR_09_112.1.fq.gz  RU_08_204.2.fq.gz

  raw-reads:
  MAVIRAD2_NoIndex_L002_R1_001.fastq.gz
  MAVIRAD2_NoIndex_L002_R2_001.fastq.gz

  stacks-logs:
  denovo       populations_base    process_radtags
  gstacks_ref  populations_strict
```

## Breakdown of the available data

This data directory contains several subdirectories describing this RADseq dataset at 
different stages of the analysis.

### `raw-reads`

The `raw-reads` directory contains the raw sequencing files (in `FASTQ` format) from 
this library. As described in [Long et al. 2024](https://doi.org/10.1093/evolut/qpae076), 
total of 130 *Manacus* samples were processed using a single-digest RADseq library 
([Baird et al. 2008](https://doi.org/10.1371/journal.pone.0003376); 
[Etter et al. 2011](https://doi.org/10.1371/journal.pone.0018561)) using the *SbfI* 
restriction enzyme and paired-end sequenced on an Illumina NovaSeq6000, generating 
2x150bp reads.

```sh
$ ls raw-reads/
  MAVIRAD2_NoIndex_L002_R1_001.fastq.gz   MAVIRAD2_NoIndex_L002_R2_001.fastq.gz
```

**Note:** For the sake of time, the files here only contain a fraction (~1 million 
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
(in `FASTQ` format) for 40 *Manacus* individuals. They are the demultiplexed, cleaned 
reads generated by running `process_radtags` on the raw data. In the excersices below, 
we are just using 40 out of the total 130 samples.

See `stacks-logs/process_radtags/process_radtags.MAVI2.log` for additional details 
about how these reads were processed.

**Note:** For the sake of time, the files here only contain a subset of the total 
reads for each individual, those corresponding to a single chromosome-level scaffold.

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

The `alignments` directory contains the aligned reads (in `bam` format) from 40 *Manacus* 
individuals. The data for these 40 samples was run previously through `process_radtags` 
and aligned to the *Manacus vitellinus* RefSeq assembly (NCBI accession 
[GCF_001715985.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001715985.3/)) using 
*BWA* `mem` ([Li 2013](https://doi.org/10.48550/arXiv.1303.3997)) and processed with the 
`view` and `sort` commands from *SAMtools* 
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

**NOTE:** For the sake of time, the `bam` files provided here have been filtered to include 
data from only one chromosome-level scaffold.

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

**NOTE:** For the sake of time, this file only contains barcodes for 40 out of the total 
130 *Manacus* samples sequenced in the library.

#### The popmap files

The `info` directory also contains two population map (i.e., "popmap") files. Each popmap is 
a tab-delimited file describing the group assignment of each sequenced individual.

The first popmap file (`popmap.tsv`), contains the information for 40 *Manacus* individuals. 
Each is assigned to a specific sampling site along a transect in Panama, as described by 
[Long et al. 2024](https://doi.org/10.1093/evolut/qpae076). In total, we have 40 samples 
across 8 total populations, each representing a different sampling site.

```sh
$ cat info/popmap.tsv
  SS_02_071<tab>020SS
  SS_02_081     020SS
  SO_03_468     030SO
  SO_03_469     030SO
  PR_09_096     090PR
  PR_09_102     090PR
  CG_10_067     100CG
  CG_10_069     100CG
  ...
```

As described in the *Stacks* [manual](https://catchenlab.life.illinois.edu/stacks/manual/#popmap), 
the first column contains the individual's ID, while the second column contains the ID of 
its respective group or population. In this example, individual `SS_02_081` is assigned to 
population `020SS` (shortened ID for transect site #2, San San Drury in 
[Long et al. 2024](https://doi.org/10.1093/evolut/qpae076)), while individual `PR_09_096` 
belongs to population `090PR` (site #9, Palma Real).

The second popmap, `popmap_catalog.tsv`, only contains a subset of 15 individuals, all 
assigned to the sample population (labelled `opt`).

```sh
$ cat info/popmap_catalog.tsv
  SS_02_071<tab>opt
  SS_02_081	    opt
  PR_09_096	    opt
  PR_09_102	    opt
  CG_10_067	    opt
  CG_10_069	    opt
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
  cov_per_sample.tsv  denovo_map.log

  stacks-logs/gstacks_ref:
  gstacks.log  gstacks.log.distribs

  stacks-logs/populations_base:
  populations.log  populations.log.distribs

  stacks-logs/populations_strict:
  populations.log  populations.log.distribs

  stacks-logs/process_radtags:
  process_radtags.MAVI2.log
```

For the sake of time, in this exercise we are running the *Stacks* pipeline on reduced 
datasets. While fast, the logs from these reduced runs might not always be the most 
informative. Thus, this `stacks-logs` directory contains logs from the real, full dataset 
and will be used to explore the analysis in more detail. 

<!---TODO: edit from here down--->

## Creating a catalog with `gstacks`

In the main `stacks-radseq`, let's create a directory to store a new *Stacks* catalog. 
This will be the output directory for `gstacks`.

```sh
$ mkdir stacks-catalog
```

Let's then move into this directory

```sh
$ cd stacks-catalog
```

Once there, we want to run the *Stacks* `gstacks` program to create a new catalog of 
RADseq loci and variant sites generated from the aligned reads of our 60 *Manacus* 
samples, as specified with the popmap file. Since our data is aligned to a genome, we 
will be running the software in reference mode (by providing the path to the `bam` 
files). Since these samples were prepared in a single-digest RADseq library and 
sequenced used paired-end reads, we all able to also remove PCR duplicates when 
processing our new RAD loci.

Here's an example of the `gstacks` command:

```sh
$ gstacks \
    -I ~/stacks-radseq/arc-radseq-data.congen23/alignments/ \
    -O . \
    -M ~/stacks-radseq/arc-radseq-data.congen23/info/popmap.tsv \
    --threads 4 \
    --rm-pcr-duplicates
```

### Inspecting the catalog

Once `gstacks` finishes running, we can inspect the coverage and PCR duplicate summary 
statistics from the `gstacks.log` file.

**NOTE:** Due to limited time, we are running `gstacks` on a set of subsampled 
alignments. The examples below show results for a run containing data for the whole 
genome. A copy of this larger catalog can be found in 
`~/stacks-radseq/arc-radseq-data.congen23/stacks_data/gstacks`.

```sh
$ cat gstacks.log | grep -B 3 -A 5 '^Genotyped'
  Removed 9864375 unpaired (forward) reads (6.3%); kept 145744842 read pairs in 130231 loci.
  Removed 102123643 read pairs whose insert length had already been seen in the
    same sample as putative PCR duplicates (70.1%); kept 43621199 read pairs.
  
  Genotyped 130231 loci:
    effective per-sample coverage: mean=9.0x, stdev=6.9x, min=1.4x, max=24.3x
    mean number of sites per locus: 655.6
    a consistent phasing was found for 1163253 of out 1204439 (96.6%) diploid loci needing phasing
  
  gstacks is done.
```

This catalog contains 130 thousand assembled loci (average length 656 bp). The average 
non-redundant coverage of 9x after removing 70% PCR duplicates. 97% of all loci were 
phased into haplotypes.

### Per-individual catalog statistics

The values above are a summary of the whole catalog. Looking at diagnostic distributions 
at an individual level might provide additional information regarding the properties 
of the catalog.

#### Alignment statistics

```sh
$ stacks-dist-extract gstacks.log.distribs bam_stats_per_sample
  sample          records   primary_kept  kept_frac  primary_kept_read2  primary_disc_mapq  primary_disc_sclip  unmapped  secondary  supplementary
  CG_10_LIB1_104  3233278   2849075       0.881      1416707             203100             74916               93885     0          12302
  CG_10_LIB1_107  2124239   1887469       0.889      937391              117824             48689               60770     0          9487
  CG_10_LIB1_110  3338987   2937839       0.880      1460636             200089             85624               99854     0          15581
  CG_10_LIB1_133  17473133  15522677      0.888      7719372             1006107            410941              463345    0          70063
  CG_10_LIB1_135  7550756   6723388       0.890      3342453             427202             174647              196429    0          29090
  CG_10_LIB2_067  6104763   5200722       0.852      2528449             435756             253311              173531    0          41443
  CG_10_LIB2_069  7023636   6030505       0.859      2928792             479803             277672              195674    0          39982
  CG_10_LIB2_070  5346987   4613356       0.863      2241003             352197             208934              142061    0          30439
  CG_10_LIB2_093  6015970   5125528       0.852      2489385             422268             252187              176579    0          39408
  ...
```

Reformatted here for readability:

| `sample` | `records` | `primary_kept` | `kept_frac` | `primary_kept_read2` | `primary_disc_mapq` | `primary_disc_sclip` | `unmapped` |
| ---------------- | ---------- | ---------- | ----- | --------- | --------- | ------- | ------- |
| `CG_10_LIB1_104` | 3,233,278  | 2,849,075  | 0.881 | 1,416,707 | 203,100   | 74,916  | 93,885  |
| `CG_10_LIB1_107` | 2,124,239  | 1,887,469  | 0.889 | 937,391   | 117,824   | 48,689  | 60,770  |
| `CG_10_LIB1_110` | 3,338,987  | 2,937,839  | 0.880 | 1,460,636 | 200,089   | 85,624  | 99,854  |
| `CG_10_LIB1_133` | 17,473,133 | 15,522,677 | 0.888 | 7,719,372 | 1,006,107 | 410,941 | 463,345 |
| `CG_10_LIB1_135` | 7,550,756  | 6,723,388  | 0.890 | 3,342,453 | 427,202   | 174,647 | 196,429 |
| `CG_10_LIB2_067` | 6,104,763  | 5,200,722  | 0.852 | 2,528,449 | 435,756   | 253,311 | 173,531 |
| `CG_10_LIB2_069` | 7,023,636  | 6,030,505  | 0.859 | 2,928,792 | 479,803   | 277,672 | 195,674 |
| `CG_10_LIB2_070` | 5,346,987  | 4,613,356  | 0.863 | 2,241,003 | 352,197   | 208,934 | 142,061 |
| `CG_10_LIB2_093` | 6,015,970  | 5,125,528  | 0.852 | 2,489,385 | 422,268   | 252,187 | 176,579 |

#### Per-sample non-redundant coverage and PCR duplicates

```sh
$ stacks-dist-extract gstacks.log.distribs effective_coverages_per_sample | grep -v '^#'
  sample          n_loci  n_used_fw_reads  mean_cov  mean_cov_ns  n_unpaired_reads  n_pcr_dupl_pairs  pcr_dupl_rate
  CG_10_LIB1_104  40247   75435            1.874     1.924        53289             1303644           0.945
  CG_10_LIB1_107  31601   49348            1.562     1.591        37140             863590            0.946
  CG_10_LIB1_110  41963   77687            1.851     1.902        57391             1342125           0.945
  CG_10_LIB1_133  78905   401437           5.088     5.465        293328            7108540           0.947
  CG_10_LIB1_135  59591   176098           2.955     3.103        127093            3077744           0.946
  CG_10_LIB2_067  94340   1385738          14.689    16.182       238750            1047785           0.431
  CG_10_LIB2_069  93509   1597878          17.088    18.665       275124            1228711           0.435
  CG_10_LIB2_070  94046   1245928          13.248    14.561       208150            918275            0.424
  CG_10_LIB2_093  94236   1365631          14.492    15.916       238788            1031724           0.430
  ...
```

Reformatted here for readability:

| `sample` | `n_loci` | `n_used_fw_reads` | `mean_cov` | `mean_cov_ns` | `n_unpaired_reads` | `n_pcr_dupl_pairs` | `pcr_dupl_rate` |
|----------------- | ------ | --------- | ------ | ------ | ------- | --------- | ----- |
| `CG_10_LIB1_104` | 40,247 | 75,435    | 1.874  | 1.924  | 53,289  | 1,303,644 | 0.945 |
| `CG_10_LIB1_107` | 31,601 | 49,348    | 1.562  | 1.591  | 37,140  | 863,590   | 0.946 |
| `CG_10_LIB1_110` | 41,963 | 77,687    | 1.851  | 1.902  | 57,391  | 1,342,125 | 0.945 |
| `CG_10_LIB1_133` | 78,905 | 401,437   | 5.088  | 5.465  | 293,328 | 7,108,540 | 0.947 |
| `CG_10_LIB1_135` | 59,591 | 176,098   | 2.955  | 3.103  | 127,093 | 3,077,744 | 0.946 |
| `CG_10_LIB2_067` | 94,340 | 1,385,738 | 14.689 | 16.182 | 238,750 | 1,047,785 | 0.431 |
| `CG_10_LIB2_069` | 93,509 | 1,597,878 | 17.088 | 18.665 | 275,124 | 1,228,711 | 0.435 |
| `CG_10_LIB2_070` | 94,046 | 1,245,928 | 13.248 | 14.561 | 208,150 | 918,275   | 0.424 |
| `CG_10_LIB2_093` | 94,236 | 1,365,631 | 14.492 | 15.916 | 238,788 | 1,031,724 | 0.430 |

#### Phasing

```sh
$ stacks-dist-extract gstacks.log.distribs phasing_rates_per_sample
  sample          n_gts  n_multisnp_hets  n_phased  misphasing_rate  n_phased_2ndpass
  CG_10_LIB1_104  39781  4025             3962      0.016            32
  CG_10_LIB1_107  31283  2270             2245      0.011            21
  CG_10_LIB1_110  41468  4166             4111      0.013            19
  CG_10_LIB1_133  77663  20838            20200     0.031            240
  CG_10_LIB1_135  58785  10114            9901      0.021            78
  CG_10_LIB2_067  91696  38802            37449     0.035            138
  CG_10_LIB2_069  91468  39171            37704     0.037            170
  CG_10_LIB2_070  91450  37588            36469     0.030            143
  CG_10_LIB2_093  91671  38126            36867     0.033            129
  ...
```

Reformatted here for readability:

| `sample` | `n_gts` | `n_multisnp_hets` | `n_phased` | `misphasing_rate` |`n_phased_2ndpass` |
| ---------------- | ------ | ------ | ------ | ----- | --- |
| `CG_10_LIB1_104` | 39,781 | 4,025  | 3,962  | 0.016 | 32  |
| `CG_10_LIB1_107` | 31,283 | 2,270  | 2,245  | 0.011 | 21  |
| `CG_10_LIB1_110` | 41,468 | 4,166  | 4,111  | 0.013 | 19  |
| `CG_10_LIB1_133` | 77,663 | 20,838 | 20,200 | 0.031 | 240 |
| `CG_10_LIB1_135` | 58,785 | 10,114 | 9,901  | 0.021 | 78  |
| `CG_10_LIB2_067` | 91,696 | 38,802 | 37,449 | 0.035 | 138 |
| `CG_10_LIB2_069` | 91,468 | 39,171 | 37,704 | 0.037 | 170 |
| `CG_10_LIB2_070` | 91,450 | 37,588 | 36,469 | 0.030 | 143 |
| `CG_10_LIB2_093` | 91,671 | 38,126 | 36,867 | 0.033 | 129 |

## Filtering the catalog and exporting genotypes

Create general directory

```sh
$ mkdir ~/stacks-radseq/filter-catalog
$ cd ~/stacks-radseq/filter-catalog
```

### General `populations` run

Get r80 loci present in the at least three of the six *Manacus* populations and 
observe alleles present at least three times (minimum count of 3, i.e., present 
in at least 2 samples). Run with the popmap containing all 60 samples.

Create the directory of this run. The name of the directory 
(`populations.p3.r80.mac3`) describes the specific filters applied to the data.

```sh
$ mkdir populations.p3.r80.mac3
$ cd populations.p3.r80.mac3 
```

Run the *Stacks* `populations` module

```sh
$ populations \
    --in-path ~/stacks-radseq/stacks-catalog \
    --out-path . \
    --popmap ~/stacks-radseq/arc-radseq-data.congen23/info/popmap.tsv \
    --threads 4 \
    --min-populations 3 \
    --min-samples-per-pop 0.80 \
    --min-mac 3
```

### Check the outputs of `populations`

**Note:** for the sake of time, we ran `populations` on the reduced catalog 
we prepared earlier. We will use a larger run (available in 
`~/stacks-radseq/arc-radseq-data.congen23/stacks_data/populations/populations.p3.r80.mac3`) 
to explore the filtering of the catalog.

Go to the large `populations` run:

```sh
$ cd ~/stacks-radseq/arc-radseq-data.congen23/stacks_data/populations/populations.p3.r80.mac3
```

List the contents of this directory:

```sh
$ ls populations/populations.p3.r80.mac3/
  populations.haplotypes.tsv  populations.log.distribs
  populations.haps.vcf        populations.snps.vcf
  populations.hapstats.tsv    populations.sumstats.tsv
  populations.log             populations.sumstats_summary.tsv
```

The `sumstats` and `hapstats` files contain the summary statistics assigned 
per-population for each SNPs and haplotype, respectively. The SNPs and 
haplotypes are also exported in VCF format.

Inspect the `populations.log` file to obtain the number of loci and variant 
sites retained after filtering:

```sh
$ cat populations.log | grep 'Kept'
  Kept 72066 loci, composed of 56612370 sites; 40540446 of those sites were filtered, 155674 variant sites remained.
```

After applying filters, this run kept 72 thousand loci, containing 156 
thousand variant sites.

Inspect the `populations.log.distribs` to obtain additional diagnostic 
distributions and per-sample missing data statistics.

#### Samples per-locus

```sh
$ stacks-dist-extract populations.log.distribs samples_per_loc_postfilters | grep -v '^#'
  n_samples   n_loci
  24          1528
  25          1250
  26          1176
  27          631
  28          408
  ...
  56          5404
  57          6720
  58          7612
  59          8920
  60          8308
```

#### Missing loci per-sample

```sh
$ stacks-dist-extract populations.log.distribs loci_per_sample | grep -v '^#'
  sample          n_loci  present_loci  missing_loci  frequency_missing
  SS_02_LIB1_076  72066   43799         28267         0.3922
  SS_02_LIB1_086  72066   55188         16878         0.2342
  SS_02_LIB1_124  72066   56081         15985         0.2218
  SS_02_LIB1_153  72066   63975         8091          0.1123
  SS_02_LIB1_179  72066   57717         14349         0.1991
  SS_02_LIB2_071  72066   64933         7133          0.0990
  SS_02_LIB2_081  72066   64924         7142          0.0991
  SS_02_LIB2_082  72066   64904         7162          0.0994
  SS_02_LIB2_085  72066   64935         7131          0.0990
  ...
```

Reformatted here for readability:

| `sample`         | `n_loci` | `present_loci` | `missing_loci` | `frequency_missing` |
| ---------------- | -------- | -------------- | -------------- | ------------------- |
| `SS_02_LIB1_076` | 72,066   | 43,799         | 28,267         | 0.3922              |
| `SS_02_LIB1_086` | 72,066   | 55,188         | 16,878         | 0.2342              |
| `SS_02_LIB1_124` | 72,066   | 56,081         | 15,985         | 0.2218              |
| `SS_02_LIB1_153` | 72,066   | 63,975         | 8,091          | 0.1123              |
| `SS_02_LIB1_179` | 72,066   | 57,717         | 14,349         | 0.1991              |
| `SS_02_LIB2_071` | 72,066   | 64,933         | 7,133          | 0.0990              |
| `SS_02_LIB2_081` | 72,066   | 64,924         | 7,142          | 0.0991              |
| `SS_02_LIB2_082` | 72,066   | 64,904         | 7,162          | 0.0994              |
| `SS_02_LIB2_085` | 72,066   | 64,935         | 7,131          | 0.0990              |

#### Missing variant sites per-sample

```sh
$ stacks-dist-extract populations.log.distribs variant_sites_per_sample | grep -v '^#'
  sample          n_sites  present_sites  missing_sites  frequency_missing
  SS_02_LIB1_076  102800   30421          72379          0.7041
  SS_02_LIB1_086  102800   82220          20580          0.2002
  SS_02_LIB1_124  102800   84021          18779          0.1827
  SS_02_LIB1_153  102800   79498          23302          0.2267
  SS_02_LIB1_179  102800   67139          35661          0.3469
  SS_02_LIB2_071  102800   102249         551            0.0054
  SS_02_LIB2_081  102800   101728         1072           0.0104
  SS_02_LIB2_082  102800   102231         569            0.0055
  SS_02_LIB2_085  102800   102380         420            0.0041
  ...
```

Reformatted here for readability

| `sample`         | `n_sites` | `present_sites` | `missing_sites` | `frequency_missing` |
| ---------------- | --------- | --------------- | --------------- | ------------------- |
| `SS_02_LIB1_076` | 102,800   |  30,421         | 72,379          | 0.7041              |
| `SS_02_LIB1_086` | 102,800   |  82,220         | 20,580          | 0.2002              |
| `SS_02_LIB1_124` | 102,800   |  84,021         | 18,779          | 0.1827              |
| `SS_02_LIB1_153` | 102,800   |  79,498         | 23,302          | 0.2267              |
| `SS_02_LIB1_179` | 102,800   |  67,139         | 35,661          | 0.3469              |
| `SS_02_LIB2_071` | 102,800   |  102,249        | 551             | 0.0054              |
| `SS_02_LIB2_081` | 102,800   |  101,728        | 1,072           | 0.0104              |
| `SS_02_LIB2_082` | 102,800   |  102,231        | 569             | 0.0055              |
| `SS_02_LIB2_085` | 102,800   |  102,380        | 420             | 0.0041              |


## Authors

**Angel G. Rivera-Colon**  
Institute of Ecology and Evolution  
University of Oregon  
Eugene, OR, USA  
<ariverac@uoregon.edu>  
<https://github.com/arcolon14>
