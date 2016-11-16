# Nextflow pipeline for ABRA (Assembly Based ReAligner)

Apply [ABRA](https://github.com/mozack/abra) to realign next generation sequencing data using localized assembly in a set of BAM files. After ABRA, the mate information is fixed using [`samtools fixmate`](http://www.htslib.org/doc/samtools.html) and BAM files are sorted and indexed using [sambamba](http://lomereiter.github.io/sambamba/).

This scripts takes a set of [BAM files](https://samtools.github.io/hts-specs/) (called `*.bam`) grouped folders as an input. There are two modes:
- When using matched tumor/normal pairs, the two samples of each pair are realigned together (see https://github.com/mozack/abra#somatic--mode). In this case the user has to provide as an input the folders containing tumor (`--tumor_bam_folder`) and normal BAM files (`--normal_bam_folder`) (it can be the same unique folder). The tumor bam file format must be (`sample` `suffix_tumor` `.bam`) with `suffix_tumor` as `_T` by default and customizable in input (`--suffix_tumor`). (e.g. `sample1_T.bam`). The normal bam file format must be (`sample` `suffix_normal` `.bam`) with `suffix_normal` as `_N` by default and customizable in input (`--suffix_normal`). (e.g. `sample1_N.bam`).
- When using only normal (or only tumor) samples, each bam is treated independently. In this case the user has to provide a single folder containing all BAM files (`bam_folder`).

In all cases BAI indexes have to be present in the same location than their BAM mates and called *.bam.bai`.

## How to install

1. Install [java](https://java.com/download/) JRE if you don't already have it.

2. Install [nextflow](http://www.nextflow.io/).

	```bash
	curl -fsSL get.nextflow.io | bash
	```
	And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
	```bash
	sudo mv nextflow /usr/local/bin
	```
	
3. Install and put in your PATH: [java](https://www.java.com/), [bedtools](http://bedtools.readthedocs.io/en/latest/), [bwa](http://bio-bwa.sourceforge.net), [sambamba](http://lomereiter.github.io/sambamba/), [samtools](http://www.htslib.org/) and download ABRA jar. Alternatively (recommended), you can simply use the docker image provided (see below).

## How to run

Simply use example:
```bash
nextflow run iarcbioinfo/abra-nf --bam_folder BAM/ --bed target.bed --ref ref.fasta --read_length 100 --abra_path /path/to/abra.jar
```

By default, BAM files produced are output in the same folder as the input folder with the `abra_sorted_fixmate.bam` suffix. One can also specify the output folder by adding the optional argument `--out_folder BAM_ABRA` to the above command line for example.

You can print the help by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/abra-nf --help
```

Instead of installing all tools in step 3 above, we recommend to use the docker image we provide containing them by simply adding `-with-docker`:
```bash
nextflow run iarcbioinfo/abra-nf -with-docker ...
```

Installing [docker](https://www.docker.com) is very system specific (but quite easy in most cases), follow  [docker documentation](https://docs.docker.com/installation/). Also follow the optional configuration step called `Create a Docker group` in their documentation.

## Detailed instructions

The exact same pipeline can be run on your computer or on a HPC cluster, by adding a [nextflow configuration file](http://www.nextflow.io/docs/latest/config.html) to choose an appropriate [executor](http://www.nextflow.io/docs/latest/executor.html). For example to work on a cluster using [SGE scheduler](https://en.wikipedia.org/wiki/Oracle_Grid_Engine), simply add a file named `nextflow.config` in the current directory (or `~/.nextflow/config` to make global changes) containing:  
```java
process.executor = 'sge'
```

Other popular schedulers such as LSF, SLURM, PBS, TORQUE etc. are also compatible. See the nextflow documentation [here](http://www.nextflow.io/docs/latest/executor.html) for more details. Also have a look at the [other parameters for the executors](http://www.nextflow.io/docs/latest/config.html#scope-executor), in particular `queueSize` that defines the number of tasks the executor will handle in a parallel manner.  
