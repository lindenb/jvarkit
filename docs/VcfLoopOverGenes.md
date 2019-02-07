# VcfLoopOverGenes

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Generates a BED file of the Genes in an annotated VCF, loop over those genes and generate a VCF for each gene, then execute an optional command.


## Usage

```
Usage: vcfloopovergenes [options] Files
  Options:
    -compress, --compress
      generate VCF.gz files
      Default: false
    --contigWinLength
      [20171018] window size when splitting per contig
      Default: 1000
    --contigWinShift
      [20171018] window shift when splitting per contig
      Default: 500
    -delete, --delete
      if a command if executed with '-exec', delete the file after the 
      completion of the command.
      Default: false
    -e, -exec, --exec
      When saving the VCF to a directory. Execute the following command line. 
      The words __PREFIX__  __CONTIG__ (or __CHROM__ ) __ID__ __NAME__ 
      __SOURCE__ __VCF__ __START__ __END__ will be replaced by their values.
      Default: <empty string>
    -g, --gene, -gene, --genes
      Loop over the gene file. If not defined VCF will be scanned for SnpEff 
      annotations and output will be a BED file with the gene names and 
      provenance.If defined, I will create a VCF for each Gene.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -j, --jobs
      When -exec is specified, use <n> jobs. A value lower than 1 means use 
      all procs available.
      Default: 1
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      For gene.bed: a File or stdout. When creating a VCF this should be a zip 
      file or an existing directory.
    -p, -prefix, --prefix
      File prefix when saving the individual VCF files.
      Default: <empty string>
    -r, --region
      An interval as the following syntax : "chrom:start-end" or 
      "chrom:middle+extend"  or "chrom:start-end+extend" or 
      "chrom:start-end+extend-percent%".A program might use a Reference 
      sequence to fix the chromosome name (e.g: 1->chr1)
      Default: <empty string>
    --snpEffNoIntergenic
      [20170711] when using SNPEFF annotations ignore intergenic variants. 
      Makes things faster if you're only working with protein-things.
      Default: false
    --splitMethod
      [20170711] How to split primary vcf
      Default: Annotations
      Possible Values: [Annotations, VariantSlidingWindow, ContigSlidingWindow]
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --variantsWinCount
      [20170711] when split per count of variants, put at most this number of 
      variants in the chunk.
      Default: 1000
    --variantsWinShift
      [20170711] when split per count of variants, shift the window by this 
      number of variants.
      Default: 500
    --version
      print version and exit

```


## Keywords

 * vcf
 * gene
 * burden


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfloopovergenes
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfLoopOverGenes.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfLoopOverGenes.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfloopovergenes** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

Generate the bed file from a VCF annotated with SnpEff

```
$ java -jar dist/vcfloopovergenes.jar -p KARAKA input.vcf.gz > genes.bed 
$ head jeter.bed
13	62807	62808	KARAKA.000000002	2V3477:ENST000002607	ANN_FeatureID	1
13	11689	20735	KARAKA.000000004	AC07.1	ANN_GeneName	30
13	75803	90595	KARAKA.000000006	ENSG000000781	ANN_GeneID	284
13	44306	68961	KARAKA.000000008	ENSG00044491	ANN_GeneID	1545
(...)
```

Generate the VCFs:


```
 $ java -jar dist/vcfloopovergenes.jar -p KARAKA -g genes.bed -o tmp input.vcf.gz
 

 $ head tmp/KARAKA.manifest.txt 
KARAKA.3_KARAKA.000000001.vcf
KARAKA.3_KARAKA.000000002.vcf
KARAKA.3_KARAKA.000000003.vcf
KARAKA.3_KARAKA.000000004.vcf
(..)
 ```

 ```
$ ls tmp/*.vcf | head
tmp/KARAKA.3_KARAKA.000000001.vcf
tmp/KARAKA.3_KARAKA.000000002.vcf
tmp/KARAKA.3_KARAKA.000000003.vcf
tmp/KARAKA.3_KARAKA.000000004.vcf
(...)
```

### Example 'Exec'

```
$ java -jar dist/vcfloopovergenes.jar -p KARAKA -g genes.bed -exec "echo __ID__ __PREFIX__ __VCF__ __CONTIG__" -o tmp input.vcf.gz 
KARAKA.000000001 KARAKA. tmp/KARAKA.000000001.vcf 3
KARAKA.000000002 KARAKA. tmp/KARAKA.000000002.vcf 3
KARAKA.000000003 KARAKA. tmp/KARAKA.000000003.vcf 3
KARAKA.000000004 KARAKA. tmp/KARAKA.000000004.vcf 3
KARAKA.000000005 KARAKA. tmp/KARAKA.000000005.vcf 3
KARAKA.000000006 KARAKA. tmp/KARAKA.000000006.vcf 3
KARAKA.000000007 KARAKA. tmp/KARAKA.000000007.vcf 3
KARAKA.000000008 KARAKA. tmp/KARAKA.000000008.vcf 3
KARAKA.000000009 KARAKA. tmp/KARAKA.000000009.vcf 3
KARAKA.000000010 KARAKA. tmp/KARAKA.000000010.vcf 3
KARAKA.000000011 KARAKA. tmp/KARAKA.000000011.vcf 3
KARAKA.000000012 KARAKA. tmp/KARAKA.000000012.vcf 3
```


### Example 'Exec'

execute the following Makefile 'count.mk' for each gene:

```make
all:${MYVCF}
	echo -n "Number of variants in $< : " && grep -vE '^#' $< | wc -l	
```

```
java -jar dist/vcfloopovergenes.jar \
	-p MATILD -gene genes.bed input.vcf.gz \
	-exec 'make -f count.mk MYVCF=__VCF__' -o tmp
```



