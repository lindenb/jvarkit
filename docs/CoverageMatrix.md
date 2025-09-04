# CoverageMatrix

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

generate a VCF file from bam coverage


## Usage

```
Usage: java -jar dist/coveragematrix.jar  [options] Files
Usage: coveragematrix [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --bin, --bin-size
      Bin size. A distance specified as a positive integer.Commas are removed. 
      The following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 1000
    --chrom, --contig
      Restrict to that contig.
    --black, --exclude
      Optional. BED Tabix indexed black-listed region
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      min mapping quality
      Default: 1
    --max-depth
      ignore position if depth > 'x'
      Default: 500
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --treshold
      DUP if 1.5-x<=depth<=1.5+x . HET_DEL if 0.5-x<=depth<=0.5+x HOM_DEL if 
      0.0-x<=depth<=0.0+x. 
      Default: 0.05
    --version
      print version and exit

```


## Keywords

 * cnv
 * bam
 * depth
 * coverage


## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone --recurse-submodules "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew coveragematrix
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200618

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/CoverageMatrix.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/CoverageMatrix.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **coveragematrix** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
 java -jar dist/coveragematrix.jar -R src/test/resources/rotavirus_rf.fa --exclude gaps.tsv.gz src/test/resources/S*.bam
```



