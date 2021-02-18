# Biostar404363

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

introduce artificial mutation SNV in bam


## Usage

```
Usage: biostar404363 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    --disable-nm
      Disable NM change. By default the value of NM is increasing to each 
      change. 
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-clip
      Ignore clipped section of the reads.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
  * -p, --position, --vcf
      VCF File containing the positions to change. if INFO/AF(allele 
      frequency) field is present, variant is inserted if rand()<= AF.
    -R, --reference
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard 
      CreateSequenceDictionary 
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * variant



## See also in Biostars

 * [https://www.biostars.org/p/404363](https://www.biostars.org/p/404363)
 * [https://www.biostars.org/p/416897](https://www.biostars.org/p/416897)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar404363
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20191023

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar404363.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar404363.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar404363Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar404363Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar404363** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## History

major bug detected after https://github.com/lindenb/jvarkit/issues/178 . Some reads with INDEL might have been ignored.

##Example

```
# look at the original bam at position 14:79838836
$ samtools tview -d T -p 14:79838836 src/test/resources/HG02260.transloc.chr9.14.bam | cut -c 1-20
     79838841  79838
NNNNNNNNNNNNNNNNNNNN
CATAGGAAAACTAAAGGCAA
cataggaaaactaaaggcaa
cataggaaaaCTAAAGGCAA
cataggaaaactaaaggcaa
CATAGGAAAACTAAAGGCAA
cataggaaaactaaaggcaa
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
           taaaggcaa

# look at the VCF , we want to change 14:79838836 with 'A' with probability = 0.5
$ cat jeter.vcf
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency among genotypes, for each ALT allele, in the same order as listed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
14	79838836	.	N	A	.	.	AF=0.5

# apply biostar404363
$  java -jar dist/biostar404363.jar -o jeter.bam -p jeter.vcf src/test/resources/HG02260.transloc.chr9.14.bam

# index the sequence
$ samtools index jeter.bam

# view the result (half of the bases were replaced because AF=0.5 in the VCF)
$ samtools tview -d T -p 14:79838836 jeter.bam | cut -c 1-20
     79838841  79838
NNNNNNNNNNNNNNNNNNNN
MATAGGAAAACTAAAGGCAA
cataggaaaactaaaggcaa
aataggaaaaCTAAAGGCAA
cataggaaaactaaaggcaa
CATAGGAAAACTAAAGGCAA
aataggaaaactaaaggcaa
AATAGGAAAACTAAAGGCAA
AATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
AATAGGAAAACTAAAGGCAA
           taaaggcaa
```

