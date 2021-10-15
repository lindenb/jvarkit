# BamToHaplotypes

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Reconstruct SNP haplotypes from reads


## Usage

```
Usage: bam2haplotypes [options] Files
  Options:
    --buffer-size
      When we're looking for variant in a lare VCF file, load the variants in 
      an interval of 'N' bases instead of doing a random access for each 
      variant. 
      Default: 1000
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-discordant-rg
      In paired mode, ignore discordant read-groups RG-ID.
      Default: false
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --out
      Output file. Optional . Default: stdout
    --paired
      Activate Paired-end mode. Variant can be supported by the read or/and is 
      mate. Input must be sorted on query name using for example 'samtools 
      collate'. 
      Default: false
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
  * -V, --vcf
      Indexed VCf file. Only diallelic SNP will be considered.
    --version
      print version and exit

```


## Keywords

 * vcf
 * phased
 * genotypes
 * bam



## See also in Biostars

 * [https://www.biostars.org/p/9493599](https://www.biostars.org/p/9493599)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bam2haplotypes
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20211015

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/phased/BamToHaplotypes.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/phased/BamToHaplotypes.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bam2haplotypes** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



# Example
```
$ java -jar dist/bam2haplotypes.jar -V src/test/resources/rotavirus_rf.vcf.gz src/test/resources/S5.bam

#CHROM	START	END	COUNT	N-VARIANTS	(POS\tALT)+
RF03	1221	1242	11	2	1221	C	1242	C
RF03	1688	1708	5	2	1688	G	1708	T
RF04	1900	1920	4	2	1900	C	1920	A
RF06	517	543	9	2	517	C	543	G
RF06	668	695	4	2	668	G	695	T
RF08	926	992	2	2	926	C	992	G
RF09	294	317	6	2	294	T	317	A
RF10	139	175	1	2	139	T	175	G
RF10	139	175	3	2	139	T	175	C
```

in **paired** mode

```
samtools collate -O -u src/test/resources/S5.bam TMP | java -jar dist/bam2haplotypes.jar --paired -V src/test/resources/rotavirus_rf.vcf.gz

#CHROM	START	END	COUNT	N-VARIANTS	(POS\tALT)+
RF02	251	578	1	2	251	A	578	G
RF03	1221	1688	1	2	1221	C	1688	G
RF03	1221	1242	7	2	1221	C	1242	C
RF03	1221	1688	1	2	1221	C	1688	G
RF03	1688	1708	1	2	1688	G	1708	T
RF03	1708	2150	1	2	1708	T	2150	T
RF03	1221	1708	1	3	1221	C	1688	G	1708	T
RF03	1221	1688	2	3	1221	C	1242	C	1688	G
RF03	1688	2150	1	3	1688	G	1708	T	2150	T
RF03	1221	1708	2	4	1221	C	1242	C	1688	G	1708	T
RF04	887	1241	1	2	887	A	1241	T
RF04	1900	1920	4	2	1900	C	1920	A
RF05	41	499	2	2	41	T	499	A
RF05	499	879	1	2	499	A	879	C
RF05	795	1297	2	2	795	A	1297	T
RF05	879	1297	2	2	879	C	1297	T
RF06	517	543	9	2	517	C	543	G
RF06	668	695	4	2	668	G	695	T
RF07	225	684	1	2	225	C	684	G
RF07	225	684	1	2	225	C	684	T
RF08	926	992	2	2	926	C	992	G
RF09	294	317	6	2	294	T	317	A
RF10	139	175	1	2	139	T	175	G
RF10	139	175	3	2	139	T	175	C
```

