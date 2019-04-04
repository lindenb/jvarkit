# VcfIn

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Only prints variants that are contained/not contained into another VCF


## Usage

```
Usage: vcfin [options] Files
  Options:
    -A, --allalt
      ALL user ALT must be found in VCF-database ALT
      Default: false
  * -D, --database
      external database uri
      Default: []
    -fi, --filterin
      Do not discard variant but add this FILTER if the variant is found in 
      the database(s)
      Default: <empty string>
    -fo, --filterout
      Do not discard variant but add this FILTER if the variant is NOT found 
      in the database(s)
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --inverse
      Print variant that are not part of the VCF-database.
      Default: false
    -cp, --only-contig-pos
      Two variants are the same if they have the same CONTIG/POS. Do NOT Look 
      at REF or ALTS. Motivation: e.g  
      http://gnomad.broadinstitute.org/dbsnp/rs11361742 in gnomad VCF: 
      [19	45454285 TAA	TA,T,TAAA ] in my VCF: [19	45454285	TA	T]. Only works 
      with --tabix option
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -t, --tabix, --tribble, --indexed
      Database is indexed with tabix or tribble. I will use random access to 
      get the data. Otherwise, VCF are assumed to be sorted the same way. 
      [20180731] I will try to convert the contig names e.g: 'chr1' -> '1'
      Default: false
    --version
      print version and exit
    -M
       max number of equivalent variants found in database(s). -1 : no limit. 
      sinclusive.With --tabix mode only. .eg: '3': the user variant must be 
      found in 3 or less VCF database.
      Default: -1
    -m
      min number of equivalent variants found in database(s), inclusive. With 
      --tabix mode only. .eg: '2': the user variant must be found in at least 
      2 VCF database.
      Default: 1

```


## Keywords

 * vcf
 * compare



## See also in Biostars

 * [https://www.biostars.org/p/287815](https://www.biostars.org/p/287815)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfin
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfIn.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfIn.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfInTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfInTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfin** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


VCF files should be sorted using the same order as the sequence dictionary (see picard SortVcf).


### Example

list variants found with gatk AND samtools, keep the variants with http://www.sequenceontology.org/browser/current_release/term/SO:0001818 , remove variants found in a previous alignment (samtools or gatk)


```
#!/bin/bash

ls Samples | while read S
do
gunzip -c NEWALIGN/{S}.gatk.vcf.gz |\
        java -jar jvarkit-git/dist/vcffilterso.jar -A SO:0001818 |\
        java -jar jvarkit-git/dist/vcfin.jar -D NEWALIGN/{S}.samtools.vcf.gz |\
        java -jar jvarkit-git/dist/vcfin.jar -i -D OLDALIGN/{S}.samtools.vcf.gz |
        java -jar jvarkit-git/dist/vcfin.jar -i -D OLDALIGN/${S}.gatk.vcf.gz |
        grep -vE '^#' |
        awk -v S=${S} -F '      ' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",S,$1,$2,$3,$4,$5,$8);}' 
done
```


#### Example 2

My list of bad variants is in the file 'bad.vcf'.


```
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	11167517	.	G	A	.	.	.
1	11167760	.	C	T	.	.	.
1	11168529	.	A	G	.	.	.
1	11169789	.	A	G	.	.	.
1	11174331	.	T	C	.	.	.
1	11174715	.	T	C	.	.	.
1	11180949	.	C	T	.	.	.
```

My main vcf file is 'input.vcf'.


```
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	11166480	.	C	T	.	.	.
1	11166541	.	G	C	.	.	.
1	11166577	.	C	T	.	.	.
1	11166713	.	T	C	.	.	.
1	11167146	.	G	A	.	.	.
1	11167158	.	C	T	.	.	.
1	11167270	.	G	T	.	.	.
1	11167517	.	G	T	.	.	.
1	11167627	.	G	C	.	.	.
1	11167760	.	C	T	.	.	.
1	11167829	.	C	T	.	.	.
1	11168529	.	A	G	.	.	.
1	11168769	.	CAAA	C	.	.	.
1	11169250	.	G	A	.	.	.
1	11169420	.	G	A	.	.	.
1	11169440	.	A	G	.	.	.
1	11169585	.	A	G	.	.	.
1	11169624	.	T	C	.	.	.
```


I want to put a FILTER in the variants if they are contained in VCF.


```
java -jar dist/vcfin.jar -A -fi InMyListOfBadVariants -D jeter2.vcf jeter1.vcf
```

output:

```
##fileformat=VCFv4.2
##FILTER=<ID=InMyListOfBadVariants,Description="Variant overlapping database.">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	11166480	.	C	T	.	.	.
1	11166541	.	G	C	.	.	.
1	11166577	.	C	T	.	.	.
1	11166713	.	T	C	.	.	.
1	11167146	.	G	A	.	.	.
1	11167158	.	C	T	.	.	.
1	11167270	.	G	T	.	.	.
1	11167517	.	G	T	.	.	.
1	11167627	.	G	C	.	.	.
1	11167760	.	C	T	.	InMyListOfBadVariants	.
1	11167829	.	C	T	.	.	.
1	11168529	.	A	G	.	InMyListOfBadVariants	.
1	11168769	.	CAAA	C	.	.	.
1	11169250	.	G	A	.	.	.
1	11169420	.	G	A	.	.	.
1	11169440	.	A	G	.	.	.
1	11169585	.	A	G	.	.	.
1	11169624	.	T	C	.	.	.

```

Please note that variant 1	11167517 is not flagged because is alternate allele is not contained in 'bad.vcf'


### History

 *  2018-08-31: database must be specified with '-D'. ContigName conversion applied for tabix/tribble data.
 *  2015-02-24: rewritten. all files must be sorted: avoid to sort on disk. Support for tabix. Option -A
 *  2015-01-26: changed option '-v' to option '-i' (-v is for version)
 *  2014: Creation

