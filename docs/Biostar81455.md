# Biostar81455

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Defining precisely the exonic genomic context based on a position .


## Usage

```
Usage: biostar81455 [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -KG, --knownGene
      UCSC knownGene File/URL. The knowGene format is a compact alternative to 
      GFF/GTF because one transcript is described using only one line.	Beware 
      chromosome names are formatted the same as your REFERENCE. A typical 
      KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
      Default: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -1
      The coordinate are one-based. The default is zero based.
      Default: false

```


## Keywords

 * bed
 * gene
 * knownGene
 * ucsc



## See also in Biostars

 * [https://www.biostars.org/p/81455](https://www.biostars.org/p/81455)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar81455
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar81455.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar81455.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar81455Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar81455Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar81455** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

tab delimited file. 2 columns: CHROM and POS

## Example

```bash
echo -e "chr22\t41258261\nchr22\t52000000\nchr22\t0" |\
	java   dist/biostar81455.jar 

chr22	41258261	uc003azg.2	41253084	41258785	POSITIVE	Exon 2	41257621	41258785	0
chr22	41258261	uc011aox.2	41253084	41305239	POSITIVE	Exon 1	41253084	41253249	-5012
chr22	41258261	uc003azi.3	41253084	41328823	POSITIVE	Exon 1	41253084	41253249	-5012
chr22	41258261	uc003azj.3	41255553	41258130	NEGATIVE	Exon 1	41255553	41258130	-131
chr22	41258261	uc010gyh.1	41258260	41282519	POSITIVE	Exon 1	41258260	41258683	0
chr22	41258261	uc011aoy.1	41258260	41363888	POSITIVE	Exon 1	41258260	41258683	0
chr22	52000000	uc011asd.2	51195513	51227614	POSITIVE	Exon 4	51227177	51227614	-772386
chr22	52000000	uc003bni.3	51195513	51238065	POSITIVE	Exon 4	51237082	51238065	-761935
chr22	52000000	uc011ase.1	51205919	51220775	NEGATIVE	Exon 1	51220615	51220775	-779225
chr22	52000000	uc003bnl.1	51205919	51222087	NEGATIVE	Exon 1	51221928	51222087	-777913
chr22	52000000	uc003bns.3	51222156	51238065	POSITIVE	Exon 3	51237082	51238065	-761935
chr22	52000000	uc003bnq.1	51222224	51227600	POSITIVE	Exon 4	51227322	51227600	-772400
chr22	52000000	uc003bnr.1	51222224	51227781	POSITIVE	Exon 4	51227319	51227781	-772219
chr22	52000000	uc010hbj.3	51222224	51238065	POSITIVE	Exon 3	51237082	51238065	-761935
chr22	0	uc002zks.4	16150259	16193004	NEGATIVE	Exon 8	16150259	16151821	16150259
chr22	0	uc002zkt.3	16162065	16172265	POSITIVE	Exon 1	16162065	16162388	16162065
chr22	0	uc002zku.3	16179617	16181004	NEGATIVE	Exon 1	16179617	16181004	16179617
chr22	0	uc002zkv.3	16187164	16193004	NEGATIVE	Exon 5	16187164	16187302	16187164	
```

