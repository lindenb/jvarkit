# Gff2KnownGene

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert GFF3/GTF format to UCSC knownGene format.


## Usage

```
Usage: gff2kg [options] Files
  Options:
    -bin, --bin
      Insert  UCSC 'bin' column as the first column.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    -trid, --trid
      Transcript identifiers in the GTF/GFF (column N°3) used to identify a 
      transcript.Multiple separated by a semicolon
      Default: transcript;mRNA;snRNA;tRNA;snoRNA
    -verbose, --verbose
      Be verbose, log messages
      Default: false
    --version
      print version and exit

```


## Keywords

 * gff
 * ,gtf
 * knownGene
 * ucsc
 * convert



## See also in Biostars

 * [https://www.biostars.org/p/276099](https://www.biostars.org/p/276099)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew gff2kg
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/Gff2KnownGene.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/Gff2KnownGene.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/Gff2KnownGeneTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/Gff2KnownGeneTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gff2kg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$  curl -s "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz" |\
	gunzip -c |\
	java -jar dist/gff2kg.jar
(...)
1826	ENST00000367917.3	chr1	+	162760522	162782607	162760590	162782210	8	162760522,162762448,162766374,162767591,162769532,162774056,162775183,162782087	162760625,162762652,162766467,162767706,162769727,162774113,162775282,162782607	gene_id=ENSG00000132196.9;transcript_id=ENST00000367917.3;gene_type=protein_coding;gene_status=KNOWN;gene_name=HSD17B7;transcript_type=protein_coding;transcript_name=HSD17B7-201;protein_id=ENSP00000356894.3;havana_gene=OTTHUMG00000034420.6;	ENST00000367917.3
(...)
```

In the UCSC (not the structure of konwGene, but we can validate intervals):

```
$ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -D hg19 -e 'select * from wgEncodeGencodeBasicV19 where name="ENST00000367917.3"' | cat
bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
1826	ENST00000367917.3	chr1	+	162760522	162782607	162760590	162782210	8	162760522,162762448,162766374,162767591,162769532,162774056,162775183,162782087,	162760625,162762652,162766467,162767706,162769727,162774113,162775282,162782607,	0	HSD17B7	cmpl	cmpl	0,2,2,2,0,0,0,0,
```

### From ensembl 

```
$	wget -O -  "ftp://ftp.ensembl.org/pub/grch37/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz" |\
	gunzip -c |\
	java -jar dist/gff2kg.jar
```

## see also

  * Ensembl vs UCSC  [https://twitter.com/yokofakun/status/743751004785545218](https://twitter.com/yokofakun/status/743751004785545218)



