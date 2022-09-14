# GtfToBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert GTF/GFF3 to BED.


## Usage

```
Usage: java -jar dist/gtf2bed.jar  [options] Files
Usage: gtf2bed [options] Files
  Options:
    -c, --columns
      comma separated columns to be displayed
      Default: gtf.source,gtf.feature,gtf.score,gtf.strand,gtf.frame,ccds_id,exon_id,exon_number,exon_version,gene_biotype,gene_id,gene_name,gene_source,gene_version,havana_transcript,havana_transcript_version,protein_id,protein_version,tag,transcript_biotype,transcript_id,transcript_name,transcript_source,transcript_version,Alias,biotype,ccdsid,constitutive,description,ensembl_end_phase,ensembl_phase,exon_id,external_name,gene_id,havana_transcript,havana_version,ID,logic_name,Name,Parent,protein_id,rank,tag,transcript_id,version
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * gtf
 * gff
 * gff3
 * bed


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew gtf2bed
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20220629

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/GtfToBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/GtfToBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gtf2bed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


```
```

