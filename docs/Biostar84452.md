# Biostar84452

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

remove clipped bases from a BAM file


## Usage

```
Usage: biostar84452 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    -t, --tag
      tag to flag samrecord as processed
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * clip



## See also in Biostars

 * [https://www.biostars.org/p/84452](https://www.biostars.org/p/84452)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar84452
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar84452.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar84452.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar84452** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Cited In

  * Suppl. material of: "The Mobile Element Locator Tool (MELT): population-scale mobile element discovery and biology". [http://genome.cshlp.org/content/27/11/1916.short](http://genome.cshlp.org/content/27/11/1916.short) Gardner et al. Genome Res. 2017. 27: 1916-1929 
  * Molly M McDonough, Lillian D Parker, Nancy Rotzel McInerney, Michael G Campana, Jesús E Maldonado; Performance of commonly requested destructive museum samples for mammalian genomic studies, Journal of Mammalogy, , gyy080, https://doi.org/10.1093/jmammal/gyy080

## Example

```bash
$  java -jar dist/biostar84452.jar samtools-0.1.18/examples/toy.sam > out.sam

@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.biostar.Biostar84452	VN:b5ebf67dd2926d8a6afadb4d1e36a4959508057f	CL:samtools-0.1.18/examples/toy.sam
(...)
r002	0	ref	9	0	2I6M1P1I1P1I4M2I	*	0	0	AAAGATAAGGGATAAA	*
(...)


$ grep r002 samtools-0.1.18/examples/toy.sam
r002	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*

```
## See also

* https://twitter.com/EugenomeUK/status/938031803612491776

