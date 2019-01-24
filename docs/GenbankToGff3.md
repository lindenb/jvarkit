# GenbankToGff3

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Experimental genbank to GFF


## Usage

```
Usage: gb2gff [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * xml
 * ncbi
 * genbank
 * convert
 * gff
 * gb


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew gb2gff
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/genbank/GenbankToGff3.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/genbank/GenbankToGff3.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gb2gff** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Warnings

 - still experimental
 - efetch doesn't always work: https://gist.github.com/lindenb/6c4f36bdf29a3108e103e3a5a0b1aff7

## Example

```
$ wget -O - -q  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AF338247.1,D38149.1&rettype=gbwithparts&retmode=xml&api_key=62b713e0cd85e6ac79699ecdfa72e85af009"  | java -jar dist/gb2gff.jar 

##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10970
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10941
##sequence-region AF338247.1 1 2032
##sequence-region D38149.1 1 1087
AF338247.1	genbank	source	1	2032	.	+	.	strain=M;db_xref=taxon:10941;note="rearranged segment 7";organism="Human rotavirus A";segment=7R;clone=M1;mol_type="genomic RNA"
AF338247.1	genbank	five_prime_UTR	1	34	.	+	.	.
AF338247.1	genbank	CDS	35	967	.	+	1	codon_start=1;product=NSP3;protein_id=AAK74117.1;transl_table=1
AF338247.1	genbank	CDS	993	1925	.	+	1	note="duplicated ORF";codon_start=1;product=NSP3;protein_id=AAK74118.1;transl_table=1
AF338247.1	genbank	three_prime_UTR	1926	2032	.	+	.	.
D38149.1	genbank	source	1	1087	.	+	.	strain=A5-16;db_xref=taxon:10970;organism="Rotavirus sp.";segment=5;mol_type="genomic RNA"
D38149.1	genbank	gene	33	185	.	+	.	gene=NSP1
D38149.1	genbank	CDS	33	185	.	+	1	codon_start=1;protein_id=BAA07347.1;gene=NSP1;transl_table=1
D38149.1	genbank	variation	141	141	.	+	.	note="bases 142-641 in D38148 deleted";gene=NSP1
```

## input 

Input is a list of *.xml genbank files (DTD/Schema: https://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd ). Otherwise assume that file is a list of filenames and unfold it.

