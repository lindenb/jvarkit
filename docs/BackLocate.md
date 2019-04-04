# BackLocate

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Mapping a mutation on a protein back to the genome.


## Usage

```
Usage: backlocate [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -k, --kg
      UCSC knownGene File/URL. The knowGene format is a compact alternative to 
      GFF/GTF because one transcript is described using only one line.	Beware 
      chromosome names are formatted the same as your REFERENCE. A typical 
      KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
      Default: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz
    -x, -X, --kgxref
      UCSC kgXRef URI. Must have at least 5 columns. $1 is knowGene-Id $5  is 
      protein identifier.
      Default: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz
    -o, --out
      Output file. Optional . Default: stdout
    -p, --printSeq
      print mRNA & protein sequences
      Default: false
  * -R, --reference
      The parameter is the path to an Indexed fasta Reference file. This fasta 
      file must be indexed with samtools faidx and with picard 
      CreateSequenceDictionary. The parameter can also be a 'key' (matching 
      the regular expression `[A-Za-z][A-Za-z0-9_\\-]*`) in a catalog file. A 
      'catalog' file is a java property file ( 
      https://docs.oracle.com/javase/tutorial/essential/environment/properties.html 
      ) where the values are the path to the fasta file.  Catalogs are 
      searched in that order : `${PWD}/fasta-ref.properties`, 
      `${HOME}/.fasta-ref.properties`, `/etc/jvarkit/fasta-ref.properties`.  
      If the key or the path are not defined by the user, they will be 
      searched in that order 1) the java property 
      -Djvarkit.fasta.reference=pathTofastaOrCatalogKey . 2) the linux 
      environement variable $FASTA_REFERENCE=pathTofastaOrCatalogKey 3) The 
      catalogs. 
      Default: <<Default Fasta Reference Supplier>>
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * prediction
 * protein



## See also in Biostars

 * [https://www.biostars.org/p/15992](https://www.biostars.org/p/15992)
 * [https://www.biostars.org/p/116366](https://www.biostars.org/p/116366)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew backlocate
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/backlocate/BackLocate.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/backlocate/BackLocate.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **backlocate** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

mutation P->M at 1090 in NOTCH2

```
$  echo -e "NOTCH2\tP1090M" | java -jar dist/backlocate.jar -R hg19.fa
(...)
[WARNING/BackLocate] 2014-11-05 12:03:08 "The reference doesn't contain chromosome chr17_ctg5_hap1"
[WARNING/BackLocate] 2014-11-05 12:03:15 "The reference doesn't contain chromosome chr4_ctg9_hap1"
[WARNING/BackLocate] 2014-11-05 12:03:16 "The reference doesn't contain chromosome chr6_apd_hap1"
[WARNING/BackLocate] 2014-11-05 12:03:16 "The reference doesn't contain chromosome chr6_cox_hap2"
[WARNING/BackLocate] 2014-11-05 12:03:16 "The reference doesn't contain chromosome chr6_dbb_hap3"
(...)
[INFO/BackLocate] 2014-11-05 12:03:18 "genes:78963"
[INFO/BackLocate] 2014-11-05 12:03:18 "loading http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz"
[INFO/BackLocate] 2014-11-05 12:03:24 "kgxref:28493"
(...)
```

```
#User.Gene	AA1	petide.pos.1	AA2	knownGene.name	knownGene.strandknownGene.AA	index0.in.rna	codon	base.in.rna	chromosome	index0.in.genomic	exon
##uc001eik.3
NOTCH2	P	1090	M	uc001eik.3	NEGATIVE	P	3267	CCA	C	chr1	120480548	Exon 20
NOTCH2	P	1090	M	uc001eik.3	NEGATIVE	P	3268	CCA	C	chr1	120480547	Exon 20
NOTCH2	P	1090	M	uc001eik.3	NEGATIVE	P	3269	CCA	A	chr1	120480546	Exon 20
##uc001eil.3
NOTCH2	P	1090	M	uc001eil.3	NEGATIVE	P	3267	CCA	C	chr1	120480548	Exon 20
NOTCH2	P	1090	M	uc001eil.3	NEGATIVE	P	3268	CCA	C	chr1	120480547	Exon 20
NOTCH2	P	1090	M	uc001eil.3	NEGATIVE	P	3269	CCA	A	chr1	120480546	Exon 20
```

```
$ echo -e "NOTCH2\tPro1090M\tInteresting" | java -jar dist/backlocate.jar -R /path/to/human_g1k_v37.fasta | grep -v "##" | java -jar dist/prettytable.jar 

+------------+-----+--------------+-----+----------------+------------------+--------------+---------------+------------+----------------------+-------------+------------+-------------------+---------+-----------------+
| #User.Gene | AA1 | petide.pos.1 | AA2 | knownGene.name | knownGene.strand | knownGene.AA | index0.in.rna | wild.codon | potential.var.codons | base.in.rna | chromosome | index0.in.genomic | exon    | extra.user.data |
+------------+-----+--------------+-----+----------------+------------------+--------------+---------------+------------+----------------------+-------------+------------+-------------------+---------+-----------------+
| NOTCH2     | Pro | 1090         | Met | uc001eik.3     | -                | P            | 3267          | CCA        | .                    | C           | 1          | 120480548         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eik.3     | -                | P            | 3268          | CCA        | .                    | C           | 1          | 120480547         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eik.3     | -                | P            | 3269          | CCA        | .                    | A           | 1          | 120480546         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eil.3     | -                | P            | 3267          | CCA        | .                    | C           | 1          | 120480548         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eil.3     | -                | P            | 3268          | CCA        | .                    | C           | 1          | 120480547         | Exon 20 | Interesting     |
| NOTCH2     | Pro | 1090         | Met | uc001eil.3     | -                | P            | 3269          | CCA        | .                    | A           | 1          | 120480546         | Exon 20 | Interesting     |
+------------+-----+--------------+-----+----------------+------------------+--------------+---------------+------------+----------------------+-------------+------------+-------------------+---------+-----------------+
```


## See also

 * http://plindenbaum.blogspot.fr/2011/03/mapping-mutation-on-protein-to-genome.html
 * https://github.com/lindenb/jvarkit/issues/14
 * https://github.com/lindenb/jvarkit/issues/13


## History

 * 2019: add extra user data
 * 2017: Moved to jcommander
 * 2014: Moved to jvarkit
 * Nov 2014 : removed all the dependencies to SQL and DAS; use a local indexed genome
 * Aug 2015 : Added a new column "potention var codon" (as https://twitter.com/_ramrs/status/631123002005061633 ) , renamed "codon" to "wild codon"

## Cited in

backlocate was cited in:

 * CRISPR-STOP: gene silencing through base-editing-induced nonsense mutations. 2017 Nat Meth. [http://dx.doi.org/10.1038/nmeth.4327](http://dx.doi.org/10.1038/nmeth.4327).
 * "Differential 3′ Processing of Specific Transcripts Expands Regulatory and Protein Diversity Across Neuronal Cell Types" Sasa Jereb, Hun-Way Hwang, Eric Van Otterloo, Eve-Ellen Govek, John J Fak, Yuan Yuan, Mary E Hatten, Robert B Darnell BioRxiv [https://www.biorxiv.org/content/biorxiv/early/2018/01/10/245886.full.pdf](https://www.biorxiv.org/content/biorxiv/early/2018/01/10/245886.full.pdf)
 
