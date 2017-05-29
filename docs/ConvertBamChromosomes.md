# ConvertBamChromosomes

Convert the names of the chromosomes in a BAM file


## Usage

```
Usage: bamrenamechr [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
  * -f, --mapping, -m
      load a custom name mapping. Format (chrom-source\tchrom-dest\n)+
    -o, --out
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: TypeImpl{name='SAM', fileExtension='sam', indexExtension='null'}
    --version
      print version and exit
    -c, -convert
      What should I do when  a converstion is not found
      Default: RAISE_EXCEPTION
      Possible Values: [RAISE_EXCEPTION, SKIP, RETURN_ORIGINAL]

```


## Keywords

 * sam
 * bam
 * chromosome
 * contig


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make bamrenamechr
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ConvertBamChromosomes.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/ConvertBamChromosomes.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamrenamechr** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example

```bash
$ cat samtools-0.1.19/examples/toy.sam

@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
r002	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*
r003	0	ref	9	30	5H6M	*	0	0	AGCTAA	*
r004	0	ref	16	30	6M14N1I5M	*	0	0	ATAGCTCTCAGC	*
r003	16	ref	29	30	6H5M	*	0	0	TAGGC	*
r001	83	ref	37	30	9M	=	7	-39	CAGCGCCAT	*
x1	0	ref2	1	30	20M	*	0	0	aggttttataaaacaaataa	????????????????????
x2	0	ref2	2	30	21M	*	0	0	ggttttataaaacaaataatt	?????????????????????
x3	0	ref2	6	30	9M4I13M	*	0	0	ttataaaacAAATaattaagtctaca	??????????????????????????
x4	0	ref2	10	30	25M	*	0	0	CaaaTaattaagtctacagagcaac	?????????????????????????
x5	0	ref2	12	30	24M	*	0	0	aaTaattaagtctacagagcaact	????????????????????????
x6	0	ref2	14	30	23M	*	0	0	Taattaagtctacagagcaacta	???????????????????????

java -jar dist/bamrenamechr.jar \
   -f <(echo -e "ref\tCHROM1\nref2\tCHROM2")
   samtools-0.1.19/examples/toy.sam

@HD	VN:1.4	SO:unsorted
@SQ	SN:CHROM1	LN:45
@SQ	SN:CHROM2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.misc.ConvertBamChromosomes	VN:dfab75cb8c06e47e9989e59df62ec8f3242934c4	CL:-f /dev/fd/63 /commun/data/packages/samtools-0.1.19/examples/toy.sam
r001	163	CHROM1	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
r002	0	CHROM1	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*
r003	0	CHROM1	9	30	5H6M	*	0	0	AGCTAA	*
r004	0	CHROM1	16	30	6M14N1I5M	*	0	0	ATAGCTCTCAGC	*
r003	16	CHROM1	29	30	6H5M	*	0	0	TAGGC	*
r001	83	CHROM1	37	30	9M	=	7	-39	CAGCGCCAT	*
x1	0	CHROM2	1	30	20M	*	0	0	AGGTTTTATAAAACAAATAA	????????????????????
x2	0	CHROM2	2	30	21M	*	0	0	GGTTTTATAAAACAAATAATT	?????????????????????
x3	0	CHROM2	6	30	9M4I13M	*	0	0	TTATAAAACAAATAATTAAGTCTACA	??????????????????????????
x4	0	CHROM2	10	30	25M	*	0	0	CAAATAATTAAGTCTACAGAGCAAC	?????????????????????????
x5	0	CHROM2	12	30	24M	*	0	0	AATAATTAAGTCTACAGAGCAACT	????????????????????????
x6	0	CHROM2	14	30	23M	*	0	0	TAATTAAGTCTACAGAGCAACTA	???????????????????????

```

## See also

* https://github.com/lindenb/jvarkit/blob/master/src/main/resources/chromnames/g1kv37_to_hg19.tsv
* https://github.com/lindenb/jvarkit/blob/master/src/main/resources/chromnames/hg19_to_g1kv37.tsv
* [[VcfRenameChromosomes]]
* http://plindenbaum.blogspot.fr/2013/07/g1kv37-vs-hg19.html



