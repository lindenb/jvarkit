# Sam2Tsv

Prints the SAM alignments as a TAB delimited file.


## Usage

```
Usage: sam2tsv [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -A, --printAlignments
      Print Alignments
      Default: false
    -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * table
 * tsv



## See also in Biostars

 * [https://www.biostars.org/p/157232](https://www.biostars.org/p/157232)
 * [https://www.biostars.org/p/59647](https://www.biostars.org/p/59647)
 * [https://www.biostars.org/p/253828](https://www.biostars.org/p/253828)
 * [https://www.biostars.org/p/264875](https://www.biostars.org/p/264875)
 * [https://www.biostars.org/p/277493](https://www.biostars.org/p/277493)


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
$ make sam2tsv
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2Tsv.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2Tsv.java)


<details>
<summary>Git History</summary>

```
Wed Oct 11 22:38:59 2017 +0200 ; add biostars 277493 ; https://github.com/lindenb/jvarkit/commit/0649c3b450c251a774a6e6085cbfc826ea150485
Wed Oct 4 17:05:37 2017 +0200 ; reading google scholar, answers to reviewers ; https://github.com/lindenb/jvarkit/commit/871a481468fbd1877f02bc171cf080c5e1d3190f
Thu Jul 27 16:58:18 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/a8aaf2d7df89f44442b36ee1120ee4dd5c1e36e6
Thu Jun 29 20:41:53 2017 +0200 ; start fix https://github.com/lindenb/jvarkit/issues/81 ; https://github.com/lindenb/jvarkit/commit/3758963956c9ceec249feeb4a076c98314756c14
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Sun May 21 17:11:09 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/aa4f02194fe00a1a842949e448661e227f16fe9f
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Mon May 15 12:10:21 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/b4895dd40d1c34f345cd2807f7a81395ba27e8ee
Thu May 11 16:20:27 2017 +0200 ; move to jcommander ; https://github.com/lindenb/jvarkit/commit/15b6fabdbdd7ce0d1e20ca51e1c1a9db8574a59e
Sun May 7 13:21:47 2017 +0200 ; rm xml ; https://github.com/lindenb/jvarkit/commit/f37088a9651fa301c024ff5566534162bed8753d
Thu Apr 20 17:17:22 2017 +0200 ; continue transition jcommander ; https://github.com/lindenb/jvarkit/commit/fcf5def101925bea9ddd001d8260cf65aa52d6a0
Wed Feb 24 11:59:51 2016 +0100 ; jsonx2json, elixir registry ; https://github.com/lindenb/jvarkit/commit/fb6af381b43b9112360587dde45d0918c2b40665
Mon Nov 30 16:53:51 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/89f3cbe043ac8c52735feec5b45e43cf873b7179
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Fri Nov 28 12:44:44 2014 +0100 ; find all coverages ; https://github.com/lindenb/jvarkit/commit/a8c96e489787bf94d752e6bbd7c091175617459b
Thu Nov 27 13:11:06 2014 +0100 ; bam compare coverage ; https://github.com/lindenb/jvarkit/commit/0be60cca2b40fa2bb2713e759271573936911aba
Wed Nov 19 15:13:16 2014 +0100 ; fix qual in sam2tsv ; https://github.com/lindenb/jvarkit/commit/4be7535a661f42f43b44d9bd513694116cce189c
Wed Nov 19 12:40:59 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/5b49c9ecb5e06c8524830b79939ad0788558cf98
Wed Nov 19 11:48:09 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/6171c728bfe1fd2e2aca921424cede79dcc40b6f
Tue Nov 18 17:06:48 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/bf0bc5db943ad39514df4676074850e0cd9cc3ef
Fri May 23 15:32:54 2014 +0200 ; continue move to htsjdk ; https://github.com/lindenb/jvarkit/commit/b5a8a3bce5ecd952abffb7aae6223d1e03a9809e
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sat Apr 5 16:02:59 2014 +0200 ; sam2tsv with qual+flag ; https://github.com/lindenb/jvarkit/commit/df519652164b3a3a10b176de53d0f4186d689895
Wed Feb 12 18:02:27 2014 +0100 ; fastq grep added ; https://github.com/lindenb/jvarkit/commit/8d109ebd8d8fd928b58289f90a970d83e3ce474e
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Sat Dec 14 16:09:43 2013 +0100 ; updated sam2tsv ; https://github.com/lindenb/jvarkit/commit/c9e7c1d22439928326b381eabf087a8bd831ddc7
Fri Oct 25 17:42:45 2013 +0200 ; close Reference Fasta (picard.100) ; https://github.com/lindenb/jvarkit/commit/9c4a6831016175308ec9a80539e2093c32e78af9
Wed Jun 19 18:30:52 2013 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/032e974cd9a068db5c8aa74e2eba723033f073a8
Thu Jun 6 16:05:49 2013 +0200 ; compare 2 bams ; https://github.com/lindenb/jvarkit/commit/4fa8928e486e47e5f0c0bf94bf49859dabc2039c
Wed Jun 5 19:01:41 2013 +0200 ; aln format ; https://github.com/lindenb/jvarkit/commit/7c5a24f5ea82dbc9c66e8d538a00e008e5f6f97e
Tue Jun 4 15:20:17 2013 +0200 ; sam2tsv ; https://github.com/lindenb/jvarkit/commit/e81d4706dd51297677ddb64dcc69aaa681eab4af
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sam2tsv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Output

Columns are:

 *  read name
 *  read flags
 *  reference name
 *  read-pos
 *  read-base
 *  read-qual
 *  ref-pos
 *  ref-base
 *  cigar-op



### Example
 


```
$ java -jar dist/sam2tsv.jar -A  \
    -r samtools-0.1.18/examples/toy.fa 
      samtools-0.1.18/examples/toy.sam
r001	163	ref	0	T	.	7	T	M
r001	163	ref	1	T	.	8	T	M
r001	163	ref	2	A	.	9	A	M
r001	163	ref	3	G	.	10	G	M
r001	163	ref	4	A	.	11	A	M
r001	163	ref	5	T	.	12	T	M
r001	163	ref	6	A	.	13	A	M
r001	163	ref	7	A	.	14	A	M
r001	163	ref	8	A	.	.	.	I
r001	163	ref	9	G	.	.	.	I
r001	163	ref	10	A	.	.	.	I
r001	163	ref	11	G	.	.	.	I
r001	163	ref	12	G	.	15	G	M
r001	163	ref	13	A	.	16	A	M
r001	163	ref	14	T	.	17	T	M
r001	163	ref	15	A	.	18	A	M
r001	163	ref	.	.	.	19	G	D
r001	163	ref	16	C	.	20	C	M
r001	163	ref	17	T	.	21	T	M
r001	163	ref	18	G	.	22	G	M
:   ref        7 TTAGATAAAGAGGATA-CTG 22      
:                ||||||||    |||| |||
:  r001        1 TTAGATAA----GATAGCTG 19      
r002	0	ref	1	A	.	.	.	I
r002	0	ref	2	A	.	.	.	I
r002	0	ref	3	A	.	9	A	M
r002	0	ref	4	G	.	10	G	M
r002	0	ref	5	A	.	11	A	M
r002	0	ref	6	T	.	12	T	M
r002	0	ref	7	A	.	13	A	M
r002	0	ref	8	A	.	14	A	M
r002	0	ref	9	G	.	.	.	I
r002	0	ref	10	G	.	.	.	I
r002	0	ref	11	G	.	15	G	M
r002	0	ref	12	A	.	16	A	M
r002	0	ref	13	T	.	17	T	M  
(...)   

```


## Example 2

sam2tsv can read data from a linux pipe.

```
samtools view -h input.bam | java -jar dist/sam2tsv.jar
```




### History

 *  Moved to a standard argc/argv command line
 *  2014-04: added qual and samflag. Fixed a bug in soft-clip
 *  2014-11: manage hard+soft clip

### Citations


Sam2tsv was cited in : 

  * "Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements" . McCoy RC, Taylor RW, Blauwkamp TA, Kelley JL, Kertesz M, et al. (2014) Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements. PLoS ONE 9(9): e106689. doi: 10.1371/journal.pone.0106689  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0106689
  * "High-Throughput Identification of Genetic Variation Impact on pre-mRNA Splicing Efficiency". Scott I Adamson, Lijun Zhan, Brenton R Graveley. doi: [https://doi.org/10.1101/191122](https://doi.org/10.1101/191122).
  * "Linkage of A-to-I RNA editing in metazoans and the impact on genome evolution "  Molecular Biology and Evolution, msx274, https://doi.org/10.1093/molbev/msx274



