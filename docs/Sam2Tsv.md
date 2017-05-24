# Sam2Tsv

Prints the SAM alignments as a TAB delimited file.


## Usage

```
Usage: sam2tsv [options] Files
  Options:
    -h, --help
      print help and exit
    -o, --output
      Output file. Optional . Default: stdout
    -A, --printAlignments
      Print Alignments
      Default: false
  * -r, -R, --reference
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2Tsv.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2Tsv.java
)
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
r002	0	ref	14	A	.	18	A	M
r002	0	ref	15	A	.	.	.	I
r002	0	ref	16	A	.	.	.	I
:   ref        8 AAAGATAAGGGATAAA 18      
:                  ||||||  ||||  
:  r002        1 --AGATAA--GATA-- 17      
r003	0	ref	0	A	.	9	A	M
r003	0	ref	1	G	.	10	G	M
r003	0	ref	2	C	.	11	A	M
r003	0	ref	3	T	.	12	T	M
r003	0	ref	4	A	.	13	A	M
r003	0	ref	5	A	.	14	A	M
:   ref        4 AGCTAA 14      
:                || |||
:  r003        1 AGATAA 6       
r004	0	ref	0	A	.	16	A	M
r004	0	ref	1	T	.	17	T	M
r004	0	ref	2	A	.	18	A	M
r004	0	ref	3	G	.	19	G	M
r004	0	ref	4	C	.	20	C	M
r004	0	ref	5	T	.	21	T	M
r004	0	ref	.	.	.	22	G	N
r004	0	ref	.	.	.	23	T	N
r004	0	ref	.	.	.	24	G	N
r004	0	ref	.	.	.	25	C	N
r004	0	ref	.	.	.	26	T	N
r004	0	ref	.	.	.	27	A	N
r004	0	ref	.	.	.	28	G	N
r004	0	ref	.	.	.	29	T	N
r004	0	ref	.	.	.	30	A	N
r004	0	ref	.	.	.	31	G	N
r004	0	ref	.	.	.	32	G	N
r004	0	ref	.	.	.	33	C	N
r004	0	ref	.	.	.	34	A	N
r004	0	ref	.	.	.	35	G	N
r004	0	ref	6	C	.	.	.	I
r004	0	ref	7	T	.	36	T	M
r004	0	ref	8	C	.	37	C	M
r004	0	ref	9	A	.	38	A	M
r004	0	ref	10	G	.	39	G	M
r004	0	ref	11	C	.	40	C	M
:   ref       16 ATAGCT--------------CTCAGC 40      
:                ||||||               |||||
:  r004        1 ATAGCTGTGCTAGTAGGCAG-TCAGC 12      
r003	16	ref	0	T	.	29	T	M
r003	16	ref	1	A	.	30	A	M
r003	16	ref	2	G	.	31	G	M
r003	16	ref	3	G	.	32	G	M
r003	16	ref	4	C	.	33	C	M
:   ref       23 TAGGC 33      
:                |||||
:  r003        1 TAGGC 5       
r001	83	ref	0	C	.	37	C	M
r001	83	ref	1	A	.	38	A	M
r001	83	ref	2	G	.	39	G	M
r001	83	ref	3	C	.	40	C	M
r001	83	ref	4	G	.	41	G	M
r001	83	ref	5	C	.	42	C	M
r001	83	ref	6	C	.	43	C	M
r001	83	ref	7	A	.	44	A	M
r001	83	ref	8	T	.	45	T	M
:   ref       37 CAGCGCCAT 45      
:                |||||||||
:  r001        1 CAGCGCCAT 9       
x1	0	ref2	0	A	30	1	a	M
x1	0	ref2	1	G	30	2	g	M
x1	0	ref2	2	G	30	3	g	M
x1	0	ref2	3	T	30	4	t	M
x1	0	ref2	4	T	30	5	t	M
x1	0	ref2	5	T	30	6	t	M
x1	0	ref2	6	T	30	7	t	M
x1	0	ref2	7	A	30	8	a	M
x1	0	ref2	8	T	30	9	t	M
x1	0	ref2	9	A	30	10	a	M
x1	0	ref2	10	A	30	11	a	M
x1	0	ref2	11	A	30	12	a	M
x1	0	ref2	12	A	30	13	a	M
x1	0	ref2	13	C	30	14	c	M
x1	0	ref2	14	A	30	15	a	M
x1	0	ref2	15	A	30	16	a	M
x1	0	ref2	16	A	30	17	t	M
x1	0	ref2	17	T	30	18	t	M
x1	0	ref2	18	A	30	19	a	M
x1	0	ref2	19	A	30	20	a	M
:  ref2        1 AGGTTTTATAAAACAAATAA 20      
:                |||||||||||||||| |||
:    x1        1 aggttttataaaacaattaa 20      
x2	0	ref2	0	G	30	2	g	M
x2	0	ref2	1	G	30	3	g	M
x2	0	ref2	2	T	30	4	t	M
x2	0	ref2	3	T	30	5	t	M
x2	0	ref2	4	T	30	6	t	M
x2	0	ref2	5	T	30	7	t	M
x2	0	ref2	6	A	30	8	a	M
x2	0	ref2	7	T	30	9	t	M
x2	0	ref2	8	A	30	10	a	M
x2	0	ref2	9	A	30	11	a	M
x2	0	ref2	10	A	30	12	a	M
x2	0	ref2	11	A	30	13	a	M
x2	0	ref2	12	C	30	14	c	M
x2	0	ref2	13	A	30	15	a	M
x2	0	ref2	14	A	30	16	a	M
x2	0	ref2	15	A	30	17	t	M
x2	0	ref2	16	T	30	18	t	M
x2	0	ref2	17	A	30	19	a	M
x2	0	ref2	18	A	30	20	a	M
x2	0	ref2	19	T	30	21	g	M
x2	0	ref2	20	T	30	22	t	M
:  ref2        2 GGTTTTATAAAACAAATAATT 22      
:                ||||||||||||||| ||| |
:    x2        1 ggttttataaaacaattaagt 21      
(...)   

```






### History

 *  Moved to a standard argc/argv command line
 *  2014-04: added qual and samflag. Fixed a bug in soft-clip
 *  2014-11: manage hard+soft clip

### Citations


Sam2tsv was cited in : 

 *  Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements* . McCoy RC, Taylor RW, Blauwkamp TA, Kelley JL, Kertesz M, et al. (2014) Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements. PLoS ONE 9(9): e106689. doi: 10.1371/journal.pone.0106689  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0106689



