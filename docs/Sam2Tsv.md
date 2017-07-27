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

 *  Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements* . McCoy RC, Taylor RW, Blauwkamp TA, Kelley JL, Kertesz M, et al. (2014) Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements. PLoS ONE 9(9): e106689. doi: 10.1371/journal.pone.0106689  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0106689



