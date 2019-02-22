# Sam2Tsv

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Prints the SAM alignments as a TAB delimited file.


## Usage

```
Usage: sam2tsv [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
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

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew sam2tsv
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2Tsv.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2Tsv.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2TsvTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2TsvTest.java)


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
 *  2019-02 : manage reads without qualities, contig name converter

### Citations


Sam2tsv was cited in : 

  * "Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements" . McCoy RC, Taylor RW, Blauwkamp TA, Kelley JL, Kertesz M, et al. (2014) Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements. PLoS ONE 9(9): e106689. doi: 10.1371/journal.pone.0106689  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0106689
  * "High-Throughput Identification of Genetic Variation Impact on pre-mRNA Splicing Efficiency". Scott I Adamson, Lijun Zhan, Brenton R Graveley. doi: [https://doi.org/10.1101/191122](https://doi.org/10.1101/191122).
  * "Linkage of A-to-I RNA editing in metazoans and the impact on genome evolution "  Molecular Biology and Evolution, msx274, https://doi.org/10.1093/molbev/msx274
  * "Vex-seq: high-throughput identification of the impact of genetic variation on pre-mRNA splicing efficiency" Genome Biology201819:71 https://doi.org/10.1186/s13059-018-1437-x
  * "Accurate detection of m6A RNA modifications in native RNA sequences" Huanle Liu, Oguzhan Begik, Morghan C Lucas, Christopher E Mason, Schraga Schwartz, John S Mattick, Martin A Smith, Eva Maria Novoa bioRxiv 525741; doi: https://doi.org/10.1101/525741 

