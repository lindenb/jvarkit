# Biostar160470

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Getting untranslated nucleotide sequences on tblastn standalone 


## Usage

```
Usage: biostar160470 [options] Files
  Options:
    -p, --bindir
      Blast binaries path
    -d, --db
      Blast db name
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

 * blasn
 * blast
 * translation
 * protein



## See also in Biostars

 * [https://www.biostars.org/p/160470](https://www.biostars.org/p/160470)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar160470
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar160470.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar160470.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar160470** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

Makefile:

```make
bin.dir=/commun/data/packages/ncbi/ncbi-blast-2.2.28+/bin

all: blastdb.nin
	cat roxan.fa | ${bin.dir}/tblastn -db blastdb -outfmt 5 | java -jar biostar160470.jar -p ${bin.dir} -d blastdb| xmllint --format - 

blastdb.nin: mysequences.fa
	${bin.dir}/makeblastdb -dbtype nucl -in $< -out blastdb
```

ouput:
```xml
              <Hsp>
                <Hsp_num>5</Hsp_num>
                <Hsp_bit-score>31.9574</Hsp_bit-score>
                <Hsp_score>71</Hsp_score>
                <Hsp_evalue>0.000226217</Hsp_evalue>
                <Hsp_query-from>520</Hsp_query-from>
                <Hsp_query-to>575</Hsp_query-to>
                <Hsp_hit-from>1711</Hsp_hit-from>
                <Hsp_hit-to>1860</Hsp_hit-to>
                <Hsp_query-frame>0</Hsp_query-frame>
                <Hsp_hit-frame>1</Hsp_hit-frame>
                <Hsp_identity>16</Hsp_identity>
                <Hsp_positive>27</Hsp_positive>
                <Hsp_gaps>6</Hsp_gaps>
                <Hsp_align-len>56</Hsp_align-len>
                <Hsp_qseq>MGEFRLCDRLQKGKACPDGDKCRCAHGQEELNEWLDRREVLKQKLAKARKDMLLCP</Hsp_qseq>
                <Hsp_hseq>VGSYYLCKDMINKQDCKYGDNCTFAYHQEEIDVWTEERK------GTLNRDLLFDP</Hsp_hseq>
                <Hsp_midline>+G + LC  +   + C  GD C  A+ QEE++ W + R+          +D+L  P</Hsp_midline>
                <Hsp_hit-DNA>GTGGGCTCCTACTACCTGTGCAAAGACATGATTAACAAGCAGGACTGTAAGTACGGGGATAACTGCACCTTCGCCTACCATCAGGAGGAGATCGACGTGTGGACCGAGGAGCGGAAG------------------CTGCTCTTCGACCCG</Hsp_hit-DNA>
              </Hsp>
              <Hsp>
                <Hsp_num>6</Hsp_num>
                <Hsp_bit-score>22.3274</Hsp_bit-score>
                <Hsp_score>46</Hsp_score>
                <Hsp_evalue>0.215374</Hsp_evalue>
                <Hsp_query-from>22</Hsp_query-from>
                <Hsp_query-to>62</Hsp_query-to>
                <Hsp_hit-from>3316</Hsp_hit-from>
                <Hsp_hit-to>3435</Hsp_hit-to>
                <Hsp_query-frame>0</Hsp_query-frame>
                <Hsp_hit-frame>1</Hsp_hit-frame>
                <Hsp_identity>15</Hsp_identity>
                <Hsp_positive>19</Hsp_positive>
                <Hsp_gaps>1</Hsp_gaps>
                <Hsp_align-len>41</Hsp_align-len>
                <Hsp_qseq>HEAPWTNLTPSWRRPTHRTTVPLAVLRNQPPRQSPACPTLP</Hsp_qseq>
                <Hsp_hseq>HQAAPSPLRPCPSSPHHRPGVRTQAHVLQPP-EAPLKPGLP</Hsp_hseq>
                <Hsp_midline>H+A  + L P    P HR  V       QPP ++P  P LP</Hsp_midline>
                <Hsp_hit-DNA>CATCAGGCAGCCCCCAGCCCCCTGAGGCCCTGTCCATCTTCTCCCCACCACCGCCCCGGTGTGCGTACCCAGGCGCACGTGCTGCAGCCCCCG---GCCCCGCTGAAACCTGGGCTGCCC</Hsp_hit-DNA>
              </Hsp>
```
