# MergeSplittedBlast

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

merge blast Hits from splitted BLAST database


## Usage

```
Usage: mergesplittedblast [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --version
      print version and exit

```


## Keywords

 * blast



## See also in Biostars

 * [https://www.biostars.org/p/90186](https://www.biostars.org/p/90186)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew mergesplittedblast
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/MergeSplittedBlast.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/MergeSplittedBlast.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **mergesplittedblast** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 

## Example

### Creating the blast database
the long sequence must be splitted using a sliding+overlapping window.
The fasta  header must be formatted as:
```
(chromName):(1-based-start-position)-(1-based-end-position):(chromLength)
```
### Makefile for test:
```make
BLASTBIN=/commun/data/packages/ncbi/ncbi-blast-2.2.28+/bin

all:blast.xml

blast.xml:query.fa database.fa
        ${BLASTBIN}/blastn -db database.fa -query query.fa -outfmt 5 -out $@ -dust no

query.fa:
        echo -e ">q1\nCCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACACAGACATCATAACAAA" > $@

database.fa:
        curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrM.fa.gz" |\
        gunzip -c | grep -v '^>' | tr -d "\n" |\
        awk '{L=length($$0);for(i=1;i+77<=L;i+=65) printf(">chrM:%d-%d:%d\n%s\n",i,i+77,L,substr($$0,i,77));}' |\
        fold -w 50 > $@ && \
        ${BLASTBIN}/makeblastdb -in $@ -dbtype nucl
```
The blast.xml produced is:
```xml
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.2.28+</BlastOutput_version>
  <BlastOutput_reference>Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), "A greedy algorithm for aligning DNA sequences", J Comput Biol 2000; 7(1-2):203-14.</BlastOutput_reference>
  <BlastOutput_db>database.fa</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>q1</BlastOutput_query-def>
  <BlastOutput_query-len>123</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>1</Parameters_sc-match>
      <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
      <Parameters_gap-open>0</Parameters_gap-open>
      <Parameters_gap-extend>0</Parameters_gap-extend>
      <Parameters_filter>m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>q1</Iteration_query-def>
      <Iteration_query-len>123</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|3</Hit_id>
          <Hit_def>chrM:196-273:16571</Hit_def>
          <Hit_accession>3</Hit_accession>
          <Hit_len>77</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>143.312</Hsp_bit-score>
              <Hsp_score>77</Hsp_score>
              <Hsp_evalue>1.29152e-37</Hsp_evalue>
              <Hsp_query-from>31</Hsp_query-from>
              <Hsp_query-to>107</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>77</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>77</Hsp_identity>
              <Hsp_positive>77</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>77</Hsp_align-len>
              <Hsp_qseq>TACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACA</Hsp_qseq>
              <Hsp_hseq>TACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACA</Hsp_hseq>
              <Hsp_midline>|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|2</Hit_id>
          <Hit_def>chrM:131-208:16571</Hit_def>
          <Hit_accession>2</Hit_accession>
          <Hit_len>77</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>78.6796</Hsp_bit-score>
              <Hsp_score>42</Hsp_score>
              <Hsp_evalue>3.69397e-18</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>42</Hsp_query-to>
              <Hsp_hit-from>36</Hsp_hit-from>
              <Hsp_hit-to>77</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>42</Hsp_identity>
              <Hsp_positive>42</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>42</Hsp_align-len>
              <Hsp_qseq>CCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTG</Hsp_qseq>
              <Hsp_hseq>CCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTG</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>3</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|4</Hit_id>
          <Hit_def>chrM:261-338:16571</Hit_def>
          <Hit_accession>4</Hit_accession>
          <Hit_len>77</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>52.8265</Hsp_bit-score>
              <Hsp_score>28</Hsp_score>
              <Hsp_evalue>2.23898e-10</Hsp_evalue>
              <Hsp_query-from>96</Hsp_query-from>
              <Hsp_query-to>123</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>28</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>28</Hsp_identity>
              <Hsp_positive>28</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>28</Hsp_align-len>
              <Hsp_qseq>CCGCTTTCCACACAGACATCATAACAAA</Hsp_qseq>
              <Hsp_hseq>CCGCTTTCCACACAGACATCATAACAAA</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>254</Statistics_db-num>
          <Statistics_db-len>19558</Statistics_db-len>
          <Statistics_hsp-len>13</Statistics_hsp-len>
          <Statistics_eff-space>1788160</Statistics_eff-space>
          <Statistics_kappa>0.46</Statistics_kappa>
          <Statistics_lambda>1.28</Statistics_lambda>
          <Statistics_entropy>0.85</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
```
### Process the XML with MergeSplittedBlast:
```bash
$ java -jar dist/mergesplittedblast.jar blast.xml |\
  xmllint --format - 
```
### Output:
```xml
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.2.28+</BlastOutput_version>
  <BlastOutput_reference>Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), "A greedy algorithm for aligning DNA sequences", J Comput Biol 2000; 7(1-2):203-14.</BlastOutput_reference>
  <BlastOutput_db>database.fa</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>q1</BlastOutput_query-def>
  <BlastOutput_query-len>123</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>1</Parameters_sc-match>
      <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
      <Parameters_gap-open>0</Parameters_gap-open>
      <Parameters_gap-extend>0</Parameters_gap-extend>
      <Parameters_filter>m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>q1</Iteration_query-def>
      <Iteration_query-len>123</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|3</Hit_id>
          <Hit_def>chrM</Hit_def>
          <Hit_accession>3</Hit_accession>
          <Hit_len>16571</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>78.6796</Hsp_bit-score>
              <Hsp_score>123</Hsp_score>
              <Hsp_evalue>3.69397e-18</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>123</Hsp_query-to>
              <Hsp_hit-from>166</Hsp_hit-from>
              <Hsp_hit-to>288</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>42</Hsp_identity>
              <Hsp_positive>42</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>123</Hsp_align-len>
              <Hsp_qseq>CCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACACAGACATCATAACAAA</Hsp_qseq>
              <Hsp_hseq>CCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACACAGACATCATAACAAA</Hsp_hseq>
              <Hsp_midline>|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
```
