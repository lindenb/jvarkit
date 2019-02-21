# PcrClipReads

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Soft clip bam files based on PCR target regions


## Usage

```
Usage: pcrclipreads [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
  * -B, --bed
      Bed file containing non-overlapping PCR fragments
    -flag, --flag
      Only run on reads having sam flag (flag). -1 = all reads. (as 
      https://github.com/lindenb/jvarkit/issues/43) 
      Default: -1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -largest, --largest
      see if a read overlaps two bed intervals use the bed region sharing the 
      longest sequence with a read. see 
      https://github.com/lindenb/jvarkit/issues/44 
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -pr, --programId
      add a program group PG to the clipped SAM records
      Default: false
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * pcr
 * bed



## See also in Biostars

 * [https://www.biostars.org/p/147136](https://www.biostars.org/p/147136)
 * [https://www.biostars.org/p/178308](https://www.biostars.org/p/178308)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew pcrclipreads
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pcr/PcrClipReads.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pcr/PcrClipReads.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pcr/PcrClipReadsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/pcr/PcrClipReadsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pcrclipreads** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




 Soft clip BAM files based on PCR target regions https://www.biostars.org/p/147136/


 *  mapping quality is set to zero if a read on mapped strand - overlap the 5' side of the PCR fragment
 *  mapping quality is set to zero if a read on mapped strand + overlap the 3' side of the PCR fragment
 *  mapping quality is set to zero if no PCR fragment is found


after processing the BAM file should be sorted, processed with samtools calmd and picard fixmate


### Example


```
echo  "seq2\t1100\t1200" > test.bed
java -jar dist/pcrclipreads.jar -B test.bed  samtools-0.1.19/examples/ex1.bam  |\
	samtools  view -q 1 -F 4 -Sbu  -  |\
	samtools  sort -o clipped.bam -  && samtools index clipped.bam

samtools tview -p seq2:1100  clipped.bam  samtools-0.1.19/examples/ex1.fa

```


### output


![img](http://i.imgur.com/bjDEnMW.jpg)



```
    1091      1101      1111      1121      1131      1141      1151      1161      1171      1181      1191
AAACAAAGGAGGTCATCATACAATGATAAAAAGATCAATTCAGCAAGAAGATATAACCATCCTACTAAATACATATGCACCTAACACAAGACTACCCAGATTCATAAAACAAATNNNNN
              ...................................                               ..................................
              ,,,                                                               ..................................
              ,,,,,                                                              .................................
              ,,,,,,,,,,,                                                        .............................N...
              ,,,,,,,,                                                             ...............................
              ,,g,,,,,,,,,,,,,                                                        ............................
              ,,,,,,,,,,,,,,,,,,,,                                                    ............................
              ,,,,,,,,,,,,,,,,,,,                                                       ..........................
              ,,,,,,,,,,,,,,,,,,,,,,                                                    ..........................
              ,,,,,,,,,,,,,,,,,,,,,,,                                                       ......................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,                                                        ..................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,                                                        ..................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                       .................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                       ................
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                        ...............
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                         ............
              ,,,,,,,,,,,,a,,,,,,,,,,,,,,,,,,,                                                             .......
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                            ......
              ,,a,,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                              ....
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                             ....
              ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                                                                .
                                                                                                                 .

```

## Cited in

 * BAMClipper: removing primers from alignments to minimize false-negative mutations in amplicon next-generation sequencing [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5431517/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5431517/)


## History

 * 20170630 : rewritten after [https://github.com/lindenb/jvarkit/issues/81](https://github.com/lindenb/jvarkit/issues/81)



