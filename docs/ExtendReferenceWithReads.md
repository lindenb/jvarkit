# ExtendReferenceWithReads

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extending ends of sequences with the help of reads


## Usage

```
Usage: extendrefwithreads [options] Files
  Options:
    -f, --callingfraction
      (0.0<float<=1.0) new base must have fraction greater than this number
      Default: 0.8
    -filter, --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -N, --mincontig
      onsider only gaps in reference with size&gt;=N
      Default: 100
    -d, --mindepth
      min depth
      Default: 1
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * read
 * fastq
 * reference
 * sam
 * bam



## See also in Biostars

 * [https://www.biostars.org/p/148089](https://www.biostars.org/p/148089)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew extendrefwithreads
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/extendref/ExtendReferenceWithReads.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/extendref/ExtendReferenceWithReads.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **extendrefwithreads** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
 ## Example

```
$  java   -jar dist/extendrefwithreads.jar \
     -R human_g1k_v37.fasta -f 0.3 \
     f1.bam f2.bam f3.bam 2> /dev/null |\
  cat -n | grep -E '(>|[atgc])' 

     1	>1
   167	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNncgattaccctaacgctcac
   168	cctaaccctcnccctntnccnncnncccnncttcttccgaTAACCCTAACCCTAACCCTA
  3791	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNatt
  3792	tatgcNctttntgctgtGATTCATGGCTGAAATCGTGTTTGACCAGCTATGTGTGTCTCT
  8691	NNNNNNNNNNNNNNNNNNNNNNNNctagGATCCTTGAAGCGCCCCCAAGGGCATCTTCTC
 64089	TGGTGAGGGAAATTAGAACCACGACAATTTGGGAACTTAGCTTCTGCCctgctccNNNNN
 66589	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNgagtAGCTGAGACTAC
 
 ```

