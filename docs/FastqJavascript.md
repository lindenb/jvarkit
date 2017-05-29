# FastqJavascript

Filters a FASTQ file using javascript( java nashorn engine). 


## Usage

```
Usage: fastqjs [options] Files
  Options:
    -R1, --R1
      for paired/interleaved input, save R1 reads in this file
    -R2, --R2
      for paired/interleaved input, save R2 reads in this file
    -X, --failing
      Save dicarded reads in that file. Optional. Default: no file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -i, --interleaved
      interleaved input
      Default: false
    -N, --limit
      limit to 'N' records -1:all
      Default: -1
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -e
      javascript expression
    -f
      javascript file

```


## Keywords

 * fastq


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
$ make fastqjs
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqJavascript.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqJavascript.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastqjs** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Motivation


The script puts 'rec' a FastqRecord, or 'pair' for an interleaved input, into the script context 



## Example

Find pairs of fastq where both reads contains a **PmeI** restriction site ( `GTTT/AAAC` )

```bash
$ paste <(gunzip -c F.fastq.gz | paste - - - -)  <(gunzip -c R.fastq.gz | paste - - - -) |\
  tr "\t" "\n" |\
  java -jar dist/fastqjs.jar -i -e 'pair.get(0).getReadString().contains("GTTTAAAC") && pair.get(1).getReadString().contains("GTTTAAAC") '

@HWI-1KL149:13:C0RNFACXX:8:1309:5373:60519 1:N:0:CTTGTA
TTCCAAAAATGTTTAAACTTTACAAATTTTCTTTCTGCAAAGGATATTTAAAACTTTGTCAAGACAAATATAAAAGTCTGTTCTTTTCATTAGTCTCTATA
+
CCCFFFFFHHHHHJJJJJJJJJBHIIJJJJJJJJJJIJIJJJJJJJJJJJJJJJJJJJJJJJIJJJJJJJJIJJIFEHIJHIJHHHHHHHFFFFFFFEEDE
@HWI-1KL149:13:C0RNFACXX:8:1309:5373:60519 2:N:0:CTTGTA
ACGCTTGATATTTGGTTTAAACATTTCTTGATTCAGAGAAGGTAGATGGTTATAGAGACTAATGAAAAGAACAGACTTTTATATTTGTCTTGACAAAGTTT
+
CCCFFFFFHHHHHJIIIJJIJJJJJIJJJJEHIJJJGIIIJJ?FHGIGIDDFHJJJJJJJIIJJIJJJJJJJJJIJHHHHHHHFFFFFFFEEEEEEDDDDC
@HWI-1KL149:13:C0RNFACXX:8:1309:6861:76085 1:N:0:CTTGTA
ACAGTATATCTATGTGAAAGTTAAAAAGAAATCGCTGTTTAGATGGAAGATGAGACCAGGTTATCATAGTTTTAGAAGAGGAGTTTAAACTTCATGCAGTG
+
CCCFDFFFHHHHHJIJJJJJIIJJJJJJJJJJJJJJJIJJJJIJJJJJJJJIIJJJJJJJFHIJJJJJJHIIJIJJHHHHHFFFFEEEEEEEDDEDDDDDD
@HWI-1KL149:13:C0RNFACXX:8:1309:6861:76085 2:N:0:CTTGTA
GGGTAGTCCACAAACAATGTGTTCATGTTGTCTCCCTCTTACTCACAAGCCTTCTGTAGCTCCCAACATTCACTGCATGAAGTTTAAACTCCTCTTCTAAA
+
@CCDDDFFHHHHHJJIJJJIJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJGIIJIJIIJJJJJJJJJIJJJIJJHHHHHHFFFFFFCEEEEEEDDDDDDED

```



