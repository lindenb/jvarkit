# SamJavascript

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filters a BAM using a javascript expression ( java nashorn engine  ).


## Usage

```
Usage: samjs [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -e, --expression
      javascript expression
    -X, --fail
      Save dicarded reads in that file
    -f, --file
      javascript file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -N, --limit
      limit to 'N' records (for debugging).
      Default: -1
    -o, --output
      Output file. Optional . Default: stdout
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
 * nashorn
 * javascript
 * filter



## See also in Biostars

 * [https://www.biostars.org/p/75168](https://www.biostars.org/p/75168)
 * [https://www.biostars.org/p/81750](https://www.biostars.org/p/81750)
 * [https://www.biostars.org/p/75354](https://www.biostars.org/p/75354)
 * [https://www.biostars.org/p/77802](https://www.biostars.org/p/77802)
 * [https://www.biostars.org/p/103052](https://www.biostars.org/p/103052)
 * [https://www.biostars.org/p/106900](https://www.biostars.org/p/106900)
 * [https://www.biostars.org/p/150530](https://www.biostars.org/p/150530)
 * [https://www.biostars.org/p/253774](https://www.biostars.org/p/253774)
 * [https://www.biostars.org/p/256615](https://www.biostars.org/p/256615)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samjs
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJavascript.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJavascript.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samjs** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

 * "bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734).


## Motivation

Filters a BAM using javascript( java nashorn engine).

For each read the script injects in the context the following values:


* **'record'** a SamRecord  [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html)
* **'header'** a SAMFileHeader  [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html)


the script should return a boolean : true accept the read, false: discard the read.

## Example

### Example 1


get a SAM where the  read OR the mate is unmapped

```bash
java -jar dist/samjs.jar  \
	-e "record.readUnmappedFlag || record.mateUnmappedFlag;" \
	ex1.bam

@HD	VN:1.4	SO:unsorted
@SQ	SN:seq1	LN:1575
@SQ	SN:seq2	LN:1584
B7_591:4:96:693:509	73	seq1	1	99	36M	*	0	0	CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG	<<<<<<<<<<<<<<<;<<<<<<<<<5<<<<<;:<;7	H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:73
EAS54_65:7:152:368:113	73	seq1	3	99	35M	*	0	0	CTAGTGGCTCATTGTAAATGTGTGGTTTAACTCGT	<<<<<<<<<<0<<<<655<<7<<<:9<<3/:<6):H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:66
EAS51_64:8:5:734:57	137	seq1	5	99	35M	*	0	0	AGTGGCTCATTGTAAATGTGTGGTTTAACTCGTCC	<<<<<<<<<<<7;71<<;<;;<7;<<3;);3*8/5H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:66
B7_591:1:289:587:906	137	seq1	6	63	36M	*	0	0	GTGGCTCATTGTAATTTTTTGTTTTAACTCTTCTCT	(-&----,----)-)-),'--)---',+-,),''*,	H0:i:0	H1:i:0	MF:i:130	NM:i:5	UQ:i:38	Aq:i:63
EAS56_59:8:38:671:758	137	seq1	9	99	35M	*	0	0	GCTCATTGTAAATGTGTGGTTTAACTCGTCCATGG	<<<<<<<<<<<<<<<;<;7<<<<<<<<7<<;:<5%H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:72
EAS56_61:6:18:467:281	73	seq1	13	99	35M	*	0	0	ATTGTAAATGTGTGGTTTAACTCGTCCCTGGCCCA	<<<<<<<<;<<<8<<<<<;8:;6/686&;(16666H0:i:0	H1:i:1	MF:i:18	NM:i:1	UQ:i:5	Aq:i:39
EAS114_28:5:296:340:699	137	seq1	13	99	36M	*	0	0	ATTGTAAATGTGTGGTTTAACTCGTCCATGGCCCAG	<<<<<;<<<;<;<<<<<<<<<<<8<8<3<8;<;<0;	H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:73
B7_597:6:194:894:408	73	seq1	15	99	35M	*	0	0	TGTAAATGTGTGGTTTAACTCGTCCATTGCCCAGC	<<<<<<<<<7<<;<<<<;<<<7;;<<<*,;;572<H0:i:0	H1:i:1	MF:i:18	NM:i:1	UQ:i:9	Aq:i:43
EAS188_4:8:12:628:973	89	seq1	18	75	35M	*	0	0	AAATGTGTGGTTTAACTCGTCCATGGCCCAGCATT	==;=:;:;;:====;=;===:=======;==;===H0:i:1	H1:i:0	MF:i:64	NM:i:0	UQ:i:0	Aq:i:0
(...)
```

### Example 2

remove reads with indels:

```
java -jar dist/samjs.jar -e 'function accept(r) { if(r.getReadUnmappedFlag()) return false; var cigar=r.getCigar();if(cigar==null) return false; for(var i=0;i< cigar.numCigarElements();++i) {if(cigar.getCigarElement(i).getOperator().isIndelOrSkippedRegion()) return false; } return true;} accept(record);' input.bam
```


