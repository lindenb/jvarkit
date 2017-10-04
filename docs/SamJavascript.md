# SamJavascript

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
      What kind of help
      Possible Values: [usage, markdown, xml]
    -N, --limit
      limit to 'N' records (for debugging).
      Default: -1
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: TypeImpl{name='SAM', fileExtension='sam', indexExtension='null'}
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

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make samjs
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJavascript.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamJavascript.java)


<details>
<summary>Git History</summary>

```
Mon Aug 7 15:05:18 2017 +0200 ; samjdk ; https://github.com/lindenb/jvarkit/commit/93cb0448be4d6deb253b21620d1da63ad2be9475
Wed Jun 7 16:15:41 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/8ed068b24f8b7a949e10999c5b219401827c03b0
Sat May 20 16:46:01 2017 +0200 ; plot case ctrl ; https://github.com/lindenb/jvarkit/commit/f03530037d4b72cb515cfdf36cfd604c3654a760
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Wed Apr 26 17:26:23 2017 +0200 ; cont jcommander ; https://github.com/lindenb/jvarkit/commit/ab6c7b760cd5376e08da24426cede7f84a6b3ae2
Fri Apr 14 17:09:04 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/881f52d5d3775325240114702b7f07148b626f4c
Fri Nov 27 15:22:25 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/04a83d5b9f0e69fd2f7087e519b0de3e2b4f9863
Tue Jul 22 18:27:05 2014 +0200 ; htsjdk version ; https://github.com/lindenb/jvarkit/commit/3780ec67df1786dc87b1d5a06c35c1c3d473446c
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Tue Feb 11 15:03:41 2014 +0100 ; fixed bug in samjs ; https://github.com/lindenb/jvarkit/commit/3c67063c3c7091c41dae684c324d99376d226dd2
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Fri Dec 20 18:56:17 2013 +0100 ; samjs force unmatched ; https://github.com/lindenb/jvarkit/commit/007454d554ce6ca16592fd19e00a17dd1d95504e
Fri Dec 20 19:16:58 2013 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/e285c4864a128ad94e0d4cd025905328e59037c6
Fri Dec 20 16:43:59 2013 +0100 ; samjs and option -X' ; https://github.com/lindenb/jvarkit/commit/742a0716f089297af9020891cd12ed7ee284e249
Tue Dec 10 13:54:40 2013 +0100 ; vcfshuffle ; https://github.com/lindenb/jvarkit/commit/813ec64f4e5105e0cbdca7b7fabea70924381896
Fri Jun 28 12:36:50 2013 +0200 ; Merge branch 'master' of https://github.com/lindenb/jvarkit ; https://github.com/lindenb/jvarkit/commit/a0fc76ad8ea4841be2094c053bb516ac077a8e01
Fri Jun 28 12:36:47 2013 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/17278939a48c5827acb24286f9241e83cf881947
Thu Jun 27 15:53:34 2013 +0200 ; input for deseq ; https://github.com/lindenb/jvarkit/commit/37b47d801ecdd2cc03f50c1b8aab15fd432500c0
Tue Jun 25 09:30:11 2013 +0200 ; move samjs to picard ; https://github.com/lindenb/jvarkit/commit/7a3e3b0bc7297e40a1511fb54fbfbb6896bc10ac
Mon May 6 18:56:46 2013 +0200 ; moving to git ; https://github.com/lindenb/jvarkit/commit/55158d13f0950f16c4a3cc3edb92a87905346ee1
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samjs** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


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



