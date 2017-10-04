# SamGrep

grep read-names in a bam file


## Usage

```
Usage: samgrep [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -V, --invert
      invert
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -f, --readfile
      file containing a list of read names
    -R, --readname
      add the read name
      Default: []
    --samoutputformat
      Sam output format.
      Default: TypeImpl{name='SAM', fileExtension='sam', indexExtension='null'}
    -n, --stopafter
      when found, remove the read from the list of names when found more that 
      'n' time (increase speed)
      Default: -1
    -x, --tee
      if output fileame specified, continue to output original input to 
      stdout. 
      Default: false
    --version
      print version and exit

```


## Keywords

 * sam
 * bam


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
$ make samgrep
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samgrep/SamGrep.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samgrep/SamGrep.java)


<details>
<summary>Git History</summary>

```
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Fri May 12 19:41:30 2017 +0200 ; fix make, empty doc ; https://github.com/lindenb/jvarkit/commit/52fcf6d46a779fd7153ebc032fae643d2e266e7e
Thu May 11 16:20:27 2017 +0200 ; move to jcommander ; https://github.com/lindenb/jvarkit/commit/15b6fabdbdd7ce0d1e20ca51e1c1a9db8574a59e
Sun May 7 13:21:47 2017 +0200 ; rm xml ; https://github.com/lindenb/jvarkit/commit/f37088a9651fa301c024ff5566534162bed8753d
Wed Apr 26 17:26:23 2017 +0200 ; cont jcommander ; https://github.com/lindenb/jvarkit/commit/ab6c7b760cd5376e08da24426cede7f84a6b3ae2
Fri Mar 25 17:18:27 2016 +0100 ; sammask ; https://github.com/lindenb/jvarkit/commit/b9c834afec6c7c9904baecd2fb2b61e57261da0f
Mon Jun 8 17:24:41 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/f9a941d604f378ff40a32666c8381cb2450c7cfa
Fri May 23 15:32:54 2014 +0200 ; continue move to htsjdk ; https://github.com/lindenb/jvarkit/commit/b5a8a3bce5ecd952abffb7aae6223d1e03a9809e
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Wed Feb 12 18:02:27 2014 +0100 ; fastq grep added ; https://github.com/lindenb/jvarkit/commit/8d109ebd8d8fd928b58289f90a970d83e3ce474e
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Sun Jan 5 16:10:56 2014 +0100 ; vcf set dict ; https://github.com/lindenb/jvarkit/commit/f023bc9b0685266627a260c67813e7b76d42bef1
Sat Jan 4 16:27:34 2014 +0100 ; vcf merge, peek vcf, samgrep cmd line, sw factory ; https://github.com/lindenb/jvarkit/commit/bd9a33f08c6e32efe9b54e1909527dc57dc55060
Wed Oct 30 13:27:25 2013 +0100 ; readme, added samgrep ; https://github.com/lindenb/jvarkit/commit/888cbb7cd4e41628858e7aaf5e0ac7979dbd6e55
Tue Jun 4 15:20:17 2013 +0200 ; sam2tsv ; https://github.com/lindenb/jvarkit/commit/e81d4706dd51297677ddb64dcc69aaa681eab4af
Mon May 6 21:49:31 2013 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/dfa89c3e08b7b6d1ba766dbdc6c7c4279f7b7a3d
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samgrep** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Examples


#### Example 1


```

java -jar  dist/samgrep.jar -R r001  -- samtools-0.1.18/examples/toy.sam 

@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
@PG     ID:0    PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep     VN:dac03b80e9fd88a15648b22550e57d10c9bed725     CL:-R r001 samtools-0.1.18/examples/toy.sam
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
r001    83      ref     37      30      9M      =       7       -39     CAGCGCCAT       *

```


#### Example 4


```

java -jar  dist/samgrep.jar -R r001 -- -n 1 samtools-0.1.18/examples/toy.sam 

@HD     VN:1.4  SO:unsorted
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
@PG     ID:0    PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep     VN:dac03b80e9fd88a15648b22550e57d10c9bed725     CL:-R r001 -n 1 samtools-0.1.18/examples/toy.sam
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112

```







### Examples

#### Example 1

```

java -jar  dist/samgrep.jar -R r001  -- samtools-0.1.18/examples/toy.sam 

@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep	VN:dac03b80e9fd88a15648b22550e57d10c9bed725	CL:-R r001 samtools-0.1.18/examples/toy.sam
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
r001	83	ref	37	30	9M	=	7	-39	CAGCGCCAT	*

```





#### Example 4



```

java -jar  dist/samgrep.jar -R r001 -- -n 1 samtools-0.1.18/examples/toy.sam 

@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.samgrep.SamGrep	VN:dac03b80e9fd88a15648b22550e57d10c9bed725	CL:-R r001 -n 1 samtools-0.1.18/examples/toy.sam
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112

```





