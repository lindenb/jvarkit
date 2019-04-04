# SamColorTag

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Add the UCSC 'YC' color tag in a BAM. See http://software.broadinstitute.org/software/igv/book/export/html/6 and http://genome.ucsc.edu/goldenPath/help/hgBamTrackHelp.html


## Usage

```
Usage: samcolortag [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -e, --expression
      javascript expression
    -f, --file
      javascript file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -E, --ignoreErrors
      Ignore javascript/color errors
      Default: false
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
 * metadata
 * javascript
 * igv
 * visualization


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samcolortag
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamColorTag.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamColorTag.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samcolortag** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## How it works

Fill the  UCSC 'YC' color tag in a BAM using a javascript expression using
the nashorn javascript engine. The YC is used by IGV or the UCSC to colorize the reads

The script engine inject 'header' ( a htsjdk.samtools.SAMFileHeader ) and 'record' (a htsjdk.samtools.SAMRecord).
The script should return 

* a java.awt.color
* a String for a named color ('blue', 'red'...)
* a hexa color #FFFFF
* a rgb color 'rgb(100,200,100)'
* null or empty string (no YC tag, the tag is cleared)

## Example

The script:

```
function rnd(n) {
    return  Math.floor(Math.random() *n) ;
	}

function randomColor()
	{
	switch(rnd(10))
		{
		case 1: return "beige";
		case 2: return "#FF00AA";
		default: return "rgb("+rnd(255)+","+rnd(255)+","+rnd(255)+")";
		}
	return true;
	}

randomColor();
```

usage:

```
$ java -jar dist/samcolortag.jar -f script.js -o in.bam  out.bam

[lindenb@kaamelot-master01 jvarkit-git]$ samtools view out.bam | head
rotavirus_1_317_5:0:0_7:0:0_2de	99	rotavirus	1	60	70M	=	248	317	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:35,107,5	MD:Z:33G4A3T14A7T4	RG:Z:S1	NM:i:5	AS:i:45	XS:i:0
rotavirus_1_535_4:0:0_4:0:0_1a6	163	rotavirus	1	60	70M	=	466	535	GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:66,172,98	MD:Z:8A18G13C6G21	RG:Z:S1	NM:i:4	AS:i:50	XS:i:0
rotavirus_1_543_5:0:0_11:0:0_390	163	rotavirus	1	60	70M	=	487	530	GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:245,245,220	MD:Z:18G1G1T22T11A12	RG:Z:S1	NM:i:5	AS:i:45	XS:i:0
rotavirus_1_578_3:0:0_7:0:0_7c	99	rotavirus	1	60	70M	=	509	578	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:67,135,235	MD:Z:43A2C5A17	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0
rotavirus_1_497_4:0:0_5:0:0_2f6	163	rotavirus	1	60	70M	=	432	497	GGCATTTAATGCTTAACAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTCTTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	YC:Z:82,150,221	MD:Z:3T10T0T48A5	RG:Z:S1	NM:i:4	AS:i:51	XS:i:0
```

Another example, using the read mapping quality

```
var c=record.getMappingQuality();
"rgb("+c+","+c+","+c+")";
```


## See also:
  * com.github.lindenb.jvarkit.util.swing.ColorUtils
  * See http://software.broadinstitute.org/software/igv/book/export/html/6
  * http://genome.ucsc.edu/goldenPath/help/hgBamTrackHelp.html

## Screenshot

<img src="https://pbs.twimg.com/media/C_eefreXUAAO-rX.jpg"/>

