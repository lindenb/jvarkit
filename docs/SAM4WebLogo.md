# SAM4WebLogo

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sequence logo for different alleles or generated from SAM/BAM 


## Usage

```
Usage: sam4weblogo [options] Files
  Options:
    -c, --clipped, --clip
      Use Clipped Bases
      Default: false
    -q, --fastq
      [20180812]print fastq-like format. Was : 
      https://github.com/lindenb/jvarkit/issues/109 
      Default: false
    -fqp, --fqp
      [20180813] fastq padding quality character
      Default: -
    -fqu, --fqu
      [20180813] fastq unknown quality character
      Default: !
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -r, --region, --interval
      Region to observe: chrom:start-end
    -o, --output
      Output file. Optional . Default: stdout
    -readFilter, --readFilter
      [20171201](moved to jexl)A JEXL Expression that will be used to filter 
      out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: 'Accept all' (Empty expression)
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * visualization
 * logo



## See also in Biostars

 * [https://www.biostars.org/p/73021](https://www.biostars.org/p/73021)
 * [https://www.biostars.org/p/368200](https://www.biostars.org/p/368200)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew sam4weblogo
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam4weblogo/SAM4WebLogo.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam4weblogo/SAM4WebLogo.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam4weblogo/SAM4WebLogoTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam4weblogo/SAM4WebLogoTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sam4weblogo** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

"Sequence logo ( http://weblogo.berkeley.edu/logo.cgi ) for different alleles or generated from SAM/BAM" http://www.biostars.org/p/73021

![ScreenShot](https://raw.github.com/lindenb/jvarkit/master/doc/sam2weblogo.png)


## Example

```bash
$ java -jar dist/sam4weblogo.jar -r seq1:80-110  sorted.bam  2> /dev/null | head -n 50
>B7_593:4:106:316:452/1
TGTTG--------------------------
>B7_593:4:106:316:452a/1
TGTTG--------------------------
>B7_593:4:106:316:452b/1
TGTTG--------------------------
>B7_589:8:113:968:19/2
TGGGG--------------------------
>B7_589:8:113:968:19a/2
TGGGG--------------------------
>B7_589:8:113:968:19b/2
TGGGG--------------------------
>EAS54_65:3:321:311:983/1
TGTGGG-------------------------
>EAS54_65:3:321:311:983a/1
TGTGGG-------------------------
>EAS54_65:3:321:311:983b/1
TGTGGG-------------------------
>B7_591:6:155:12:674/2
TGTGGGGG-----------------------
>B7_591:6:155:12:674a/2
TGTGGGGG-----------------------
>B7_591:6:155:12:674b/2
TGTGGGGG-----------------------
>EAS219_FC30151:7:51:1429:1043/2
TGTGGGGGGCGCCG-----------------
>EAS219_FC30151:7:51:1429:1043a/2
TGTGGGGGGCGCCG-----------------
>EAS219_FC30151:7:51:1429:1043b/2
TGTGGGGGGCGCCG-----------------
>B7_591:5:42:540:501/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410/1
TGGGGGGGGCGCAGT----------------
>B7_591:5:42:540:501a/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410a/1
TGGGGGGGGCGCAGT----------------
>B7_591:5:42:540:501b/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410b/1
TGGGGGGGGCGCAGT----------------
```

### fastq-like output

```
$ java -jar dist/sam4weblogo.jar -r 'RF01:100-130' src/test/resources/S1.bam -q -c

@RF01_44_622_1:0:0_1:0:0_3a/1
TATTCTTCCAATAG-----------------
+
22222222222222                 
@RF01_44_499_0:0:0_3:0:0_7b/2
TATTCTTCCAATAG-----------------
+
22222222222222                 
@RF01_67_565_0:0:0_2:0:0_67/2
TATTCTTCCAATAGTGAATTAGAGAATAGAT
+
2222222222222222222222222222222
@RF01_94_620_1:0:0_2:0:0_15/2
TATTCTTCCAATAGTGAATTAGAGAATAGAT
+
2222222222222222222222222222222
@RF01_102_665_1:0:0_1:0:0_71/1
--TTCTTCCAATAGTGAATTAGAGAATAGAT
+
  22222222222222222222222222222
@RF01_110_504_2:0:0_1:0:0_5d/2
----------ATAGTGAATTAGATAATAGAT
+
          222222222222222222222
@RF01_121_598_1:0:0_3:0:0_6e/2
---------------------GAGAATAGAT
+
                     2222222222
```


## See also

* https://www.biostars.org/p/103052/
* http://www.sciencedirect.com/science/article/pii/S1874778715300210

