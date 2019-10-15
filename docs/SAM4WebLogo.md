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
    --format, -F
      output format.
      Default: fasta
      Possible Values: [fasta, fastq, tabular]
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
      Region to observe. A source of intervals. The following suffixes are 
      recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise 
      it could be an empty string (no interval) or a list of plain interval 
      separated by '[ \t\n;,]'
      Default: (unspecified)
    --naming
      How to print a Read. A format 'a la C-printf'. %% :% , %n:read name, %s: 
      read bases, %q: read quals, %f : read flags,%m: mapq, %c: contig, %b: 
      start, %B: unclipped start, %e: end, %E: unclipped end,%I: read group 
      id, %N: sample name,%S: SAM String.
      Default: %n (%f) %N
    --no-insert
      Do not show insertions
      Default: false
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
    -R, --reference
      For Reading CRAM. Indexed fasta Reference file. This file must be 
      indexed with samtools faidx and with picard CreateSequenceDictionary
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

## History

On October 14th, I tried to implement the insertions. I haven't tested this feature in depth.

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
$ java -jar dist/sam4weblogo.jar -r 'RF01:100-130' src/test/resources/S1.bam --format fastq -c

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
## tabular output

```
$ java -jar /home/lindenb/src/jvarkit-git/dist/sam4weblogo.jar  -r 'chr4:15648762-15648862'    src/test/resources/retrocopy01.bwa.bam --format tabular -c
CTTGAACCCAGGAGGGGGAGGTGCCAGGGAGCCGAGATCATGCCACTGCACCCCAGCCTGGGCAACAAACCAAACCTCCATCC-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1111:24870:51570/1
CTTGAACCCAGGAGGCAGAGGTGGCAGTGAGCAGAAAACAAGCCACTGCACCCCAGCCGGGGAAAAAAAACAAGACCCCATCTaA-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1107:19532:10468/1
CTTGAACCCAGGAGGCAGAGGTTGCAGTGAGCCGAGATCATGCCACTGCACTCCAGCCTGGGCAACAAAGCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2121:23043:26536/2
ctggacCCCAGGAGGCAGAGGTTGCAGGTACCCGAAATAATGCCCCTGCCCCCCAGCCTGGGCAACAAAGCAAGACCCCATCT-A-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1212:10642:64966/2
-------------GGCAGAGGTTGCAGTGAGCCGAGATCATGCCACTGCACTCCAGCCTGGGCAACAAAGCAAGACTCCATCT-CaAAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1119:27661:8095/2
--------------------------------CGAGATCATGCCACTGCACTCCAGCCTGGGCAACAAAGCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2224:26189:58778/2
-------------------------------------------CACTGCACTCCAGCCTGGGCAACAAAGCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2103:11972:45611/1
ctggacccaaggaggcagaggttgaagggagcagaaacaaggccacggccccCCAGCCTGGGAAACAAAGCAAACCCCAATCT-A-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2121:29701:6108/1
---------------------------------------------------------CTGGGCAACAAAGCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1114:11992:43044/1
---------------------------------------------------------CTGGGCAACAAAGCAAGACTCCATCT-CaAAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2109:9019:23477/2
----------------------------------------------------------------ACAAAGCAAGACTCCATCT-CaAAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2110:15706:60765/2
---------------------------------------------------------------------GCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1218:27600:65793/1
-------------------------------------------------------------------------GACNCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1209:1326:20805/1
-----------------------------------------------------------------------------CCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2201:25895:8095/2
gggccgatggtgaaaagtggcagggcgcaacaaaaaatgcccgtcacaccaaccgggggcaaaaaaaaaaaaccccaaccaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1102:32319:57829/2
-----------------------------------------------------agcagggaaaacaaaccaaaacacaacata-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1106:29437:72033/1
--------------------------------------------------------------------------------------AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1109:27813:4104/1
------------------------------------------ccacgaccccaccacggaaaaaaaaaaaaaacacaaaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1112:17137:65951/1
tttaaaccaagacgccagacgtccaagatatacgactacccgccactgcaaaaccacctgtgcacaaacacacaaataaaact-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1117:10693:25042/2
------------------------------------------------------caccagggaaaaaaaaaaaaactcaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1118:24962:18502/1
ttgacaccagaagcaaaagttaaaagtcaaccgaaatcagcaaaataatctacaacgtaggaaaaagccaaagtctcattcta-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1119:15432:30791/1
tgaaagacggggggcagaggttgcaagagcccaaaaacatgccccggcacccaaacgtgggcaacaaagaaaaacccaatcaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1123:1590:53891/1
----------------------------------------------------------------------------------a-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1215:26585:67586/1
--------------------------------------------------------------------------------------AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1216:18142:49900/1
---------------------------------------atgccactgaacccaagccagggcaaaaaagaaaaaaccaatcc-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1219:25043:11048/2
-------------------------------------------------------------------aaaaaaacacaaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2113:14215:55086/1
aataaacccaaaagaaaaagtgcgcagtaggccaaaatcaagcaaatgaaatcaagccggggaaaaaaacaaagaccccacat-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2113:14418:23618/1
-------------------------------------------aaccgaaacccaacccagaaaaaaaaaaaaaacacaaaca-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2120:21420:20629/2
--------------------------------------aaagcacagacaaaacaaccagacaaacaaaaaaaaaacaaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2212:28564:14301/1
tacgaaccagggaggctcgggctgaggaagccgaacttcattcctcccaccgcaagactggcaaataaatcaagtcccaattt-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2216:20324:64685/1
--------------------------------------------accaccccccacccgggccacaaaaaaaaaaccaaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2216:23429:26607/1
---------------------------------------------------------------------------------------AAAAAAAAAAAGAGAA X2:1:H7T57CCXY:8:2124:1306:6249/1
----------------------------------------------------------------------------------------AAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1222:16579:28488/1
-----------------------------------------------------------------------------------------AAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1222:17513:27890/1
-------------------------------------------------------------------------------------------AAAAAAAAAAAA X2:1:H7T57CCXY:8:2109:14316:7902/1
-------------------------------------------------------------------------------------------------AAAAAA X2:1:H7T57CCXY:8:2108:14357:62329/1
--------------------------------------------------------------------------------------------------AAAAA X2:1:H7T57CCXY:8:1124:18223:49162/1
```

## See also

* https://www.biostars.org/p/103052/
* http://www.sciencedirect.com/science/article/pii/S1874778715300210

