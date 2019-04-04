# Biostar76892

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

fix strand of two paired reads close but on the same strand. 


## Usage

```
Usage: biostar76892 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -d, --maxc
      distance beween two reads.A distance specified as a positive 
      integer.Comma are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb 
      Default: 30
    -osf, --osf
      only save pairs of reads which have been corrected by this program
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



## See also in Biostars

 * [https://www.biostars.org/p/76892](https://www.biostars.org/p/76892)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar76892
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar76892.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar76892.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar76892** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

Before fixing:
```bash
samtools view src.bam | grep "HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906"
HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906        177     3       1264832 37      101M    =       1264940 109     AGGTGGTGAAGCATGAGATGTAGGGAGAGCTGCTTTAAAACCCAGCACAAGGCTGGTT
GTACTGGCTCACACCTGTAATCCCAGGTCTTTGGGAGGCTGAG    """#""""#"""#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""#""""""#""##"""""""#"#"   X0:i:1  X1:i:0  BD:
Z:ABBACCABBABAAABABAABAACBBABAA@AAABAABBAABB@AAAAAAABA@ABAA@BAA@@BAAAAAAAA@@BAAAABA@ABAAABAACBBACBAABAA       MD:Z:0C0C0C1C0C0G1G0T1T0T0G0C0C0T0T0C0T0G0T0A0C0A2T2C0A0T0G0
C1T0G0T0G2G0T0T0T0G1G0T1T1C0C0A0A0G0T0G0C0G0A1T0G1G0C0T1T0T2C0G0T0G0T0G1C1C0A1G1G0T1C0C1T0T2C0A0     RG:Z:idp63088   XG:i:0  BI:Z:ABBAEDCCCBCBBABABAAAA@CBBAB@A@@A@AAAAAAA
BA@@AAA@A@BA@@BAA@A@@@@BAA@A@@@A@@AAAAABA@@BAAABBACBBACBBABAA       AM:i:37 NM:i:75 SM:i:37 XM:i:1  XO:i:0  MQ:i:37 XT:A:U
HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906        113     3       1264940 37      101M    =       1264832 -109    GGTATCTCCATGCTCGAAGCCCTGACCTACTGTATTGCCCCGAAAGTCTTCCCTGCTG
TGGCTGCATCTTTTCCACGTGGATAATCTTGGTTCATCTCTAG    """##"""""""""""""""""""""""""#"""""""""""""""""""""""""""#"""""""""""""""#""""""""""""#""#""##"#"##"   X0:i:1  X1:i:0  BD:
Z:BBAABBBBAAABBBCBAABCBA@BAAAAAAABAAAAACCCBABAAAAAAACBAAAAABABA@AA@AAABBAAAAACB@BBAAAAAAAABBBAABBBBAAAA       MD:Z:0T0T1G0C0A1G0T1C0C0A1G0T0G0C0A0T0G0T0G0T0G2A0T4G0C0A0A0
T0G0T0G0C1G0G0T0G1C0A0G0T0T0G0C0A4C1A0T0G0C0G0T0G2G0G1C0G0T0G0A1C0G0T0G1G0C2T2T0C0G0T0G0T0A0T1       RG:Z:idp63088   XG:i:0  BI:Z:BABADDCCBBBCBBCBAABCBA@AABAAA@AAAAA@BBBB
BAAA@AA@AABA@@A@@A@BA@@A@AA@AAAAAAABB@BAAAAAAAA@CBAAABBBBAAAA       AM:i:37 NM:i:74 SM:i:37 XM:i:0  XO:i:0  MQ:i:37 XT:A:U
```
Fixing:
```bash
 java -jar dist/biostar76892.jar ONLYSAVEFIXED=true \
 	IN=src.bam \
 	OUT=fix.bam  \
 	VALIDATION_STRINGENCY=LENIENT
```
result

```bash
samtools view fix.bam | grep "HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906"
HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906        163     3       1264832 37      101M    =       1264940 109     AGGTGGTGAAGCATGAGATGTAGGGAGAGCTGCTTTAAAACCCAGCACAAGGCTGGTT
GTACTGGCTCACACCTGTAATCCCAGGTCTTTGGGAGGCTGAG    """#""""#"""#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""#""""""#""##"""""""#"#"   X0:i:1  X1:i:0  BD:
Z:ABBACCABBABAAABABAABAACBBABAA@AAABAABBAABB@AAAAAAABA@ABAA@BAA@@BAAAAAAAA@@BAAAABA@ABAAABAACBBACBAABAA       MD:Z:0C0C0C1C0C0G1G0T1T0T0G0C0C0T0T0C0T0G0T0A0C0A2T2C0A0T0G0
C1T0G0T0G2G0T0T0T0G1G0T1T1C0C0A0A0G0T0G0C0G0A1T0G1G0C0T1T0T2C0G0T0G0T0G1C1C0A1G1G0T1C0C1T0T2C0A0     RG:Z:idp63088   XG:i:0  BI:Z:ABBAEDCCCBCBBABABAAAA@CBBAB@A@@A@AAAAAAA
BA@@AAA@A@BA@@BAA@A@@@@BAA@A@@@A@@AAAAABA@@BAAABBACBBACBBABAA       AM:i:37 NM:i:75 SM:i:37 XM:i:1  XO:i:0  MQ:i:37 XT:A:U  rv:i:1
HWI-1KL149:6:C0KTCACXX:5:2107:2283:35906        83      3       1264940 37      101M    =       1264832 -109    GGTATCTCCATGCTCGAAGCCCTGACCTACTGTATTGCCCCGAAAGTCTTCCCTGCTG
TGGCTGCATCTTTTCCACGTGGATAATCTTGGTTCATCTCTAG    """##"""""""""""""""""""""""""#"""""""""""""""""""""""""""#"""""""""""""""#""""""""""""#""#""##"#"##"   X0:i:1  X1:i:0  BD:
Z:BBAABBBBAAABBBCBAABCBA@BAAAAAAABAAAAACCCBABAAAAAAACBAAAAABABA@AA@AAABBAAAAACB@BBAAAAAAAABBBAABBBBAAAA       MD:Z:0T0T1G0C0A1G0T1C0C0A1G0T0G0C0A0T0G0T0G0T0G2A0T4G0C0A0A0
T0G0T0G0C1G0G0T0G1C0A0G0T0T0G0C0A4C1A0T0G0C0G0T0G2G0G1C0G0T0G0A1C0G0T0G1G0C2T2T0C0G0T0G0T0A0T1       RG:Z:idp63088   XG:i:0  BI:Z:BABADDCCBBBCBBCBAABCBA@AABAAA@AAAAA@BBBB
BAAA@AA@AABA@@A@@A@BA@@A@AA@AAAAAAABB@BAAAAAAAA@CBAAABBBBAAAA       AM:i:37 NM:i:74 SM:i:37 XM:i:0  XO:i:0  MQ:i:37 XT:A:U  rv:i:1
```

