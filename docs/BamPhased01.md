# BamPhased01

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extract Reads from a SAM/BAM file supporting at least two variants in a VCF file.


## Usage

```
Usage: bamphased01 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    --buffer-size
      When we're looking for variant in a lare VCF file, load the variants in 
      an interval of 'N' bases instead of doing a random access for each 
      variant. 
      Default: 1000
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore-discordant-rg
      In paired mode, ignore discordant read-groups RG-ID.
      Default: false
    --mapq0
      If set. Do not remove the reads failing the test but set the MAPQ to 0.
      Default: false
    --min-supporting
      Min number of variants that should be supported by one read.
      Default: 2
    -o, --out
      Output file. Optional . Default: stdout
    --paired
      Activate Paired-end mode. Variant can be supported by the read or/and is 
      mate. Input must be sorted on query name using for example 'samtools 
      collate'. 
      Default: false
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --tag
      Tag in metadata of read containing the variants positions. Ignored if 
      empty 
      Default: XP
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
  * -V, --vcf
      Indexed VCf file
    --version
      print version and exit

```


## Keywords

 * vcf
 * phased
 * genotypes
 * bam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bamphased01
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210218

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/phased/BamPhased01.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/phased/BamPhased01.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamphased01** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


for @isamtalves

Only diallelic SNV are supported.

## Example

```
$ java -jar dist/bamphased01.jar \
	-V src/test/resources/rotavirus_rf.vcf.gz \
	src/test/resources/S4.bam 


@HD	VN:1.6	SO:coordinate
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:S4	SM:S4	LB:L4	CN:Nantes
@PG	ID:0	CL:-V src/test/resources/rotavirus_rf.vcf.gz src/test/resources/S4.bam	VN:04c54fe	PN:bamphased01
@CO	bamphased01. compilation:20210218183501 githash:04c54fe htsjdk:2.23.0 date:20210218183539. cmd:-V src/test/resources/rotavirus_rf.vcf.gz src/test/resources/S4.bam
RF10_109_650_2:2:0_2:0:0_13	99	RF10	109	60	70M	=	581	542	CGATACTCGAGGATCCAGGGATGGCGTATTATCCTTATCTAGCAACTGTCCTAACAGTTTTGTTCAGGTT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:4	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:51	XS:i:0
RF10_128_592_1:2:0_1:0:0_27	163	RF10	128	60	70M	=	523	465	GATGGCGTATTATCCTTAAATAGCATCTGTCCTAACAGTTTTGTTCAGGTTGCACAAAGCATCTATTCCA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:3	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:55	XS:i:0
RF10_133_733_1:2:0_1:0:0_f	99	RF10	133	60	70M	=	664	601	CGTATTATCCTTATATAGCAACTGTCCTAACAGTTTTGTTCAGGTTGCACAAAGCATCTATTCCAACAAT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:3	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:55	XS:i:0
RF10_137_727_1:2:0_2:0:0_8	163	RF10	137	60	70M	=	658	591	TTATCCTTATATAGCATCTGTCCTAACAGTTTTGTTCAGGTTGCACAAAGCATCTATTGCAACAATGAAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:3	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:57	XS:i:0
RF10_138_562_1:2:0_1:0:0_1e	163	RF10	138	60	70M	=	493	425	TATCCTTATATAGCATCTGTCCTAACAGTTATGTTCAGGTTGCACAAAGCATCTATTCCAACAATGAAAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:3	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:58	XS:i:0
```


paired mode:

```
$ samtools sort -n -T xx src/test/resources/S2.bam | \
	java -jar dist/bamphased01.jar -V src/test/resources/rotavirus_rf.vcf.gz   --min-supporting 2 --paired --ignore-discordant-rg

@HD	VN:1.6	SO:queryname
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:S2	SM:S2	LB:L2	CN:Nantes
@PG	ID:0	PN:bamphased01	VN:5d263eb	CL:-V src/test/resources/rotavirus_rf.vcf.gz --min-supporting 2 --paired --ignore-discordant-rg
@CO	bamphased01. compilation:20210220193801 githash:5d263eb htsjdk:2.23.0 date:20210220193911. cmd:-V src/test/resources/rotavirus_rf.vcf.gz --min-supporting 2 --paired --ignore-discordant-rg
RF03_2142_2592_1:1:0_2:1:0_1c	99	RF03	2142	60	70M	=	2523	451	AACTATTTTATTGGAATTAAGTTCAAAAATATACCCTATGAATATGATGATAAAGTACCCCATCTTACAT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:2	XP:Z:RF03_2201_G_C;RF03_2573_A_G	AS:i:60	XS:i:0
RF03_2142_2592_1:1:0_2:1:0_1c	147	RF03	2523	60	70M	=	2142	-451	ATAAAAGGCGACACACTGTTAGATATGACTGAGTGAGCTAAAAACTTAACGCACTGGTCACATCGTGACC	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF03_2201_G_C;RF03_2573_A_G	AS:i:55	XS:i:0
RF05_813_1307_2:1:0_1:1:0_10	99	RF05	813	60	70M	=	1238	495	AAAGCGTGAACGCATATAGTTGATGCTAGAAATTATATTAGTATTATGAACTCATCGTATACTGAGAATT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF05_879_C_A;RF05_1297_T_G	AS:i:56	XS:i:0
RF05_813_1307_2:1:0_1:1:0_10	147	RF05	1238	60	70M	=	813	-495	CATAAATGGAACGGAGTATGTATTATTAGACTATGAAGTGAACTGGGAAGTGAGGGGACGAGTCAGGCAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:2	XP:Z:RF05_879_C_A;RF05_1297_T_G	AS:i:60	XS:i:0
RF05_829_1327_3:1:0_4:1:0_16	99	RF05	829	60	70M	=	1258	499	TCGTTGATGCTCGAAATTATATTAGTATTATGAAGTCATCGTATACTGAGAATTACAGTGTGTCACAAAG	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:4	XP:Z:RF05_879_C_A;RF05_1297_T_G	AS:i:53	XS:i:0
RF05_829_1327_3:1:0_4:1:0_16	147	RF05	1258	60	70M	=	829	-499	TATTATTAGACTATCAAGTGAACTGGGAAGTGAGGGGACGAGTCATGCATAACAGGGATGCGAAAGTACC	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:5	XP:Z:RF05_879_C_A;RF05_1297_T_G	AS:i:45	XS:i:0
RF05_843_1306_0:1:0_2:1:0_23	83	RF05	1237	60	70M	=	843	-464	ACATAAATCGAACCGAGTATGTATTATTAGACTATGAAGTGAACTGGGAAGTGAGGGGACGAGTCATGCA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF05_1297_T_G;RF05_879_C_A	AS:i:55	XS:i:0
RF05_843_1306_0:1:0_2:1:0_23	163	RF05	843	60	70M	=	1237	464	AATTATATTAGTATTATGAACTCATCGTATACTGAGAATTACAGTGTGTCACAAAGATGTAAATTGTTTA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:1	XP:Z:RF05_1297_T_G;RF05_879_C_A	AS:i:65	XS:i:0
RF05_857_1394_2:1:0_1:1:0_2f	99	RF05	857	60	70M	=	1325	538	TATGAACTCATCGTATACTGAGAATTACAGTGTGTCACAAAGATGTACATTGATTACTAAGTATAAATTT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF05_879_C_A;RF05_1339_A_C	AS:i:55	XS:i:0
RF05_857_1394_2:1:0_1:1:0_2f	147	RF05	1325	60	70M	=	857	-538	TCCAAGAATTTTGACTATGAATGATACAAAGAAGATACTGAGTGCAATGATATTTGACTGGTTTGACACA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:2	XP:Z:RF05_879_C_A;RF05_1339_A_C	AS:i:64	XS:i:0
RF05_863_1322_3:1:0_0:1:0_0	83	RF05	1253	60	70M	=	863	-460	GTATGTATTATTAGACTATGAAGTGAACTGGGAAGTGAGGGGACGAGTCATGCAAAACATGGATGGGAAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:1	XP:Z:RF05_1297_T_G;RF05_879_C_A	AS:i:65	XS:i:0
RF05_863_1322_3:1:0_0:1:0_0	163	RF05	863	60	70M	=	1253	460	CTCATCGTATACTGAGAATTACAGTTTGTCACAAAGATGTAAAATGTTTACTAAGTATAAATATGGGATT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:4	XP:Z:RF05_1297_T_G;RF05_879_C_A	AS:i:50	XS:i:0
RF05_874_1375_3:1:0_1:1:0_4d	99	RF05	874	60	70M	=	1306	502	CTGAGAATTACAGTGTGTCACAAAGATTTAAATTGTTTACTAAGTATACATTTGGGATTGTATCAAGATA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:4	XP:Z:RF05_879_C_A;RF05_1339_A_C	AS:i:54	XS:i:0
RF05_874_1375_3:1:0_1:1:0_4d	147	RF05	1306	60	70M	=	874	-502	AAAACATGGATGGGAAAGTACCAAGAATTTTGACTATGAATGATACAAAGAAGATACTCAGTGCAATGAT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:2	XP:Z:RF05_879_C_A;RF05_1339_A_C	AS:i:60	XS:i:0
RF05_916_1354_3:0:0_1:2:0_53	147	RF05	1285	60	70M	=	916	-439	AAGTGAGGGGAAGAGTCATGCAAAACATGGATGGGAAAGTACCAAGAATTTTGACTATGAATGATACAAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF05_1297_T_G;RF05_1339_A_C	AS:i:55	XS:i:0
```

